from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import subprocess
import time

import numpy as np
import pandas as pd


MANIFEST = Path(
    "data/derived/current_wgs_26_manifest.csv"
)

BED = Path(
    "data/derived/coordinates/all_ta_hg38.sorted.bed"
)

CACHE_ROOT = Path(
    ".cache/current-wgs"
)

RESULT_ROOT = Path(
    "results/current-wgs"
)

EXPECTED_INTERVALS = 66362
EXPECTED_BROKEN = 5377
EXPECTED_NONBROKEN = 60985


def sha256(path):
    h = hashlib.sha256()

    with path.open("rb") as f:
        for chunk in iter(
            lambda: f.read(1024 * 1024),
            b"",
        ):
            h.update(chunk)

    return h.hexdigest()


def run(cmd, stdout=None):
    subprocess.run(
        [str(x) for x in cmd],
        stdout=stdout,
        check=True,
    )


def count_lines(path):
    with path.open("rb") as f:
        return sum(1 for _ in f)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Reconstruct broken-TA FPBM from current public "
            "CCLE WGS using chromosome-parallel remote access."
        )
    )

    parser.add_argument(
        "model_id",
        help="DepMap model ID, e.g. ACH-000971",
    )

    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Concurrent chromosome fetches (default: 4)",
    )

    parser.add_argument(
        "--force",
        action="store_true",
        help="Replace an existing completed result",
    )

    return parser.parse_args()


args = parse_args()

if args.workers < 1:
    raise RuntimeError("--workers must be at least 1")

model_id = args.model_id

manifest = pd.read_csv(MANIFEST)

match = manifest[
    manifest["ModelID"] == model_id
]

if len(match) != 1:
    raise RuntimeError(
        f"Expected exactly one manifest row for {model_id}; "
        f"found {len(match)}"
    )

model = match.iloc[0]

seq_id = model["SequencingID"]
bam_url = model["BAMUrl"]
bai_url = model["BAIUrl"]
published_fpbm = float(model["PublishedFPBM"])

result_dir = RESULT_ROOT / model_id
result_dir.mkdir(
    parents=True,
    exist_ok=True,
)

final_bedcov = (
    result_dir
    / f"{model_id}_all_ta.bedcov.tsv"
)

result_json = (
    result_dir
    / "result.json"
)

if result_json.exists() and not args.force:
    raise RuntimeError(
        f"Completed result already exists: {result_json}\n"
        "Use --force only if you intentionally want to replace it."
    )

work = CACHE_ROOT / model_id
bed_dir = work / "beds"
bam_dir = work / "bams"
cov_dir = work / "bedcov"

for directory in (
    work,
    bed_dir,
    bam_dir,
    cov_dir,
):
    directory.mkdir(
        parents=True,
        exist_ok=True,
    )

remote_bai = (
    work
    / f"{seq_id}.remote.bai"
)


# ------------------------------------------------------------
# Canonical TA BED
# ------------------------------------------------------------

bed = pd.read_csv(
    BED,
    sep="\t",
    header=None,
    names=[
        "chrom",
        "start",
        "end",
        "locus",
    ],
)

if len(bed) != EXPECTED_INTERVALS:
    raise RuntimeError(
        f"Expected {EXPECTED_INTERVALS} intervals, "
        f"found {len(bed)}"
    )

if bed["locus"].duplicated().any():
    raise RuntimeError(
        "Canonical BED contains duplicate locus identifiers"
    )

broken_mask = (
    bed["locus"].str.startswith("brokenTA_")
    | bed["locus"].str.startswith("B_")
)

nonbroken_mask = (
    bed["locus"].str.startswith("NB_")
)

if int(broken_mask.sum()) != EXPECTED_BROKEN:
    raise RuntimeError(
        "Unexpected broken-TA locus count"
    )

if int(nonbroken_mask.sum()) != EXPECTED_NONBROKEN:
    raise RuntimeError(
        "Unexpected non-broken TA locus count"
    )

bed_sha256 = sha256(BED)


# ------------------------------------------------------------
# Download remote BAM index once.
# ------------------------------------------------------------

if not remote_bai.exists():
    print("Downloading BAM index...")

    run(
        [
            "curl",
            "--fail",
            "--location",
            "--silent",
            "--show-error",
            "-o",
            remote_bai,
            bai_url,
        ]
    )

if remote_bai.stat().st_size == 0:
    raise RuntimeError(
        "Downloaded BAM index is empty"
    )


# ------------------------------------------------------------
# Build one job per chromosome/reference contig.
# ------------------------------------------------------------

jobs = []

for chrom, group in bed.groupby(
    "chrom",
    sort=False,
):
    safe = chrom.replace("/", "_")

    chrom_bed = (
        bed_dir / f"{safe}.bed"
    )

    chrom_bam = (
        bam_dir / f"{safe}.bam"
    )

    chrom_bai = Path(
        str(chrom_bam) + ".bai"
    )

    chrom_cov = (
        cov_dir / f"{safe}.bedcov.tsv"
    )

    group.to_csv(
        chrom_bed,
        sep="\t",
        header=False,
        index=False,
    )

    jobs.append(
        {
            "chrom": chrom,
            "loci": len(group),
            "bed": chrom_bed,
            "bam": chrom_bam,
            "bai": chrom_bai,
            "cov": chrom_cov,
        }
    )


def valid_coverage(job):
    path = job["cov"]

    if not path.exists():
        return False

    return count_lines(path) == job["loci"]


def process_job(job):
    chrom = job["chrom"]

    # Resume completed chromosome work.
    if valid_coverage(job):
        job["bam"].unlink(
            missing_ok=True
        )
        job["bai"].unlink(
            missing_ok=True
        )

        return {
            "chrom": chrom,
            "loci": job["loci"],
            "seconds": 0.0,
            "resumed": True,
        }

    # Remove partial outputs from interrupted work.
    job["cov"].unlink(
        missing_ok=True
    )
    job["bam"].unlink(
        missing_ok=True
    )
    job["bai"].unlink(
        missing_ok=True
    )

    started = time.perf_counter()

    run(
        [
            "samtools",
            "view",
            "-M",
            "-L",
            job["bed"],
            "-b",
            "-o",
            job["bam"],
            "-X",
            bam_url,
            remote_bai,
        ]
    )

    run(
        [
            "samtools",
            "quickcheck",
            "-v",
            job["bam"],
        ]
    )

    run(
        [
            "samtools",
            "index",
            job["bam"],
        ]
    )

    with job["cov"].open("w") as out:
        run(
            [
                "samtools",
                "bedcov",
                job["bed"],
                job["bam"],
            ],
            stdout=out,
        )

    if not valid_coverage(job):
        raise RuntimeError(
            f"{chrom}: coverage row count mismatch"
        )

    # Coverage is now durable; large temporary alignment
    # files are no longer needed.
    job["bam"].unlink(
        missing_ok=True
    )
    job["bai"].unlink(
        missing_ok=True
    )

    elapsed = (
        time.perf_counter()
        - started
    )

    return {
        "chrom": chrom,
        "loci": job["loci"],
        "seconds": elapsed,
        "resumed": False,
    }


# ------------------------------------------------------------
# Run chromosome jobs concurrently.
# ------------------------------------------------------------

print()
print("CURRENT-WGS BROKEN-TA RECONSTRUCTION")
print("=" * 72)
print("ModelID:      ", model_id)
print("Cell line:    ", model["CellLineName"])
print("SequencingID: ", seq_id)
print("Workers:      ", args.workers)
print("Intervals:    ", len(bed))
print("BED SHA256:   ", bed_sha256)
print()

overall_started = time.perf_counter()

completed = []

with ThreadPoolExecutor(
    max_workers=args.workers
) as executor:

    futures = {
        executor.submit(
            process_job,
            job,
        ): job
        for job in jobs
    }

    for future in as_completed(futures):
        result = future.result()
        completed.append(result)

        marker = (
            "RESUME"
            if result["resumed"]
            else "FETCH "
        )

        print(
            f"[{len(completed):>2}/{len(jobs)}] "
            f"{marker} "
            f"{result['chrom']:<24} "
            f"loci={result['loci']:>5} "
            f"time={result['seconds']/60:>6.2f} min",
            flush=True,
        )

wall_seconds = (
    time.perf_counter()
    - overall_started
)


# ------------------------------------------------------------
# Combine chromosome coverage into canonical locus order.
# ------------------------------------------------------------

columns = [
    "chrom",
    "start",
    "end",
    "locus",
    "bedcov",
]

parts = []

for job in jobs:
    if not valid_coverage(job):
        raise RuntimeError(
            f"Incomplete chromosome output: {job['chrom']}"
        )

    part = pd.read_csv(
        job["cov"],
        sep="\t",
        header=None,
        names=columns,
    )

    parts.append(part)

coverage = pd.concat(
    parts,
    ignore_index=True,
)

if len(coverage) != EXPECTED_INTERVALS:
    raise RuntimeError(
        f"Combined coverage has {len(coverage)} rows"
    )

if coverage["locus"].duplicated().any():
    raise RuntimeError(
        "Combined coverage contains duplicate loci"
    )

canonical_order = pd.Series(
    np.arange(len(bed)),
    index=bed["locus"],
)

coverage["_order"] = (
    coverage["locus"]
    .map(canonical_order)
)

if coverage["_order"].isna().any():
    raise RuntimeError(
        "Coverage contains loci not present in canonical BED"
    )

coverage = (
    coverage
    .sort_values("_order")
    .drop(columns="_order")
    .reset_index(drop=True)
)

expected_coordinates = (
    bed[
        [
            "chrom",
            "start",
            "end",
            "locus",
        ]
    ]
    .reset_index(drop=True)
)

actual_coordinates = (
    coverage[
        [
            "chrom",
            "start",
            "end",
            "locus",
        ]
    ]
    .reset_index(drop=True)
)

if not expected_coordinates.equals(
    actual_coordinates
):
    raise RuntimeError(
        "Combined coverage coordinates differ from canonical BED"
    )

coverage.to_csv(
    final_bedcov,
    sep="\t",
    header=False,
    index=False,
)

coverage_sha256 = sha256(
    final_bedcov
)


# ------------------------------------------------------------
# Reconstruct FPBM using the published normalization.
# ------------------------------------------------------------

coverage["length"] = (
    coverage["end"]
    - coverage["start"]
)

coverage["FPB"] = (
    coverage["bedcov"]
    / (
        coverage["length"]
        * 2.0
    )
)

total_fpb = coverage["FPB"].sum()

if not np.isfinite(total_fpb) or total_fpb <= 0:
    raise RuntimeError(
        f"Invalid total FPB: {total_fpb}"
    )

coverage["FPBM"] = (
    coverage["FPB"]
    / total_fpb
    * 1_000_000
)

broken = coverage[
    coverage["locus"].str.startswith("brokenTA_")
    | coverage["locus"].str.startswith("B_")
]

nonbroken = coverage[
    coverage["locus"].str.startswith("NB_")
]

if len(broken) != EXPECTED_BROKEN:
    raise RuntimeError(
        f"Expected {EXPECTED_BROKEN} broken loci, "
        f"found {len(broken)}"
    )

if len(nonbroken) != EXPECTED_NONBROKEN:
    raise RuntimeError(
        f"Expected {EXPECTED_NONBROKEN} non-broken loci, "
        f"found {len(nonbroken)}"
    )

broken_fpbm = float(
    broken["FPBM"].mean()
)

nonbroken_fpbm = float(
    nonbroken["FPBM"].mean()
)

ratio = (
    broken_fpbm
    / published_fpbm
)


# ------------------------------------------------------------
# Durable provenance/result record.
# ------------------------------------------------------------

result = {
    "schema_version": 1,
    "model_id": model_id,
    "cell_line_name": model["CellLineName"],
    "lineage": model["Lineage"],
    "sequencing_id": seq_id,
    "profile_id": model["ProfileID"],
    "wrn_gene_effect_26q3": float(
        model["WRNGeneEffect"]
    ),
    "published_fpbm": published_fpbm,
    "current_broken_ta_mean_fpbm": broken_fpbm,
    "current_nonbroken_ta_mean_fpbm": nonbroken_fpbm,
    "current_to_published_ratio": ratio,
    "interval_count": EXPECTED_INTERVALS,
    "broken_interval_count": EXPECTED_BROKEN,
    "nonbroken_interval_count": EXPECTED_NONBROKEN,
    "canonical_bed": str(BED),
    "canonical_bed_sha256": bed_sha256,
    "coverage_file": str(final_bedcov),
    "coverage_sha256": coverage_sha256,
    "workers": args.workers,
    "wall_seconds": wall_seconds,
    "method": (
        "samtools view -M chromosome-parallel remote BAM access; "
        "local samtools bedcov; published FPB/FPBM normalization"
    ),
}

result_json.write_text(
    json.dumps(
        result,
        indent=2,
        sort_keys=True,
    )
    + "\n"
)


# ------------------------------------------------------------
# Successful extraction: remove resumable cache.
# Final coverage remains locally under results/current-wgs/.
# ------------------------------------------------------------

shutil.rmtree(
    work,
    ignore_errors=True,
)


print()
print("RESULT")
print("=" * 72)
print(
    f"Current broken-TA mean FPBM:  {broken_fpbm:.4f}"
)
print(
    f"Current non-broken mean FPBM: {nonbroken_fpbm:.4f}"
)
print(
    f"Published FPBM:                {published_fpbm:.4f}"
)
print(
    f"Current / published ratio:     {ratio:.3f}"
)
print()
print(
    f"Wall time:                     {wall_seconds/60:.2f} min"
)
print(
    f"Coverage SHA256:               {coverage_sha256}"
)
print()
print(
    f"Coverage: {final_bedcov}"
)
print(
    f"Result:   {result_json}"
)
print()
print("PASS: current-WGS reconstruction complete.")
