from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
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

GROUND_TRUTH = Path(
    "results/km12/KM12_all_ta.bedcov.tsv"
)

WORK = Path(
    "results/km12/chrom-parallel-benchmark"
)

BED_DIR = WORK / "beds"
BAM_DIR = WORK / "bams"
COV_DIR = WORK / "bedcov"

LOCAL_REMOTE_BAI = WORK / "CDS-vHArJs.remote.bai"

COMBINED = (
    WORK
    / "KM12_all_ta.chrom_parallel.bedcov.tsv"
)

MODEL_ID = "ACH-000969"
WORKERS = 4


def run(cmd, stdout=None):
    subprocess.run(
        [str(x) for x in cmd],
        stdout=stdout,
        check=True,
    )


manifest = pd.read_csv(MANIFEST)

row = manifest[
    manifest["ModelID"] == MODEL_ID
]

if len(row) != 1:
    raise RuntimeError(
        f"Expected one manifest row for {MODEL_ID}"
    )

row = row.iloc[0]

bam_url = row["BAMUrl"]
bai_url = row["BAIUrl"]

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

if len(bed) != 66362:
    raise RuntimeError(
        f"Expected 66362 loci, found {len(bed)}"
    )

for directory in (
    WORK,
    BED_DIR,
    BAM_DIR,
    COV_DIR,
):
    directory.mkdir(
        parents=True,
        exist_ok=True,
    )


print("KM12 CHROMOSOME-PARALLEL REMOTE EXTRACTION")
print("=" * 72)
print("Model:", MODEL_ID)
print("SequencingID:", row["SequencingID"])
print("Intervals:", len(bed))
print("Workers:", WORKERS)
print()


# ------------------------------------------------------------
# Download the remote BAM index once.
#
# All four workers share this local BAI so they do not each
# retrieve/cache another copy.
# ------------------------------------------------------------

if not LOCAL_REMOTE_BAI.exists():
    print("Downloading BAM index once...")

    run(
        [
            "curl",
            "--fail",
            "--location",
            "--silent",
            "--show-error",
            "-o",
            LOCAL_REMOTE_BAI,
            bai_url,
        ]
    )

print(
    "Local remote-index copy: "
    f"{LOCAL_REMOTE_BAI.stat().st_size / 1024**2:.1f} MiB"
)

print()


# ------------------------------------------------------------
# Create one canonical BED per chromosome.
# ------------------------------------------------------------

chrom_jobs = []

for chrom, group in bed.groupby(
    "chrom",
    sort=False,
):
    safe = chrom.replace("/", "_")

    chrom_bed = (
        BED_DIR / f"{safe}.bed"
    )

    chrom_bam = (
        BAM_DIR / f"{safe}.bam"
    )

    chrom_bedcov = (
        COV_DIR / f"{safe}.bedcov.tsv"
    )

    group.to_csv(
        chrom_bed,
        sep="\t",
        header=False,
        index=False,
    )

    chrom_jobs.append(
        (
            chrom,
            len(group),
            chrom_bed,
            chrom_bam,
            chrom_bedcov,
        )
    )


print(
    f"Chromosome jobs: {len(chrom_jobs)}"
)
print()


# ------------------------------------------------------------
# Fetch one chromosome job.
# ------------------------------------------------------------

def process_chromosome(job):
    (
        chrom,
        loci,
        chrom_bed,
        chrom_bam,
        chrom_bedcov,
    ) = job

    started = time.perf_counter()

    run(
        [
            "samtools",
            "view",
            "-M",
            "-L",
            chrom_bed,
            "-b",
            "-o",
            chrom_bam,
            "-X",
            bam_url,
            LOCAL_REMOTE_BAI,
        ]
    )

    run(
        [
            "samtools",
            "quickcheck",
            "-v",
            chrom_bam,
        ]
    )

    run(
        [
            "samtools",
            "index",
            chrom_bam,
        ]
    )

    with chrom_bedcov.open("w") as out:
        run(
            [
                "samtools",
                "bedcov",
                chrom_bed,
                chrom_bam,
            ],
            stdout=out,
        )

    elapsed = (
        time.perf_counter()
        - started
    )

    size_mb = (
        chrom_bam.stat().st_size
        / 1024**2
    )

    return {
        "chrom": chrom,
        "loci": loci,
        "seconds": elapsed,
        "bam_mb": size_mb,
        "bedcov": chrom_bedcov,
    }


# ------------------------------------------------------------
# Four concurrent remote readers.
# ------------------------------------------------------------

overall_started = time.perf_counter()

results = []

with ThreadPoolExecutor(
    max_workers=WORKERS
) as executor:

    futures = {
        executor.submit(
            process_chromosome,
            job,
        ): job[0]
        for job in chrom_jobs
    }

    for future in as_completed(futures):
        result = future.result()
        results.append(result)

        completed = len(results)

        print(
            f"[{completed:>2}/{len(chrom_jobs)}] "
            f"{result['chrom']:<6} "
            f"loci={result['loci']:>5}  "
            f"BAM={result['bam_mb']:>6.1f} MiB  "
            f"time={result['seconds']/60:>6.2f} min"
        )


elapsed = (
    time.perf_counter()
    - overall_started
)


# ------------------------------------------------------------
# Combine chromosome coverage results in canonical BED order.
# ------------------------------------------------------------

columns = [
    "chrom",
    "start",
    "end",
    "locus",
    "bedcov",
]

parts = []

for result in results:
    part = pd.read_csv(
        result["bedcov"],
        sep="\t",
        header=None,
        names=columns,
    )

    parts.append(part)

new = pd.concat(
    parts,
    ignore_index=True,
)

if len(new) != 66362:
    raise RuntimeError(
        f"Combined coverage has {len(new)} rows"
    )

order = pd.Series(
    np.arange(len(bed)),
    index=bed["locus"],
)

new["_order"] = (
    new["locus"]
    .map(order)
)

if new["_order"].isna().any():
    raise RuntimeError(
        "Combined result contains unknown loci"
    )

new = (
    new
    .sort_values("_order")
    .drop(columns="_order")
)

new.to_csv(
    COMBINED,
    sep="\t",
    header=False,
    index=False,
)


# ------------------------------------------------------------
# Exact regression against previous KM12 ground truth.
# ------------------------------------------------------------

old = pd.read_csv(
    GROUND_TRUTH,
    sep="\t",
    header=None,
    names=columns,
)

old = (
    old
    .set_index("locus")
    .sort_index()
)

check = (
    new
    .set_index("locus")
    .sort_index()
)

coordinates_match = old[
    ["chrom", "start", "end"]
].equals(
    check[
        ["chrom", "start", "end"]
    ]
)

coverage_difference = (
    old["bedcov"]
    - check["bedcov"]
).abs()

coverage_match = bool(
    (coverage_difference == 0).all()
)

max_difference = (
    coverage_difference.max()
)

total_bam_mb = sum(
    r["bam_mb"]
    for r in results
)


print()
print("REGRESSION RESULT")
print("=" * 72)
print("Rows:", len(check))
print(
    "Coordinates identical:",
    coordinates_match,
)
print(
    "Coverage identical:",
    coverage_match,
)
print(
    "Maximum absolute coverage difference:",
    max_difference,
)

print()
print("PERFORMANCE")
print("=" * 72)
print(
    f"Workers:          {WORKERS}"
)
print(
    f"Wall time:        {elapsed/60:.2f} min"
)
print(
    f"Total subset BAM: {total_bam_mb:.1f} MiB"
)
print(
    "Previous serial multi-region benchmark: "
    "112.00 min"
)

speedup = (
    112.00
    / (elapsed / 60)
)

print(
    f"Wall-time speedup: {speedup:.2f}x"
)


if not (
    coordinates_match
    and coverage_match
):
    raise RuntimeError(
        "Chromosome-parallel extraction "
        "does not reproduce KM12 ground truth"
    )


print()
print(
    "PASS: chromosome-parallel extraction "
    "exactly reproduces KM12."
)
