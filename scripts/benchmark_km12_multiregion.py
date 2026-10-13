from pathlib import Path
import hashlib
import subprocess
import time

import pandas as pd


MANIFEST = Path("data/derived/current_wgs_26_manifest.csv")
BED = Path("data/derived/coordinates/all_ta_hg38.sorted.bed")
GROUND_TRUTH = Path("results/km12/KM12_all_ta.bedcov.tsv")

WORK = Path("results/km12/multiregion-benchmark")
SUBSET_BAM = WORK / "KM12_ta_reads.bam"
SUBSET_BAI = WORK / "KM12_ta_reads.bam.bai"
NEW_BEDCOV = WORK / "KM12_all_ta.multiregion.bedcov.tsv"

MODEL_ID = "ACH-000969"


def run(cmd):
    print()
    print("$", " ".join(str(x) for x in cmd))
    started = time.perf_counter()

    subprocess.run(
        [str(x) for x in cmd],
        check=True,
    )

    elapsed = time.perf_counter() - started

    print(f"Elapsed: {elapsed / 60:.2f} minutes")

    return elapsed


def sha256(path):
    h = hashlib.sha256()

    with path.open("rb") as f:
        for chunk in iter(
            lambda: f.read(1024 * 1024),
            b"",
        ):
            h.update(chunk)

    return h.hexdigest()


manifest = pd.read_csv(MANIFEST)

row = manifest[
    manifest["ModelID"] == MODEL_ID
]

if len(row) != 1:
    raise RuntimeError(
        f"Expected one manifest row for {MODEL_ID}, found {len(row)}"
    )

row = row.iloc[0]

bam_url = row["BAMUrl"]
bai_url = row["BAIUrl"]

if not BED.exists():
    raise RuntimeError(f"Missing BED: {BED}")

if not GROUND_TRUTH.exists():
    raise RuntimeError(
        f"Missing KM12 ground truth: {GROUND_TRUTH}"
    )

bed_lines = sum(1 for _ in BED.open())

if bed_lines != 66362:
    raise RuntimeError(
        f"Expected 66362 intervals, found {bed_lines}"
    )

WORK.mkdir(
    parents=True,
    exist_ok=True,
)

print("KM12 MULTI-REGION REMOTE EXTRACTION BENCHMARK")
print("=" * 72)
print("Model:", MODEL_ID)
print("SequencingID:", row["SequencingID"])
print("Intervals:", bed_lines)
print("BED SHA256:", sha256(BED))
print()

# ------------------------------------------------------------
# Fetch only alignments overlapping the union of TA intervals.
#
# -M uses samtools' multi-region iterator.
# -L supplies the canonical sorted TA BED.
# -X supplies BAM and BAI explicitly because the remote index
# lives at a separate URL.
# ------------------------------------------------------------

fetch_seconds = run(
    [
        "samtools",
        "view",
        "-M",
        "-L",
        BED,
        "-b",
        "-o",
        SUBSET_BAM,
        "-X",
        bam_url,
        bai_url,
    ]
)

if not SUBSET_BAM.exists():
    raise RuntimeError(
        "samtools view did not create subset BAM"
    )

print()
print(
    "Subset BAM size: "
    f"{SUBSET_BAM.stat().st_size / 1024**2:.1f} MiB"
)

# Verify the subset BAM is healthy.
run(
    [
        "samtools",
        "quickcheck",
        "-v",
        SUBSET_BAM,
    ]
)

# The multi-region iterator should preserve coordinate order.
# Indexing is also a practical verification of that assumption.
index_seconds = run(
    [
        "samtools",
        "index",
        SUBSET_BAM,
    ]
)

if not SUBSET_BAI.exists():
    raise RuntimeError(
        "samtools index did not create local BAI"
    )

# ------------------------------------------------------------
# Run bedcov locally over the extracted reads.
# ------------------------------------------------------------

print()
print(
    "$ samtools bedcov "
    f"{BED} {SUBSET_BAM} > {NEW_BEDCOV}"
)

started = time.perf_counter()

with NEW_BEDCOV.open("w") as out:
    subprocess.run(
        [
            "samtools",
            "bedcov",
            BED,
            SUBSET_BAM,
        ],
        stdout=out,
        check=True,
    )

bedcov_seconds = time.perf_counter() - started

print(
    f"Elapsed: {bedcov_seconds / 60:.2f} minutes"
)

# ------------------------------------------------------------
# Exact numerical regression against the existing KM12 result.
# Ground truth may be in unsorted BED order, so compare by
# locus identifier rather than row number.
# ------------------------------------------------------------

columns = [
    "chrom",
    "start",
    "end",
    "locus",
    "bedcov",
]

old = pd.read_csv(
    GROUND_TRUTH,
    sep="\t",
    header=None,
    names=columns,
)

new = pd.read_csv(
    NEW_BEDCOV,
    sep="\t",
    header=None,
    names=columns,
)

if len(old) != 66362 or len(new) != 66362:
    raise RuntimeError(
        f"Unexpected row counts: old={len(old)}, new={len(new)}"
    )

old = old.set_index("locus").sort_index()
new = new.set_index("locus").sort_index()

if set(old.index) != set(new.index):
    raise RuntimeError(
        "Locus sets differ between benchmark and ground truth"
    )

coordinates_match = (
    old[
        ["chrom", "start", "end"]
    ].equals(
        new[
            ["chrom", "start", "end"]
        ]
    )
)

coverage_match = (
    old["bedcov"].to_numpy()
    == new["bedcov"].to_numpy()
).all()

max_abs_difference = (
    old["bedcov"] - new["bedcov"]
).abs().max()

print()
print("REGRESSION RESULT")
print("=" * 72)
print("Rows:", len(new))
print("Coordinates identical:", coordinates_match)
print("Coverage identical:", coverage_match)
print(
    "Maximum absolute coverage difference:",
    max_abs_difference,
)

print()
print("TIMING")
print("=" * 72)
print(
    f"Remote multi-region fetch: {fetch_seconds / 60:.2f} min"
)
print(
    f"Local BAM indexing:        {index_seconds / 60:.2f} min"
)
print(
    f"Local bedcov:              {bedcov_seconds / 60:.2f} min"
)
print(
    "Total:                     "
    f"{(fetch_seconds + index_seconds + bedcov_seconds) / 60:.2f} min"
)

if not coordinates_match or not coverage_match:
    raise RuntimeError(
        "Optimized extraction does not reproduce KM12 ground truth"
    )

print()
print("PASS: optimized extraction exactly reproduces KM12.")
