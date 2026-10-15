from pathlib import Path
import json

import numpy as np
import pandas as pd
from scipy.stats import linregress, pearsonr, spearmanr


COHORT = Path("results/current-wgs/cohort_summary.csv")
MANIFEST = Path("data/derived/current_wgs_26_manifest.csv")
INVENTORY = Path("data/source/public_wgs_inventory.csv")
PROFILES = Path("data/source/OmicsProfiles.csv")
COORD_PROV = Path(
    "data/derived/coordinates/build_provenance.json"
)

RESULT_ROOT = Path("results/current-wgs")

OUTPUT = Path(
    "results/current-wgs/calibration_audit.txt"
)

DIAGNOSTICS = Path(
    "results/current-wgs/calibration_audit.csv"
)

JHUEM = "ACH-000928"
ROBUST_Z_SCALE = 0.6744897501960817


def correlations(x, y):
    pr, pp = pearsonr(x, y)
    sr, sp = spearmanr(x, y)

    return pr, pp, sr, sp


def describe(label, df):
    published = df["PublishedFPBM"].to_numpy(float)
    current = df["CurrentBrokenTAFPBM"].to_numpy(float)
    ratio = df["CurrentToPublishedRatio"].to_numpy(float)

    fit = linregress(
        published,
        current,
    )

    through_origin = (
        np.dot(published, current)
        / np.dot(published, published)
    )

    pr, pp, sr, sp = correlations(
        published,
        current,
    )

    ratio_pr, ratio_pp, ratio_sr, ratio_sp = correlations(
        published,
        ratio,
    )

    median = float(np.median(ratio))
    mad = float(
        np.median(
            np.abs(ratio - median)
        )
    )

    relative_to_median = (
        ratio / median - 1.0
    )

    return {
        "label": label,
        "n": len(df),
        "pearson": pr,
        "pearson_p": pp,
        "spearman": sr,
        "spearman_p": sp,
        "intercept": fit.intercept,
        "slope": fit.slope,
        "r2": fit.rvalue ** 2,
        "through_origin_scale": through_origin,
        "ratio_mean": float(np.mean(ratio)),
        "ratio_median": median,
        "ratio_sd": float(np.std(ratio, ddof=1)),
        "ratio_min": float(np.min(ratio)),
        "ratio_max": float(np.max(ratio)),
        "ratio_mad": mad,
        "within_2pct_median": int(
            np.sum(
                np.abs(relative_to_median)
                <= 0.02
            )
        ),
        "within_3pct_median": int(
            np.sum(
                np.abs(relative_to_median)
                <= 0.03
            )
        ),
        "ratio_vs_published_pearson":
            ratio_pr,
        "ratio_vs_published_pearson_p":
            ratio_pp,
        "ratio_vs_published_spearman":
            ratio_sr,
        "ratio_vs_published_spearman_p":
            ratio_sp,
    }


cohort = pd.read_csv(COHORT)
manifest = pd.read_csv(MANIFEST)

if len(cohort) != 26:
    raise RuntimeError(
        f"Expected 26 cohort rows, found {len(cohort)}"
    )

if not cohort["Status"].eq("complete").all():
    raise RuntimeError(
        "Cohort contains incomplete models"
    )

if len(manifest) != 26:
    raise RuntimeError(
        f"Expected 26 manifest rows, found {len(manifest)}"
    )


# ------------------------------------------------------------
# Calibration statistics
# ------------------------------------------------------------

data = cohort.merge(
    manifest[
        [
            "ModelID",
            "SequencingID",
            "ProfileID",
            "BAMSizeGiB",
        ]
    ],
    on="ModelID",
    how="inner",
    validate="one_to_one",
)

if len(data) != 26:
    raise RuntimeError(
        "Cohort/manifest merge did not produce 26 rows"
    )

ratio_median = float(
    data["CurrentToPublishedRatio"].median()
)

ratio_mad = float(
    (
        data["CurrentToPublishedRatio"]
        - ratio_median
    )
    .abs()
    .median()
)

if ratio_mad <= 0:
    raise RuntimeError(
        "Ratio MAD is zero; robust outlier score undefined"
    )

data["RatioRobustZ"] = (
    ROBUST_Z_SCALE
    * (
        data["CurrentToPublishedRatio"]
        - ratio_median
    )
    / ratio_mad
)

data["PercentFromMedianScale"] = (
    (
        data["CurrentToPublishedRatio"]
        / ratio_median
    )
    - 1.0
) * 100.0

full = describe(
    "ALL_26",
    data,
)

without_jhuem = describe(
    "WITHOUT_JHUEM1",
    data[
        data["ModelID"] != JHUEM
    ].reset_index(drop=True),
)


# ------------------------------------------------------------
# Crude BAM-size relationship
# ------------------------------------------------------------

bam_pr, bam_pp = pearsonr(
    data["BAMSizeGiB"],
    data["CurrentToPublishedRatio"],
)

bam_sr, bam_sp = spearmanr(
    data["BAMSizeGiB"],
    data["CurrentToPublishedRatio"],
)

no_j = data[
    data["ModelID"] != JHUEM
]

bam_pr_no_j, bam_pp_no_j = pearsonr(
    no_j["BAMSizeGiB"],
    no_j["CurrentToPublishedRatio"],
)

bam_sr_no_j, bam_sp_no_j = spearmanr(
    no_j["BAMSizeGiB"],
    no_j["CurrentToPublishedRatio"],
)


# ------------------------------------------------------------
# JHUEM-1 lineage context
# ------------------------------------------------------------

jhuem = data[
    data["ModelID"] == JHUEM
]

if len(jhuem) != 1:
    raise RuntimeError(
        "Expected exactly one JHUEM-1 row"
    )

jhuem = jhuem.iloc[0]

uterus_peers = data[
    (data["Lineage"] == jhuem["Lineage"])
    & (data["ModelID"] != JHUEM)
]

uterus_peer_median = float(
    uterus_peers[
        "CurrentToPublishedRatio"
    ].median()
)

uterus_peer_min = float(
    uterus_peers[
        "CurrentToPublishedRatio"
    ].min()
)

uterus_peer_max = float(
    uterus_peers[
        "CurrentToPublishedRatio"
    ].max()
)


# ------------------------------------------------------------
# Public WGS candidate audit
# ------------------------------------------------------------

inventory = pd.read_csv(INVENTORY)

wgs = inventory[
    inventory["datatype"]
    .astype(str)
    .str.lower()
    .eq("wgs")
    &
    inventory["reference_genome"]
    .astype(str)
    .str.lower()
    .eq("hg38")
    &
    inventory["model_id"].isin(
        data["ModelID"]
    )
].copy()

wgs = wgs.dropna(
    subset=[
        "sample_id",
        "cram_bam_url",
        "crai_bai_url",
    ]
)

candidate_counts = (
    wgs.groupby("model_id")
    .size()
    .rename("PublicHg38WGSCandidates")
)

data = data.merge(
    candidate_counts,
    left_on="ModelID",
    right_index=True,
    how="left",
)

data["PublicHg38WGSCandidates"] = (
    data["PublicHg38WGSCandidates"]
    .fillna(0)
    .astype(int)
)

jhuem_candidates = wgs[
    wgs["model_id"] == JHUEM
].copy()


# ------------------------------------------------------------
# Profile-default audit
# ------------------------------------------------------------

profiles = pd.read_csv(PROFILES)

profile_wgs = profiles[
    profiles["DataType"]
    .astype(str)
    .str.lower()
    .eq("wgs")
].copy()

selected_profile = profile_wgs[
    profile_wgs["SequencingID"].isin(
        data["SequencingID"]
    )
][
    [
        "SequencingID",
        "IsDefaultEntryForModel",
        "IsDefaultEntryForMC",
    ]
].drop_duplicates()

data = data.merge(
    selected_profile,
    on="SequencingID",
    how="left",
)


# ------------------------------------------------------------
# Result-provenance consistency
# ------------------------------------------------------------

result_meta = []

for model_id in data["ModelID"]:
    path = (
        RESULT_ROOT
        / model_id
        / "result.json"
    )

    result = json.loads(
        path.read_text()
    )

    result_meta.append(
        {
            "ModelID": model_id,
            "ResultJSONSchema":
                result.get(
                    "schema_version"
                ),
            "CanonicalBEDSHA256":
                result.get(
                    "canonical_bed_sha256"
                ),
            "ExtractorSHA256":
                result.get(
                    "extractor_sha256"
                ),
            "FetchRetriesResult":
                result.get(
                    "fetch_retry_count"
                ),
            "ResumedJobsResult":
                result.get(
                    "resumed_job_count"
                ),
        }
    )

result_meta = pd.DataFrame(
    result_meta
)

data = data.merge(
    result_meta,
    on="ModelID",
    how="left",
    validate="one_to_one",
)

known_bed_shas = sorted(
    data["CanonicalBEDSHA256"]
    .dropna()
    .unique()
)

known_extractor_shas = sorted(
    data["ExtractorSHA256"]
    .dropna()
    .unique()
)


# ------------------------------------------------------------
# Coordinate-retention context
# ------------------------------------------------------------

coord = json.loads(
    COORD_PROV.read_text()
)

counts = coord["counts"]

source_loci = (
    counts["broken_hg19"]
    + counts["nonbroken_hg19"]
)

mapped_loci = (
    counts["broken_hg38"]
    + counts["nonbroken_hg38"]
)

retention = (
    mapped_loci
    / source_loci
)


# ------------------------------------------------------------
# Diagnostics output
# ------------------------------------------------------------

diagnostics = data[
    [
        "ModelID",
        "CellLineName",
        "Lineage",
        "PublishedFPBM",
        "CurrentBrokenTAFPBM",
        "CurrentToPublishedRatio",
        "RatioRobustZ",
        "PercentFromMedianScale",
        "BAMSizeGiB",
        "SequencingID",
        "ProfileID",
        "PublicHg38WGSCandidates",
        "IsDefaultEntryForModel",
        "IsDefaultEntryForMC",
        "FetchRetries",
        "ResumedJobs",
        "ResultJSONSchema",
        "CanonicalBEDSHA256",
        "ExtractorSHA256",
    ]
].sort_values(
    "RatioRobustZ",
    key=lambda x: x.abs(),
    ascending=False,
)

diagnostics.to_csv(
    DIAGNOSTICS,
    index=False,
)


# ------------------------------------------------------------
# Human-readable report
# ------------------------------------------------------------

lines = []


def emit(text=""):
    lines.append(text)
    print(text)


def report(stats):
    emit(stats["label"])
    emit("-" * 72)
    emit(
        f"Models: {stats['n']}"
    )
    emit(
        f"Published/current Pearson:  "
        f"{stats['pearson']:.6f}"
    )
    emit(
        f"Published/current Spearman: "
        f"{stats['spearman']:.6f}"
    )
    emit(
        f"Linear fit: Current = "
        f"{stats['intercept']:.6f} + "
        f"{stats['slope']:.6f} * Published"
    )
    emit(
        f"R^2: {stats['r2']:.6f}"
    )
    emit(
        f"Through-origin scale: "
        f"{stats['through_origin_scale']:.6f}"
    )
    emit(
        f"Ratio mean:   "
        f"{stats['ratio_mean']:.6f}"
    )
    emit(
        f"Ratio median: "
        f"{stats['ratio_median']:.6f}"
    )
    emit(
        f"Ratio SD:     "
        f"{stats['ratio_sd']:.6f}"
    )
    emit(
        f"Ratio MAD:    "
        f"{stats['ratio_mad']:.6f}"
    )
    emit(
        f"Ratio range:  "
        f"{stats['ratio_min']:.6f} to "
        f"{stats['ratio_max']:.6f}"
    )
    emit(
        "Within +/-2% of median scale: "
        f"{stats['within_2pct_median']}"
        f"/{stats['n']}"
    )
    emit(
        "Within +/-3% of median scale: "
        f"{stats['within_3pct_median']}"
        f"/{stats['n']}"
    )
    emit(
        "Ratio vs published FPBM Pearson: "
        f"{stats['ratio_vs_published_pearson']:.6f}"
    )
    emit(
        "Ratio vs published FPBM Spearman: "
        f"{stats['ratio_vs_published_spearman']:.6f}"
    )
    emit()


emit("CURRENT-WGS FPBM CALIBRATION AUDIT")
emit("=" * 72)
emit()

report(full)
report(without_jhuem)

emit("JHUEM-1")
emit("-" * 72)
emit(
    f"Ratio: {jhuem['CurrentToPublishedRatio']:.6f}"
)
emit(
    f"Robust ratio z: "
    f"{jhuem['RatioRobustZ']:.6f}"
)
emit(
    f"Percent from cohort median scale: "
    f"{jhuem['PercentFromMedianScale']:+.2f}%"
)
emit(
    f"BAM size: {jhuem['BAMSizeGiB']:.2f} GiB"
)
emit(
    f"Fetch retries: {jhuem['FetchRetries']}"
)
emit(
    f"Resumed jobs: {jhuem['ResumedJobs']}"
)
emit(
    f"Public hg38 WGS candidates: "
    f"{len(jhuem_candidates)}"
)
emit(
    f"Selected SequencingID: "
    f"{jhuem['SequencingID']}"
)
emit(
    f"Selected ProfileID: "
    f"{jhuem['ProfileID']}"
)
emit()

emit("Other uterus models")
emit("-" * 72)
emit(
    f"Peers: {len(uterus_peers)}"
)
emit(
    f"Median ratio: {uterus_peer_median:.6f}"
)
emit(
    f"Range: {uterus_peer_min:.6f} "
    f"to {uterus_peer_max:.6f}"
)
emit()

emit("WGS candidate selection")
emit("-" * 72)
emit(
    "Models with one public hg38 WGS candidate: "
    f"{int((data['PublicHg38WGSCandidates'] == 1).sum())}/26"
)
emit(
    "Models with >1 public hg38 WGS candidate: "
    f"{int((data['PublicHg38WGSCandidates'] > 1).sum())}/26"
)
emit()

if len(jhuem_candidates):
    emit("JHUEM-1 candidate rows:")
    for row in jhuem_candidates.itertuples():
        emit(
            f"  {row.sample_id} "
            f"profile={row.profile_id} "
            f"size={int(row.cram_bam_size)/1024**3:.2f} GiB"
        )
    emit()

emit("BAM-size relationship to calibration ratio")
emit("-" * 72)
emit(
    f"All 26 Pearson:  {bam_pr:.6f} "
    f"(p={bam_pp:.6g})"
)
emit(
    f"All 26 Spearman: {bam_sr:.6f} "
    f"(p={bam_sp:.6g})"
)
emit(
    f"Without JHUEM Pearson:  "
    f"{bam_pr_no_j:.6f} "
    f"(p={bam_pp_no_j:.6g})"
)
emit(
    f"Without JHUEM Spearman: "
    f"{bam_sr_no_j:.6f} "
    f"(p={bam_sp_no_j:.6g})"
)
emit()

emit("Coordinate context")
emit("-" * 72)
emit(
    f"Original hg19 loci: {source_loci}"
)
emit(
    f"Mapped hg38 loci:   {mapped_loci}"
)
emit(
    f"Coordinate retention: "
    f"{retention*100:.4f}%"
)
emit(
    "Coordinate retention is descriptive only; "
    "it does not establish the cause of calibration shift."
)
emit()

emit("Result provenance")
emit("-" * 72)
emit(
    f"Known canonical BED SHA values: "
    f"{len(known_bed_shas)}"
)
for value in known_bed_shas:
    emit(f"  {value}")

emit(
    f"Known extractor SHA values: "
    f"{len(known_extractor_shas)}"
)
for value in known_extractor_shas:
    emit(f"  {value}")

emit()
emit(f"Diagnostics: {DIAGNOSTICS}")

OUTPUT.write_text(
    "\n".join(lines) + "\n"
)
