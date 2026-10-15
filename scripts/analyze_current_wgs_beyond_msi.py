from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr, t


COHORT = Path("results/current-wgs/cohort_summary.csv")
SIGNATURES = Path("data/source/OmicsGlobalSignatures.csv")

OUTPUT = Path(
    "results/current-wgs/beyond_msi_analysis.txt"
)

DIAGNOSTICS = Path(
    "results/current-wgs/beyond_msi_diagnostics.csv"
)

N_PERM = 100_000
SEED = 42


cohort = pd.read_csv(COHORT)

if len(cohort) != 26:
    raise RuntimeError(
        f"Expected 26 cohort models, found {len(cohort)}"
    )

if not cohort["Status"].eq("complete").all():
    raise RuntimeError(
        "Cohort contains incomplete models"
    )


signatures = pd.read_csv(SIGNATURES)

required_signature_columns = [
    "ModelID",
    "IsDefaultEntryForModel",
    "MSIScore",
]

missing = [
    column
    for column in required_signature_columns
    if column not in signatures.columns
]

if missing:
    raise RuntimeError(
        f"Missing signature columns: {missing}"
    )


# Canonical 26Q3 MSI-score selection:
# 1. model-default rows only
# 2. numeric, non-null MSIScore
# 3. exactly one score per cohort model
default_flag = (
    signatures["IsDefaultEntryForModel"]
    .astype(str)
    .str.strip()
    .str.lower()
)

signatures["MSIScore"] = pd.to_numeric(
    signatures["MSIScore"],
    errors="coerce",
)

msi = signatures[
    default_flag.eq("yes")
    & signatures["MSIScore"].notna()
    & signatures["ModelID"].isin(
        cohort["ModelID"]
    )
][
    [
        "ModelID",
        "MSIScore",
    ]
].copy()

if len(msi) != 26:
    raise RuntimeError(
        f"Expected 26 default non-null MSI rows, found {len(msi)}"
    )

if msi["ModelID"].nunique() != 26:
    duplicates = (
        msi.groupby("ModelID")
        .size()
        .loc[lambda x: x > 1]
    )

    raise RuntimeError(
        "Expected one MSI score per model. "
        f"Duplicates: {duplicates.to_dict()}"
    )


data = cohort.merge(
    msi,
    on="ModelID",
    how="inner",
    validate="one_to_one",
)

required_analysis_columns = [
    "CurrentBrokenTAFPBM",
    "PublishedFPBM",
    "WRNGeneEffect",
    "MSIScore",
]

data = data.dropna(
    subset=required_analysis_columns
).copy()

if len(data) != 26:
    raise RuntimeError(
        f"Expected 26 analyzable models, found {len(data)}"
    )


def correlation(x, y):
    pr, pp = pearsonr(x, y)
    sr, sp = spearmanr(x, y)

    return {
        "pearson_r": pr,
        "pearson_p": pp,
        "spearman_r": sr,
        "spearman_p": sp,
    }


def residualize(values, covariate):
    design = np.column_stack(
        [
            np.ones(len(values)),
            covariate,
        ]
    )

    coef = np.linalg.lstsq(
        design,
        values,
        rcond=None,
    )[0]

    return values - design @ coef


current_raw = correlation(
    data["CurrentBrokenTAFPBM"],
    data["WRNGeneEffect"],
)

published_raw = correlation(
    data["PublishedFPBM"],
    data["WRNGeneEffect"],
)

msi_wrn = correlation(
    data["MSIScore"],
    data["WRNGeneEffect"],
)

current_msi = correlation(
    data["CurrentBrokenTAFPBM"],
    data["MSIScore"],
)

published_msi = correlation(
    data["PublishedFPBM"],
    data["MSIScore"],
)


msi_values = data["MSIScore"].to_numpy(float)
current_values = data[
    "CurrentBrokenTAFPBM"
].to_numpy(float)
published_values = data[
    "PublishedFPBM"
].to_numpy(float)
wrn_values = data[
    "WRNGeneEffect"
].to_numpy(float)


current_resid = residualize(
    current_values,
    msi_values,
)

published_resid = residualize(
    published_values,
    msi_values,
)

wrn_resid = residualize(
    wrn_values,
    msi_values,
)


current_partial_r = pearsonr(
    current_resid,
    wrn_resid,
)[0]

partial_df = len(data) - 3

current_partial_t = (
    current_partial_r
    * np.sqrt(
        partial_df
        / (1 - current_partial_r ** 2)
    )
)

current_partial_p = (
    2
    * t.sf(
        abs(current_partial_t),
        df=partial_df,
    )
)

current_rank_r, current_rank_p = spearmanr(
    current_resid,
    wrn_resid,
)

published_partial_r = pearsonr(
    published_resid,
    wrn_resid,
)[0]

published_partial_t = (
    published_partial_r
    * np.sqrt(
        partial_df
        / (1 - published_partial_r ** 2)
    )
)

published_partial_p = (
    2
    * t.sf(
        abs(published_partial_t),
        df=partial_df,
    )
)

published_rank_r, published_rank_p = spearmanr(
    published_resid,
    wrn_resid,
)


# Two-sided permutation test for reconstructed current FPBM.
rng = np.random.default_rng(SEED)

observed = abs(current_partial_r)
extreme = 0

for _ in range(N_PERM):
    shuffled = rng.permutation(
        current_resid
    )

    permuted_r = np.corrcoef(
        shuffled,
        wrn_resid,
    )[0, 1]

    if abs(permuted_r) >= observed:
        extreme += 1

permutation_p = (
    extreme + 1
) / (
    N_PERM + 1
)


data["CurrentFPBMResidualAfterMSI"] = (
    current_resid
)

data["PublishedFPBMResidualAfterMSI"] = (
    published_resid
)

data["WRNResidualAfterMSI"] = (
    wrn_resid
)

data["ResidualProduct"] = (
    current_resid * wrn_resid
)

diagnostics = data[
    [
        "ModelID",
        "CellLineName",
        "Lineage",
        "MSIScore",
        "CurrentBrokenTAFPBM",
        "PublishedFPBM",
        "WRNGeneEffect",
        "CurrentFPBMResidualAfterMSI",
        "PublishedFPBMResidualAfterMSI",
        "WRNResidualAfterMSI",
        "ResidualProduct",
    ]
].copy()

diagnostics.to_csv(
    DIAGNOSTICS,
    index=False,
)


lines = []


def emit(text=""):
    lines.append(text)
    print(text)


emit("CURRENT-WGS FPBM BEYOND OFFICIAL MSI SCORE")
emit("=" * 72)
emit(f"Models: {len(data)}")
emit(f"Cohort: {COHORT}")
emit(f"MSI source: {SIGNATURES}")
emit()
emit("MSI row-selection rule")
emit("-" * 72)
emit('IsDefaultEntryForModel == "Yes"')
emit("MSIScore must be numeric and non-null")
emit("Exactly one selected MSI score per cohort model")
emit()
emit("Raw reconstructed FPBM vs current WRN")
emit("-" * 72)
emit(
    f"Pearson r:  {current_raw['pearson_r']:.6f}   "
    f"p={current_raw['pearson_p']:.6g}"
)
emit(
    f"Spearman r: {current_raw['spearman_r']:.6f}   "
    f"p={current_raw['spearman_p']:.6g}"
)
emit()
emit("Official MSI score vs current WRN")
emit("-" * 72)
emit(
    f"Pearson r:  {msi_wrn['pearson_r']:.6f}   "
    f"p={msi_wrn['pearson_p']:.6g}"
)
emit(
    f"Spearman r: {msi_wrn['spearman_r']:.6f}   "
    f"p={msi_wrn['spearman_p']:.6g}"
)
emit()
emit("Reconstructed FPBM vs official MSI score")
emit("-" * 72)
emit(
    f"Pearson r:  {current_msi['pearson_r']:.6f}   "
    f"p={current_msi['pearson_p']:.6g}"
)
emit(
    f"Spearman r: {current_msi['spearman_r']:.6f}   "
    f"p={current_msi['spearman_p']:.6g}"
)
emit()
emit("After removing MSI score from BOTH")
emit("-" * 72)
emit(
    f"Partial Pearson:   {current_partial_r:.6f}   "
    f"p={current_partial_p:.6g}"
)
emit(
    f"Partial Pearson df: {partial_df}"
)
emit(
    f"Residual Spearman: {current_rank_r:.6f}   "
    f"p={current_rank_p:.6g}"
)
emit(
    f"Permutation p:     {permutation_p:.6f}"
)
emit()
emit("Published-FPBM reference using same MSI rows")
emit("-" * 72)
emit(
    f"Raw Pearson:       {published_raw['pearson_r']:.6f}"
)
emit(
    f"Raw Spearman:      {published_raw['spearman_r']:.6f}"
)
emit(
    f"FPBM-MSI Pearson:  {published_msi['pearson_r']:.6f}"
)
emit(
    f"Partial Pearson:   {published_partial_r:.6f}   "
    f"p={published_partial_p:.6g}"
)
emit(
    f"Residual Spearman: {published_rank_r:.6f}   "
    f"p={published_rank_p:.6g}"
)
emit()
emit("Interpretation of residual direction")
emit("-" * 72)
emit(
    "Positive residual correlation means models with "
    "more broken-TA coverage than MSI predicts also tend "
    "to have less WRN dependency than MSI predicts."
)
emit(
    "Equivalently, more broken-TA coverage loss than MSI "
    "predicts tends to track with stronger WRN dependency."
)
emit()
emit(
    f"Permutation seed: {SEED}"
)
emit(
    f"Permutations: {N_PERM}"
)
emit(
    f"Diagnostics: {DIAGNOSTICS}"
)

OUTPUT.write_text(
    "\n".join(lines) + "\n"
)
