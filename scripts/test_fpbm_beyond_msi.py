from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr


SOURCE = Path("data/source/WRN_2020_SourceData_Fig4.xlsx")
GENE_EFFECT = Path("data/source/CRISPRGeneEffect.csv")
SIGNATURES = Path("data/source/OmicsGlobalSignatures.csv")

WRN_COLUMN = "WRN (7486)"
N_PERM = 100_000
SEED = 42


published = pd.read_excel(
    SOURCE,
    sheet_name="4d",
    engine="openpyxl",
)

gene_effect = pd.read_csv(GENE_EFFECT)

id_column = (
    "ModelID"
    if "ModelID" in gene_effect.columns
    else gene_effect.columns[0]
)

current = gene_effect[
    [id_column, WRN_COLUMN]
].copy()

current.columns = ["ModelID", "WRN_26Q3"]
current["WRN_26Q3"] = pd.to_numeric(
    current["WRN_26Q3"],
    errors="coerce",
)

signatures = pd.read_csv(SIGNATURES)

signatures["MSIScore"] = pd.to_numeric(
    signatures["MSIScore"],
    errors="coerce",
)

if "IsDefaultEntryForModel" in signatures.columns:
    default = (
        signatures["IsDefaultEntryForModel"]
        .astype(str)
        .str.lower()
    )

    keep = default.isin(["true", "1"])

    if keep.any():
        signatures = signatures[keep]

signatures = (
    signatures[
        ["ModelID", "MSIScore"]
    ]
    .dropna()
    .groupby(
        "ModelID",
        as_index=False,
    )
    .median(numeric_only=True)
)

data = (
    published
    .rename(
        columns={"DepMap_ID": "ModelID"}
    )
    .merge(
        current,
        on="ModelID",
        how="inner",
    )
    .merge(
        signatures,
        on="ModelID",
        how="inner",
    )
)

data = data[
    data["CCLE_MSI"] == "MSI"
].dropna(
    subset=[
        "FPBM_mean",
        "WRN_26Q3",
        "MSIScore",
    ]
).copy()

assert len(data) == 26


print("BROKEN-TA PHENOTYPE BEYOND OFFICIAL MSI")
print("=" * 60)
print("MSI models:", len(data))
print()


# Raw FPBM relationship
r, p = pearsonr(
    data["FPBM_mean"],
    data["WRN_26Q3"],
)

rs, ps = spearmanr(
    data["FPBM_mean"],
    data["WRN_26Q3"],
)

print("Raw FPBM vs 26Q3 WRN")
print(f"Pearson:   {r:.3f}   p={p:.4g}")
print(f"Spearman:  {rs:.3f}   p={ps:.4g}")
print()


# Official MSI relationship
r_msi, p_msi = pearsonr(
    data["MSIScore"],
    data["WRN_26Q3"],
)

print("Official MSI score vs 26Q3 WRN")
print(f"Pearson:   {r_msi:.3f}   p={p_msi:.4g}")
print()


# Residualize both FPBM and WRN against MSI score.
msi = data["MSIScore"].to_numpy(float)
fpbm = data["FPBM_mean"].to_numpy(float)
wrn = data["WRN_26Q3"].to_numpy(float)

design = np.column_stack(
    [
        np.ones(len(data)),
        msi,
    ]
)

fpbm_coef = np.linalg.lstsq(
    design,
    fpbm,
    rcond=None,
)[0]

wrn_coef = np.linalg.lstsq(
    design,
    wrn,
    rcond=None,
)[0]

fpbm_resid = fpbm - design @ fpbm_coef
wrn_resid = wrn - design @ wrn_coef

partial_r, partial_p = pearsonr(
    fpbm_resid,
    wrn_resid,
)

rank_r, rank_p = spearmanr(
    fpbm_resid,
    wrn_resid,
)


# Permutation test
rng = np.random.default_rng(SEED)

observed = abs(partial_r)
extreme = 0

for _ in range(N_PERM):
    shuffled = rng.permutation(fpbm_resid)

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


print("After removing official MSI score from BOTH")
print("=" * 60)
print(
    f"Partial Pearson:   "
    f"{partial_r:.3f}   p={partial_p:.4g}"
)
print(
    f"Residual Spearman: "
    f"{rank_r:.3f}   p={rank_p:.4g}"
)
print(
    f"Permutation p:     "
    f"{permutation_p:.5f}"
)

print()
print("Interpretation:")
print("Positive residual correlation means:")
print("higher broken-TA coverage than MSI predicts")
print("        ↓")
print("less WRN dependency than MSI predicts")
print()
print("Equivalently:")
print("more broken-TA coverage loss than MSI predicts")
print("        ↓")
print("more WRN dependency than MSI predicts")
