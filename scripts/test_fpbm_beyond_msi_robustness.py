from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr


SOURCE = Path("data/source/WRN_2020_SourceData_Fig4.xlsx")
GENE_EFFECT = Path("data/source/CRISPRGeneEffect.csv")
SIGNATURES = Path("data/source/OmicsGlobalSignatures.csv")

WRN_COLUMN = "WRN (7486)"
BOOTSTRAPS = 20_000
SEED = 42


def load_data():
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
        .rename(columns={"DepMap_ID": "ModelID"})
        .merge(current, on="ModelID")
        .merge(signatures, on="ModelID")
    )

    data = data[
        data["CCLE_MSI"] == "MSI"
    ].dropna(
        subset=[
            "FPBM_mean",
            "WRN_26Q3",
            "MSIScore",
        ]
    ).reset_index(drop=True)

    assert len(data) == 26

    return data


def residual_corr(df):
    msi = df["MSIScore"].to_numpy(float)
    fpbm = df["FPBM_mean"].to_numpy(float)
    wrn = df["WRN_26Q3"].to_numpy(float)

    design = np.column_stack(
        [
            np.ones(len(df)),
            msi,
        ]
    )

    fpbm_resid = fpbm - design @ np.linalg.lstsq(
        design,
        fpbm,
        rcond=None,
    )[0]

    wrn_resid = wrn - design @ np.linalg.lstsq(
        design,
        wrn,
        rcond=None,
    )[0]

    return (
        pearsonr(fpbm_resid, wrn_resid)[0],
        spearmanr(fpbm_resid, wrn_resid)[0],
    )


data = load_data()

full_pearson, full_spearman = residual_corr(data)

print("ROBUSTNESS OF BROKEN-TA SIGNAL BEYOND MSI")
print("=" * 65)
print(f"N: {len(data)}")
print(f"Full partial Pearson:  {full_pearson:.3f}")
print(f"Full residual Spearman:{full_spearman:.3f}")


# Leave-one-out influence analysis
loo = []

for i in range(len(data)):
    subset = data.drop(index=i).reset_index(drop=True)

    pr, sr = residual_corr(subset)

    loo.append(
        {
            "removed": data.loc[i, "ModelID"],
            "pearson": pr,
            "spearman": sr,
            "FPBM": data.loc[i, "FPBM_mean"],
            "WRN": data.loc[i, "WRN_26Q3"],
            "MSI": data.loc[i, "MSIScore"],
        }
    )

loo = pd.DataFrame(loo)

print()
print("Leave-one-out range")
print("=" * 65)
print(
    f"Pearson:  {loo['pearson'].min():.3f} "
    f"to {loo['pearson'].max():.3f}"
)
print(
    f"Spearman: {loo['spearman'].min():.3f} "
    f"to {loo['spearman'].max():.3f}"
)

print()
print("Most influential removals")
print("=" * 65)

order = (
    loo["pearson"] - full_pearson
).abs().sort_values(
    ascending=False
).index

for _, row in loo.loc[order].head(8).iterrows():
    print(
        f"{row['removed']:<12} "
        f"Pearson_without={row['pearson']:>6.3f}  "
        f"Spearman_without={row['spearman']:>6.3f}  "
        f"FPBM={row['FPBM']:>6.3f}  "
        f"WRN={row['WRN']:>7.3f}  "
        f"MSI={row['MSI']:>7.3f}"
    )


# Bootstrap partial Pearson
rng = np.random.default_rng(SEED)
boots = []

for _ in range(BOOTSTRAPS):
    idx = rng.integers(
        0,
        len(data),
        len(data),
    )

    sample = data.iloc[idx].reset_index(drop=True)

    try:
        r, _ = residual_corr(sample)

        if np.isfinite(r):
            boots.append(r)
    except Exception:
        pass

low, median, high = np.percentile(
    boots,
    [2.5, 50, 97.5],
)

print()
print("Bootstrap partial Pearson")
print("=" * 65)
print(f"Median:  {median:.3f}")
print(f"95% CI:  {low:.3f} to {high:.3f}")
print(f"Bootstrap samples: {len(boots)}")
