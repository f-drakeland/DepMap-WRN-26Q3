from pathlib import Path

import pandas as pd
from scipy.stats import pearsonr, spearmanr


SOURCE = Path("data/source/WRN_2020_SourceData_Fig4.xlsx")
GENE_EFFECT = Path("data/source/CRISPRGeneEffect.csv")
WRN_COLUMN = "WRN (7486)"


def report(label, df):
    fpbm_pr, fpbm_pp = pearsonr(
        df["FPBM_mean"],
        df["WRN_26Q3"],
    )
    fpbm_sr, fpbm_sp = spearmanr(
        df["FPBM_mean"],
        df["WRN_26Q3"],
    )

    wrn_pr, wrn_pp = pearsonr(
        df["avg_WRN_dep"],
        df["WRN_26Q3"],
    )
    wrn_sr, wrn_sp = spearmanr(
        df["avg_WRN_dep"],
        df["WRN_26Q3"],
    )

    print(label)
    print("=" * 60)
    print(f"N: {len(df)}")

    print("\nPublished broken-TA FPBM vs 26Q3 WRN")
    print(f"Pearson:   {fpbm_pr:.3f}   p={fpbm_pp:.3g}")
    print(f"Spearman:  {fpbm_sr:.3f}   p={fpbm_sp:.3g}")

    print("\nPublished WRN dependency vs 26Q3 WRN")
    print(f"Pearson:   {wrn_pr:.3f}   p={wrn_pp:.3g}")
    print(f"Spearman:  {wrn_sr:.3f}   p={wrn_sp:.3g}")
    print()


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

current.columns = [
    "DepMap_ID",
    "WRN_26Q3",
]

current["WRN_26Q3"] = pd.to_numeric(
    current["WRN_26Q3"],
    errors="coerce",
)

data = published.merge(
    current,
    on="DepMap_ID",
    how="inner",
).dropna(
    subset=[
        "FPBM_mean",
        "avg_WRN_dep",
        "WRN_26Q3",
    ]
)

assert len(data) == 204
assert (data["CCLE_MSI"] == "MSI").sum() == 26
assert (data["CCLE_MSI"] == "MSS").sum() == 178

report("ALL OVERLAPPING MODELS", data)
report("MSI", data[data["CCLE_MSI"] == "MSI"])
report("MSS", data[data["CCLE_MSI"] == "MSS"])

print("Overlap summary")
print("=" * 60)
print("Published models:", len(published))
print("Matched to 26Q3:", len(data))
print("MSI matched:", (data["CCLE_MSI"] == "MSI").sum())
print("MSS matched:", (data["CCLE_MSI"] == "MSS").sum())
