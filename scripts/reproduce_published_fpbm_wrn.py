from pathlib import Path

import pandas as pd
from scipy.stats import pearsonr, spearmanr


SOURCE = Path("data/source/WRN_2020_SourceData_Fig4.xlsx")


def report(label, df):
    pearson_r, pearson_p = pearsonr(
        df["FPBM_mean"],
        df["avg_WRN_dep"],
    )
    spearman_r, spearman_p = spearmanr(
        df["FPBM_mean"],
        df["avg_WRN_dep"],
    )

    print(label)
    print("=" * 60)
    print(f"N:          {len(df)}")
    print(f"Pearson:    {pearson_r:.3f}   p={pearson_p:.3g}")
    print(f"Spearman:   {spearman_r:.3f}   p={spearman_p:.3g}")
    print()


data = pd.read_excel(
    SOURCE,
    sheet_name="4d",
    engine="openpyxl",
).dropna(
    subset=[
        "DepMap_ID",
        "CCLE_MSI",
        "avg_WRN_dep",
        "FPBM_mean",
    ]
)

# Guard against silently analyzing a different source table.
assert len(data) == 238
assert (data["CCLE_MSI"] == "MSI").sum() == 31
assert (data["CCLE_MSI"] == "MSS").sum() == 207

report("ALL CELL LINES", data)
report("MSI", data[data["CCLE_MSI"] == "MSI"])
report("MSS", data[data["CCLE_MSI"] == "MSS"])
