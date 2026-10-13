import pandas as pd
import openpyxl

COVERAGE_FILE = "results/km12/KM12_all_ta.bedcov.tsv"
SOURCE_WORKBOOK = "data/source/WRN_2020_SourceData_Fig4.xlsx"
MODEL_ID = "ACH-000969"

x = pd.read_csv(
    COVERAGE_FILE,
    sep="\t",
    header=None,
    names=["chrom", "start", "end", "locus", "bedcov"],
)

x["length"] = x["end"] - x["start"]

# Published method:
# FPB = coverage / (interval length * 2)
x["FPB"] = x["bedcov"] / (x["length"] * 2.0)

total_fpb = x["FPB"].sum()
x["FPBM"] = (x["FPB"] / total_fpb) * 1_000_000

# Support both naming schemes used during exploration.
broken = x[
    x["locus"].str.startswith("brokenTA_")
    | x["locus"].str.startswith("B_")
].copy()

nonbroken = x[
    x["locus"].str.startswith("NB_")
].copy()

if len(broken) != 5377:
    raise RuntimeError(
        f"Expected 5377 broken loci, found {len(broken)}"
    )

if len(nonbroken) != 60985:
    raise RuntimeError(
        f"Expected 60985 non-broken loci, found {len(nonbroken)}"
    )

wb = openpyxl.load_workbook(
    SOURCE_WORKBOOK,
    read_only=True,
    data_only=True,
)
ws = wb["4d"]

published = None

for row in ws.iter_rows(min_row=2, values_only=True):
    if row[0] == MODEL_ID:
        published = float(row[4])
        break

if published is None:
    raise RuntimeError(f"{MODEL_ID} not found in published source data")

current = broken["FPBM"].mean()

print("KM12 CURRENT-WGS BROKEN-TA RECONSTRUCTION")
print("=" * 60)
print(f"Total TA intervals:              {len(x)}")
print(f"Broken intervals:                {len(broken)}")
print(f"Non-broken intervals:            {len(nonbroken)}")
print()
print(f"Current broken-TA mean FPBM:     {current:.4f}")
print(f"Current non-broken mean FPBM:    {nonbroken['FPBM'].mean():.4f}")
print(f"Published KM12 FPBM_mean:        {published:.4f}")
print()
print(f"Current / published ratio:       {current/published:.3f}")
