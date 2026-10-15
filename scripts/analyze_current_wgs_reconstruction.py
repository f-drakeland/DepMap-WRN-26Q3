from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import linregress, pearsonr, spearmanr


INPUT = Path("results/current-wgs/cohort_summary.csv")
OUTPUT = Path("results/current-wgs/reconstruction_analysis.txt")
DIAGNOSTICS = Path("results/current-wgs/reconstruction_diagnostics.csv")


df = pd.read_csv(INPUT)

required = [
    "ModelID",
    "CellLineName",
    "PublishedFPBM",
    "CurrentBrokenTAFPBM",
    "CurrentToPublishedRatio",
]

missing = [
    column
    for column in required
    if column not in df.columns
]

if missing:
    raise RuntimeError(
        f"Missing required columns: {missing}"
    )

if len(df) != 26:
    raise RuntimeError(
        f"Expected 26 cohort models, found {len(df)}"
    )

if not df["Status"].eq("complete").all():
    raise RuntimeError(
        "Cohort contains incomplete models"
    )

analysis = df.dropna(
    subset=[
        "PublishedFPBM",
        "CurrentBrokenTAFPBM",
    ]
).copy()

if len(analysis) != 26:
    raise RuntimeError(
        f"Expected 26 analyzable models, found {len(analysis)}"
    )

published = analysis["PublishedFPBM"]
current = analysis["CurrentBrokenTAFPBM"]

pearson_r, pearson_p = pearsonr(
    published,
    current,
)

spearman_r, spearman_p = spearmanr(
    published,
    current,
)

fit = linregress(
    published,
    current,
)

analysis["Difference"] = (
    analysis["CurrentBrokenTAFPBM"]
    - analysis["PublishedFPBM"]
)

analysis["AbsoluteDifference"] = (
    analysis["Difference"].abs()
)

analysis["FitPredictedCurrent"] = (
    fit.intercept
    + fit.slope
    * analysis["PublishedFPBM"]
)

analysis["FitResidual"] = (
    analysis["CurrentBrokenTAFPBM"]
    - analysis["FitPredictedCurrent"]
)

analysis["AbsoluteFitResidual"] = (
    analysis["FitResidual"].abs()
)

ratios = analysis["CurrentToPublishedRatio"]

loo_rows = []

for idx in analysis.index:
    keep = analysis.index != idx

    r, _ = pearsonr(
        analysis.loc[keep, "PublishedFPBM"],
        analysis.loc[keep, "CurrentBrokenTAFPBM"],
    )

    loo_rows.append(
        (
            analysis.loc[idx, "ModelID"],
            analysis.loc[idx, "CellLineName"],
            r,
        )
    )

loo = pd.DataFrame(
    loo_rows,
    columns=[
        "RemovedModelID",
        "RemovedCellLine",
        "PearsonR",
    ],
)

loo_min = loo.loc[
    loo["PearsonR"].idxmin()
]

loo_max = loo.loc[
    loo["PearsonR"].idxmax()
]

largest_difference = analysis.loc[
    analysis["AbsoluteDifference"].idxmax()
]

largest_residual = analysis.loc[
    analysis["AbsoluteFitResidual"].idxmax()
]

diagnostics = analysis[
    [
        "ModelID",
        "CellLineName",
        "Lineage",
        "PublishedFPBM",
        "CurrentBrokenTAFPBM",
        "CurrentNonbrokenTAFPBM",
        "CurrentToPublishedRatio",
        "Difference",
        "FitPredictedCurrent",
        "FitResidual",
        "AbsoluteFitResidual",
    ]
].sort_values(
    "AbsoluteFitResidual",
    ascending=False,
)

diagnostics.to_csv(
    DIAGNOSTICS,
    index=False,
)

lines = []

def emit(text=""):
    lines.append(text)
    print(text)


emit("CURRENT-WGS FPBM RECONSTRUCTION")
emit("=" * 72)
emit(f"Models: {len(analysis)}")
emit(f"Input: {INPUT}")
emit()
emit("Published FPBM vs reconstructed current-WGS FPBM")
emit("-" * 72)
emit(
    f"Pearson r:  {pearson_r:.6f}   p={pearson_p:.6g}"
)
emit(
    f"Spearman r: {spearman_r:.6f}   p={spearman_p:.6g}"
)
emit(
    f"R^2:        {fit.rvalue ** 2:.6f}"
)
emit()
emit("Linear calibration")
emit("-" * 72)
emit(
    f"Current = {fit.intercept:.6f} + "
    f"{fit.slope:.6f} * Published"
)
emit(
    f"Slope p: {fit.pvalue:.6g}"
)
emit()
emit("Current / published ratio")
emit("-" * 72)
emit(
    f"Mean:   {ratios.mean():.6f}"
)
emit(
    f"Median: {ratios.median():.6f}"
)
emit(
    f"Min:    {ratios.min():.6f}"
)
emit(
    f"Max:    {ratios.max():.6f}"
)
emit()
emit("Current - published difference")
emit("-" * 72)
emit(
    f"Mean:   {analysis['Difference'].mean():.6f}"
)
emit(
    f"Median: {analysis['Difference'].median():.6f}"
)
emit(
    f"MAE:    {analysis['AbsoluteDifference'].mean():.6f}"
)
emit()
emit("Largest raw difference")
emit("-" * 72)
emit(
    f"{largest_difference['ModelID']} "
    f"{largest_difference['CellLineName']}"
)
emit(
    f"Published: {largest_difference['PublishedFPBM']:.6f}"
)
emit(
    f"Current:   {largest_difference['CurrentBrokenTAFPBM']:.6f}"
)
emit(
    f"Difference:{largest_difference['Difference']:+.6f}"
)
emit(
    f"Ratio:     {largest_difference['CurrentToPublishedRatio']:.6f}"
)
emit()
emit("Largest regression residual")
emit("-" * 72)
emit(
    f"{largest_residual['ModelID']} "
    f"{largest_residual['CellLineName']}"
)
emit(
    f"Residual: {largest_residual['FitResidual']:+.6f}"
)
emit()
emit("Leave-one-out Pearson sensitivity")
emit("-" * 72)
emit(
    f"Range: {loo['PearsonR'].min():.6f} "
    f"to {loo['PearsonR'].max():.6f}"
)
emit(
    "Lowest after removing: "
    f"{loo_min['RemovedModelID']} "
    f"{loo_min['RemovedCellLine']}"
)
emit(
    "Highest after removing: "
    f"{loo_max['RemovedModelID']} "
    f"{loo_max['RemovedCellLine']}"
)
emit()
emit(f"Diagnostics: {DIAGNOSTICS}")

OUTPUT.write_text(
    "\n".join(lines) + "\n"
)
