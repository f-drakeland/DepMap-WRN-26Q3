from pathlib import Path

import pandas as pd
from scipy.stats import linregress, pearsonr, spearmanr


INPUT = Path("results/current-wgs/cohort_summary.csv")
OUTPUT = Path("results/current-wgs/wrn_association_analysis.txt")
DIAGNOSTICS = Path("results/current-wgs/wrn_association_diagnostics.csv")


df = pd.read_csv(INPUT)

required = [
    "ModelID",
    "CellLineName",
    "PublishedFPBM",
    "CurrentBrokenTAFPBM",
    "WRNGeneEffect",
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
        f"Expected 26 models, found {len(df)}"
    )

if not df["Status"].eq("complete").all():
    raise RuntimeError(
        "Cohort contains incomplete models"
    )

analysis = df.dropna(
    subset=[
        "PublishedFPBM",
        "CurrentBrokenTAFPBM",
        "WRNGeneEffect",
    ]
).copy()

if len(analysis) != 26:
    raise RuntimeError(
        f"Expected 26 analyzable models, found {len(analysis)}"
    )


def correlations(x, y):
    pearson_r, pearson_p = pearsonr(x, y)
    spearman_r, spearman_p = spearmanr(x, y)

    return {
        "pearson_r": pearson_r,
        "pearson_p": pearson_p,
        "spearman_r": spearman_r,
        "spearman_p": spearman_p,
    }


current = correlations(
    analysis["CurrentBrokenTAFPBM"],
    analysis["WRNGeneEffect"],
)

published = correlations(
    analysis["PublishedFPBM"],
    analysis["WRNGeneEffect"],
)

fit = linregress(
    analysis["CurrentBrokenTAFPBM"],
    analysis["WRNGeneEffect"],
)

analysis["PredictedWRNGeneEffect"] = (
    fit.intercept
    + fit.slope
    * analysis["CurrentBrokenTAFPBM"]
)

analysis["WRNResidual"] = (
    analysis["WRNGeneEffect"]
    - analysis["PredictedWRNGeneEffect"]
)

analysis["AbsoluteWRNResidual"] = (
    analysis["WRNResidual"].abs()
)

loo_rows = []

for idx in analysis.index:
    keep = analysis.index != idx

    pearson_r, _ = pearsonr(
        analysis.loc[keep, "CurrentBrokenTAFPBM"],
        analysis.loc[keep, "WRNGeneEffect"],
    )

    spearman_r, _ = spearmanr(
        analysis.loc[keep, "CurrentBrokenTAFPBM"],
        analysis.loc[keep, "WRNGeneEffect"],
    )

    loo_rows.append(
        {
            "RemovedModelID":
                analysis.loc[idx, "ModelID"],
            "RemovedCellLine":
                analysis.loc[idx, "CellLineName"],
            "PearsonR":
                pearson_r,
            "SpearmanR":
                spearman_r,
        }
    )

loo = pd.DataFrame(loo_rows)

pearson_low = loo.loc[
    loo["PearsonR"].idxmin()
]

pearson_high = loo.loc[
    loo["PearsonR"].idxmax()
]

spearman_low = loo.loc[
    loo["SpearmanR"].idxmin()
]

spearman_high = loo.loc[
    loo["SpearmanR"].idxmax()
]

largest_residual = analysis.loc[
    analysis["AbsoluteWRNResidual"].idxmax()
]

diagnostics = analysis[
    [
        "ModelID",
        "CellLineName",
        "Lineage",
        "PublishedFPBM",
        "CurrentBrokenTAFPBM",
        "WRNGeneEffect",
        "PredictedWRNGeneEffect",
        "WRNResidual",
        "AbsoluteWRNResidual",
    ]
].sort_values(
    "AbsoluteWRNResidual",
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


emit("CURRENT-WGS FPBM VS CURRENT WRN DEPENDENCY")
emit("=" * 72)
emit(f"Models: {len(analysis)}")
emit(f"Input: {INPUT}")
emit()
emit("Interpretation of direction")
emit("-" * 72)
emit(
    "Lower FPBM = greater broken-TA coverage loss."
)
emit(
    "More-negative WRN GeneEffect = stronger WRN dependency."
)
emit(
    "Therefore a positive correlation is the expected direction."
)
emit()
emit("Reconstructed current-WGS FPBM vs current WRN")
emit("-" * 72)
emit(
    f"Pearson r:  {current['pearson_r']:.6f}   "
    f"p={current['pearson_p']:.6g}"
)
emit(
    f"Spearman r: {current['spearman_r']:.6f}   "
    f"p={current['spearman_p']:.6g}"
)
emit()
emit("Published FPBM vs the same current WRN values")
emit("-" * 72)
emit(
    f"Pearson r:  {published['pearson_r']:.6f}   "
    f"p={published['pearson_p']:.6g}"
)
emit(
    f"Spearman r: {published['spearman_r']:.6f}   "
    f"p={published['spearman_p']:.6g}"
)
emit()
emit("Descriptive coefficient difference")
emit("-" * 72)
emit(
    "Current minus published Pearson: "
    f"{current['pearson_r'] - published['pearson_r']:+.6f}"
)
emit(
    "Current minus published Spearman: "
    f"{current['spearman_r'] - published['spearman_r']:+.6f}"
)
emit(
    "These coefficient differences are descriptive only; "
    "no test of dependent correlations is claimed here."
)
emit()
emit("Current FPBM linear fit to WRN")
emit("-" * 72)
emit(
    f"WRN GeneEffect = {fit.intercept:.6f} + "
    f"{fit.slope:.6f} * Current FPBM"
)
emit(
    f"R^2: {fit.rvalue ** 2:.6f}"
)
emit()
emit("Largest WRN residual from current-FPBM fit")
emit("-" * 72)
emit(
    f"{largest_residual['ModelID']} "
    f"{largest_residual['CellLineName']}"
)
emit(
    f"Current FPBM: {largest_residual['CurrentBrokenTAFPBM']:.6f}"
)
emit(
    f"Observed WRN: {largest_residual['WRNGeneEffect']:.6f}"
)
emit(
    f"Residual:     {largest_residual['WRNResidual']:+.6f}"
)
emit()
emit("Leave-one-out sensitivity: current FPBM vs WRN")
emit("-" * 72)
emit(
    f"Pearson range:  "
    f"{loo['PearsonR'].min():.6f} to "
    f"{loo['PearsonR'].max():.6f}"
)
emit(
    "Lowest Pearson after removing: "
    f"{pearson_low['RemovedModelID']} "
    f"{pearson_low['RemovedCellLine']}"
)
emit(
    "Highest Pearson after removing: "
    f"{pearson_high['RemovedModelID']} "
    f"{pearson_high['RemovedCellLine']}"
)
emit(
    f"Spearman range: "
    f"{loo['SpearmanR'].min():.6f} to "
    f"{loo['SpearmanR'].max():.6f}"
)
emit(
    "Lowest Spearman after removing: "
    f"{spearman_low['RemovedModelID']} "
    f"{spearman_low['RemovedCellLine']}"
)
emit(
    "Highest Spearman after removing: "
    f"{spearman_high['RemovedModelID']} "
    f"{spearman_high['RemovedCellLine']}"
)
emit()
emit(f"Diagnostics: {DIAGNOSTICS}")

OUTPUT.write_text(
    "\n".join(lines) + "\n"
)
