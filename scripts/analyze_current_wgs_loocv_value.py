from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, wilcoxon
from sklearn.linear_model import LinearRegression


COHORT = Path("results/current-wgs/cohort_summary.csv")
SIGNATURES = Path("data/source/OmicsGlobalSignatures.csv")

OUTPUT = Path(
    "results/current-wgs/loocv_prediction_analysis.txt"
)

PREDICTIONS = Path(
    "results/current-wgs/loocv_predictions.csv"
)

N_PERM = 100_000
SEED = 42


def load_data():
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
            f"Expected 26 selected MSI rows, found {len(msi)}"
        )

    if msi["ModelID"].nunique() != 26:
        raise RuntimeError(
            "Expected exactly one MSI score per cohort model"
        )

    data = cohort.merge(
        msi,
        on="ModelID",
        how="inner",
        validate="one_to_one",
    )

    data = data.dropna(
        subset=[
            "CurrentBrokenTAFPBM",
            "PublishedFPBM",
            "WRNGeneEffect",
            "MSIScore",
        ]
    ).reset_index(drop=True)

    if len(data) != 26:
        raise RuntimeError(
            f"Expected 26 analyzable models, found {len(data)}"
        )

    return data


def loocv_predict(df, feature_columns):
    y = df[
        "WRNGeneEffect"
    ].to_numpy(float)

    x = df[
        feature_columns
    ].to_numpy(float)

    pred = np.zeros(
        len(df),
        dtype=float,
    )

    for i in range(len(df)):
        train = (
            np.arange(len(df)) != i
        )

        model = LinearRegression().fit(
            x[train],
            y[train],
        )

        pred[i] = model.predict(
            x[[i]]
        )[0]

    return pred


def metrics(y, pred):
    error = y - pred
    squared_error = error ** 2

    rmse = np.sqrt(
        squared_error.mean()
    )

    mae = np.abs(
        error
    ).mean()

    r = pearsonr(
        y,
        pred,
    )[0]

    denominator = np.sum(
        (y - y.mean()) ** 2
    )

    predictive_r2 = (
        1
        - np.sum(squared_error)
        / denominator
    )

    return {
        "error": error,
        "squared_error": squared_error,
        "rmse": rmse,
        "mae": mae,
        "pearson_r": r,
        "predictive_r2": predictive_r2,
    }


def paired_comparison(
    baseline_error,
    candidate_error,
):
    improvement = (
        baseline_error
        - candidate_error
    )

    observed = improvement.mean()

    try:
        _, wilcoxon_p = wilcoxon(
            baseline_error,
            candidate_error,
            alternative="greater",
        )
    except ValueError:
        wilcoxon_p = np.nan

    rng = np.random.default_rng(
        SEED
    )

    signs = rng.choice(
        [-1.0, 1.0],
        size=(
            N_PERM,
            len(improvement),
        ),
    )

    permuted = (
        signs
        * improvement
    ).mean(axis=1)

    permutation_p = (
        np.count_nonzero(
            permuted >= observed
        )
        + 1
    ) / (
        N_PERM + 1
    )

    return {
        "mean_mse_improvement":
            observed,
        "wilcoxon_p":
            wilcoxon_p,
        "permutation_p":
            permutation_p,
        "models_improved":
            int(
                np.sum(
                    candidate_error
                    < baseline_error
                )
            ),
    }


def analyze_set(label, df):
    y = df[
        "WRNGeneEffect"
    ].to_numpy(float)

    pred_msi = loocv_predict(
        df,
        [
            "MSIScore",
        ],
    )

    pred_current = loocv_predict(
        df,
        [
            "MSIScore",
            "CurrentBrokenTAFPBM",
        ],
    )

    pred_published = loocv_predict(
        df,
        [
            "MSIScore",
            "PublishedFPBM",
        ],
    )

    msi = metrics(
        y,
        pred_msi,
    )

    current = metrics(
        y,
        pred_current,
    )

    published = metrics(
        y,
        pred_published,
    )

    current_compare = paired_comparison(
        msi["squared_error"],
        current["squared_error"],
    )

    published_compare = paired_comparison(
        msi["squared_error"],
        published["squared_error"],
    )

    rows = df[
        [
            "ModelID",
            "CellLineName",
            "Lineage",
            "MSIScore",
            "CurrentBrokenTAFPBM",
            "PublishedFPBM",
            "WRNGeneEffect",
        ]
    ].copy()

    rows.insert(
        0,
        "AnalysisSet",
        label,
    )

    rows[
        "PredictionMSIOnly"
    ] = pred_msi

    rows[
        "PredictionMSIPlusCurrentFPBM"
    ] = pred_current

    rows[
        "PredictionMSIPlusPublishedFPBM"
    ] = pred_published

    rows[
        "SquaredErrorMSIOnly"
    ] = msi["squared_error"]

    rows[
        "SquaredErrorMSIPlusCurrentFPBM"
    ] = current["squared_error"]

    rows[
        "SquaredErrorMSIPlusPublishedFPBM"
    ] = published["squared_error"]

    return {
        "label": label,
        "n": len(df),
        "msi": msi,
        "current": current,
        "published": published,
        "current_compare":
            current_compare,
        "published_compare":
            published_compare,
        "rows": rows,
    }


data = load_data()

all_models = analyze_set(
    "ALL_26",
    data,
)

without_km12 = analyze_set(
    "WITHOUT_KM12",
    data[
        data["ModelID"]
        != "ACH-000969"
    ].reset_index(drop=True),
)

pd.concat(
    [
        all_models["rows"],
        without_km12["rows"],
    ],
    ignore_index=True,
).to_csv(
    PREDICTIONS,
    index=False,
)


lines = []


def emit(text=""):
    lines.append(text)
    print(text)


def report(result):
    msi = result["msi"]
    current = result["current"]
    published = result["published"]

    current_compare = (
        result["current_compare"]
    )

    published_compare = (
        result["published_compare"]
    )

    emit(result["label"])
    emit("=" * 72)
    emit(f"Models: {result['n']}")
    emit()

    emit("MSI-only")
    emit("-" * 72)
    emit(
        f"RMSE:          {msi['rmse']:.6f}"
    )
    emit(
        f"MAE:           {msi['mae']:.6f}"
    )
    emit(
        f"Prediction r:  {msi['pearson_r']:.6f}"
    )
    emit(
        f"Predictive R2: {msi['predictive_r2']:.6f}"
    )
    emit()

    emit("MSI + reconstructed current-WGS FPBM")
    emit("-" * 72)
    emit(
        f"RMSE:          {current['rmse']:.6f}"
    )
    emit(
        f"RMSE change:   "
        f"{current['rmse'] - msi['rmse']:+.6f}"
    )
    emit(
        f"Percent change:"
        f" {(current['rmse'] - msi['rmse']) / msi['rmse'] * 100:+.2f}%"
    )
    emit(
        f"MAE:           {current['mae']:.6f}"
    )
    emit(
        f"Prediction r:  {current['pearson_r']:.6f}"
    )
    emit(
        f"r change:      "
        f"{current['pearson_r'] - msi['pearson_r']:+.6f}"
    )
    emit(
        f"Predictive R2: {current['predictive_r2']:.6f}"
    )
    emit(
        "Models with lower squared error: "
        f"{current_compare['models_improved']}"
        f"/{result['n']}"
    )
    emit(
        "Mean MSE improvement: "
        f"{current_compare['mean_mse_improvement']:+.6f}"
    )
    emit(
        f"Wilcoxon p:    "
        f"{current_compare['wilcoxon_p']:.6g}"
    )
    emit(
        f"Sign-flip p:   "
        f"{current_compare['permutation_p']:.6f}"
    )
    emit()

    emit("MSI + published FPBM reference")
    emit("-" * 72)
    emit(
        f"RMSE:          {published['rmse']:.6f}"
    )
    emit(
        f"RMSE change:   "
        f"{published['rmse'] - msi['rmse']:+.6f}"
    )
    emit(
        f"Percent change:"
        f" {(published['rmse'] - msi['rmse']) / msi['rmse'] * 100:+.2f}%"
    )
    emit(
        f"MAE:           {published['mae']:.6f}"
    )
    emit(
        f"Prediction r:  {published['pearson_r']:.6f}"
    )
    emit(
        f"r change:      "
        f"{published['pearson_r'] - msi['pearson_r']:+.6f}"
    )
    emit(
        f"Predictive R2: {published['predictive_r2']:.6f}"
    )
    emit(
        "Models with lower squared error: "
        f"{published_compare['models_improved']}"
        f"/{result['n']}"
    )
    emit(
        "Mean MSE improvement: "
        f"{published_compare['mean_mse_improvement']:+.6f}"
    )
    emit(
        f"Wilcoxon p:    "
        f"{published_compare['wilcoxon_p']:.6g}"
    )
    emit(
        f"Sign-flip p:   "
        f"{published_compare['permutation_p']:.6f}"
    )
    emit()


emit("LOOCV VALUE OF FPBM BEYOND MSI")
emit("=" * 72)
emit(
    "Outcome: current WRN GeneEffect"
)
emit(
    "Lower RMSE is better."
)
emit(
    "Negative RMSE change means FPBM improved prediction."
)
emit()
emit(
    "Paired Wilcoxon and sign-flip p-values compare "
    "held-out squared errors. Because LOOCV folds share "
    "training observations, these are sensitivity measures "
    "rather than independent-fold hypothesis tests."
)
emit()

report(all_models)
report(without_km12)

emit(f"Permutation seed: {SEED}")
emit(f"Permutations: {N_PERM}")
emit(f"Predictions: {PREDICTIONS}")

OUTPUT.write_text(
    "\n".join(lines) + "\n"
)
