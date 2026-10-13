from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, wilcoxon
from sklearn.linear_model import LinearRegression


SOURCE = Path("data/source/WRN_2020_SourceData_Fig4.xlsx")
GENE_EFFECT = Path("data/source/CRISPRGeneEffect.csv")
SIGNATURES = Path("data/source/OmicsGlobalSignatures.csv")

WRN_COLUMN = "WRN (7486)"
N_PERM = 100_000
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


def loocv(df):
    y = df["WRN_26Q3"].to_numpy(float)

    pred_msi = np.zeros(len(df))
    pred_both = np.zeros(len(df))

    x_msi = df[
        ["MSIScore"]
    ].to_numpy(float)

    x_both = df[
        ["MSIScore", "FPBM_mean"]
    ].to_numpy(float)

    for i in range(len(df)):
        train = np.arange(len(df)) != i

        model_msi = LinearRegression().fit(
            x_msi[train],
            y[train],
        )

        model_both = LinearRegression().fit(
            x_both[train],
            y[train],
        )

        pred_msi[i] = model_msi.predict(
            x_msi[[i]]
        )[0]

        pred_both[i] = model_both.predict(
            x_both[[i]]
        )[0]

    error_msi = (y - pred_msi) ** 2
    error_both = (y - pred_both) ** 2

    rmse_msi = np.sqrt(error_msi.mean())
    rmse_both = np.sqrt(error_both.mean())

    r_msi = pearsonr(
        y,
        pred_msi,
    )[0]

    r_both = pearsonr(
        y,
        pred_both,
    )[0]

    return (
        y,
        pred_msi,
        pred_both,
        error_msi,
        error_both,
        rmse_msi,
        rmse_both,
        r_msi,
        r_both,
    )


def test(label, df):
    (
        y,
        pred_msi,
        pred_both,
        error_msi,
        error_both,
        rmse_msi,
        rmse_both,
        r_msi,
        r_both,
    ) = loocv(df)

    improvement = error_msi - error_both

    print(label)
    print("=" * 65)
    print("N:", len(df))
    print(f"MSI-only RMSE:       {rmse_msi:.4f}")
    print(f"MSI+FPBM RMSE:       {rmse_both:.4f}")
    print(f"RMSE change:         {rmse_both-rmse_msi:+.4f}")
    print(
        "Percent RMSE change: "
        f"{(rmse_both-rmse_msi)/rmse_msi*100:+.2f}%"
    )
    print()
    print(f"MSI-only r:          {r_msi:.3f}")
    print(f"MSI+FPBM r:          {r_both:.3f}")
    print(f"Correlation change:  {r_both-r_msi:+.3f}")
    print()

    _, wilcoxon_p = wilcoxon(
        error_msi,
        error_both,
        alternative="greater",
    )

    print(f"Wilcoxon p:          {wilcoxon_p:.4g}")

    observed = improvement.mean()

    rng = np.random.default_rng(SEED)
    extreme = 0

    for _ in range(N_PERM):
        signs = rng.choice(
            [-1, 1],
            size=len(improvement),
        )

        permuted = np.mean(
            improvement * signs
        )

        if permuted >= observed:
            extreme += 1

    permutation_p = (
        extreme + 1
    ) / (
        N_PERM + 1
    )

    print(f"Permutation p:       {permutation_p:.5f}")
    print(f"Mean MSE improvement:{observed:.4f}")
    print()


data = load_data()

test(
    "ALL 26 MSI MODELS",
    data,
)

test(
    "WITHOUT ACH-000969",
    data[
        data["ModelID"] != "ACH-000969"
    ].reset_index(drop=True),
)
