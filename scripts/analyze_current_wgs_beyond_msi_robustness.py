from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr, t


COHORT = Path("results/current-wgs/cohort_summary.csv")
SIGNATURES = Path("data/source/OmicsGlobalSignatures.csv")

OUTPUT = Path(
    "results/current-wgs/beyond_msi_robustness.txt"
)

LOO_OUTPUT = Path(
    "results/current-wgs/beyond_msi_leave_one_out.csv"
)

BOOTSTRAPS = 20_000
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


def partial_stats(df, fpbm_column):
    fpbm = df[fpbm_column].to_numpy(float)
    wrn = df["WRNGeneEffect"].to_numpy(float)
    msi = df["MSIScore"].to_numpy(float)

    fpbm_resid = residualize(
        fpbm,
        msi,
    )

    wrn_resid = residualize(
        wrn,
        msi,
    )

    pearson_r = pearsonr(
        fpbm_resid,
        wrn_resid,
    )[0]

    dfree = len(df) - 3

    pearson_t = (
        pearson_r
        * np.sqrt(
            dfree
            / (1 - pearson_r ** 2)
        )
    )

    pearson_p = (
        2
        * t.sf(
            abs(pearson_t),
            df=dfree,
        )
    )

    spearman_r, spearman_p = spearmanr(
        fpbm_resid,
        wrn_resid,
    )

    return {
        "pearson_r": pearson_r,
        "pearson_p": pearson_p,
        "spearman_r": spearman_r,
        "spearman_p": spearman_p,
    }


def partial_r_only(df):
    fpbm = df[
        "CurrentBrokenTAFPBM"
    ].to_numpy(float)

    wrn = df[
        "WRNGeneEffect"
    ].to_numpy(float)

    msi = df[
        "MSIScore"
    ].to_numpy(float)

    fpbm_resid = residualize(
        fpbm,
        msi,
    )

    wrn_resid = residualize(
        wrn,
        msi,
    )

    if (
        np.std(fpbm_resid) == 0
        or np.std(wrn_resid) == 0
    ):
        return np.nan

    return np.corrcoef(
        fpbm_resid,
        wrn_resid,
    )[0, 1]


data = load_data()

full = partial_stats(
    data,
    "CurrentBrokenTAFPBM",
)

published_full = partial_stats(
    data,
    "PublishedFPBM",
)


# Leave-one-out influence.
loo_rows = []

for i in range(len(data)):
    subset = (
        data.drop(index=i)
        .reset_index(drop=True)
    )

    current = partial_stats(
        subset,
        "CurrentBrokenTAFPBM",
    )

    published = partial_stats(
        subset,
        "PublishedFPBM",
    )

    loo_rows.append(
        {
            "RemovedModelID":
                data.loc[i, "ModelID"],
            "RemovedCellLine":
                data.loc[i, "CellLineName"],
            "Lineage":
                data.loc[i, "Lineage"],
            "CurrentFPBM":
                data.loc[
                    i,
                    "CurrentBrokenTAFPBM",
                ],
            "PublishedFPBM":
                data.loc[
                    i,
                    "PublishedFPBM",
                ],
            "WRNGeneEffect":
                data.loc[
                    i,
                    "WRNGeneEffect",
                ],
            "MSIScore":
                data.loc[
                    i,
                    "MSIScore",
                ],
            "CurrentPartialPearson":
                current["pearson_r"],
            "CurrentPartialPearsonP":
                current["pearson_p"],
            "CurrentResidualSpearman":
                current["spearman_r"],
            "PublishedPartialPearson":
                published["pearson_r"],
            "DeltaFromFullCurrent":
                (
                    current["pearson_r"]
                    - full["pearson_r"]
                ),
        }
    )

loo = pd.DataFrame(loo_rows)

loo["AbsoluteInfluence"] = (
    loo["DeltaFromFullCurrent"].abs()
)

loo = loo.sort_values(
    "AbsoluteInfluence",
    ascending=False,
).reset_index(drop=True)

loo.to_csv(
    LOO_OUTPUT,
    index=False,
)


lowest_pearson = loo.loc[
    loo["CurrentPartialPearson"].idxmin()
]

highest_pearson = loo.loc[
    loo["CurrentPartialPearson"].idxmax()
]

lowest_spearman = loo.loc[
    loo["CurrentResidualSpearman"].idxmin()
]

highest_spearman = loo.loc[
    loo["CurrentResidualSpearman"].idxmax()
]


# Nonparametric model-level bootstrap.
rng = np.random.default_rng(SEED)

boot = []

for _ in range(BOOTSTRAPS):
    idx = rng.integers(
        0,
        len(data),
        len(data),
    )

    sample = (
        data.iloc[idx]
        .reset_index(drop=True)
    )

    r = partial_r_only(sample)

    if np.isfinite(r):
        boot.append(r)

boot = np.asarray(
    boot,
    dtype=float,
)

if len(boot) == 0:
    raise RuntimeError(
        "No valid bootstrap samples"
    )

low, median, high = np.percentile(
    boot,
    [
        2.5,
        50,
        97.5,
    ],
)

positive_fraction = np.mean(
    boot > 0
)


lines = []


def emit(text=""):
    lines.append(text)
    print(text)


emit("CURRENT-WGS FPBM BEYOND MSI: ROBUSTNESS")
emit("=" * 72)
emit(f"Models: {len(data)}")
emit(f"Bootstrap samples requested: {BOOTSTRAPS}")
emit(f"Bootstrap seed: {SEED}")
emit()
emit("Full 26-model result")
emit("-" * 72)
emit(
    f"Partial Pearson:   {full['pearson_r']:.6f}   "
    f"p={full['pearson_p']:.6g}"
)
emit(
    f"Residual Spearman: {full['spearman_r']:.6f}   "
    f"p={full['spearman_p']:.6g}"
)
emit()
emit("Published-FPBM reference")
emit("-" * 72)
emit(
    f"Partial Pearson:   "
    f"{published_full['pearson_r']:.6f}   "
    f"p={published_full['pearson_p']:.6g}"
)
emit(
    f"Residual Spearman: "
    f"{published_full['spearman_r']:.6f}   "
    f"p={published_full['spearman_p']:.6g}"
)
emit()
emit("Leave-one-out sensitivity")
emit("-" * 72)
emit(
    f"Partial Pearson range: "
    f"{loo['CurrentPartialPearson'].min():.6f} "
    f"to {loo['CurrentPartialPearson'].max():.6f}"
)
emit(
    f"Residual Spearman range: "
    f"{loo['CurrentResidualSpearman'].min():.6f} "
    f"to {loo['CurrentResidualSpearman'].max():.6f}"
)
emit(
    "Positive partial Pearson after removal: "
    f"{int((loo['CurrentPartialPearson'] > 0).sum())}"
    f"/{len(loo)}"
)
emit(
    "Nominal p < 0.05 after removal: "
    f"{int((loo['CurrentPartialPearsonP'] < 0.05).sum())}"
    f"/{len(loo)}"
)
emit()
emit(
    "Lowest Pearson after removing: "
    f"{lowest_pearson['RemovedModelID']} "
    f"{lowest_pearson['RemovedCellLine']} "
    f"({lowest_pearson['CurrentPartialPearson']:.6f})"
)
emit(
    "Highest Pearson after removing: "
    f"{highest_pearson['RemovedModelID']} "
    f"{highest_pearson['RemovedCellLine']} "
    f"({highest_pearson['CurrentPartialPearson']:.6f})"
)
emit(
    "Lowest Spearman after removing: "
    f"{lowest_spearman['RemovedModelID']} "
    f"{lowest_spearman['RemovedCellLine']} "
    f"({lowest_spearman['CurrentResidualSpearman']:.6f})"
)
emit(
    "Highest Spearman after removing: "
    f"{highest_spearman['RemovedModelID']} "
    f"{highest_spearman['RemovedCellLine']} "
    f"({highest_spearman['CurrentResidualSpearman']:.6f})"
)
emit()
emit("Most influential removals")
emit("-" * 72)

for _, row in loo.head(8).iterrows():
    emit(
        f"{row['RemovedModelID']} "
        f"{row['RemovedCellLine']}: "
        f"r_without={row['CurrentPartialPearson']:.6f} "
        f"delta={row['DeltaFromFullCurrent']:+.6f}"
    )

emit()
emit("Bootstrap partial Pearson")
emit("-" * 72)
emit(
    f"Valid samples: {len(boot)}"
)
emit(
    f"Median: {median:.6f}"
)
emit(
    f"Mean:   {boot.mean():.6f}"
)
emit(
    f"95% percentile CI: "
    f"{low:.6f} to {high:.6f}"
)
emit(
    f"Fraction > 0: {positive_fraction:.6f}"
)
emit()
emit(f"Leave-one-out table: {LOO_OUTPUT}")

OUTPUT.write_text(
    "\n".join(lines) + "\n"
)
