from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr


SOURCE = Path("data/source")
DERIVED = Path("data/derived")

mapping = pd.read_csv(
    DERIVED / "repeat_wgs_49_models.csv"
).set_index("ModelID")

repeats = pd.read_csv(
    SOURCE / "OmicsMicrosatelliteRepeats.csv"
)

signatures = pd.read_csv(
    SOURCE / "OmicsGlobalSignatures.csv"
)

# ------------------------------------------------------------
# Cohort structure
# ------------------------------------------------------------

assert len(mapping) == 49
assert mapping.index.is_unique

expected_lineages = {
    "Bowel": 16,
    "Uterus": 19,
    "Lymphoid": 14,
}

actual_lineages = (
    mapping["Lineage"]
    .value_counts()
    .to_dict()
)

assert actual_lineages == expected_lineages

missing_seq = [
    seq
    for seq in mapping["SequencingID"]
    if seq not in repeats.columns
]

assert not missing_seq, missing_seq

# ------------------------------------------------------------
# Reconstruct repeat disruption:
#
# abs(observed repeat count - reference repeat count)
#
# Historical per-model generic burden used all 2,828 loci and
# ignored missing calls when taking the mean.
# ------------------------------------------------------------

reference = (
    repeats["Repeat"]
    .astype(str)
    .str.extract(r"^(\d+)\[")[0]
    .astype(float)
    .to_numpy()
)

assert len(reference) == 2828
assert np.isfinite(reference).all()

shift = pd.DataFrame(
    index=mapping.index,
    columns=np.arange(len(repeats)),
    dtype=float,
)

for model_id, row in mapping.iterrows():
    observed = pd.to_numeric(
        repeats[row["SequencingID"]],
        errors="coerce",
    ).to_numpy(float)

    shift.loc[model_id] = np.abs(
        observed - reference
    )

burden = shift.mean(
    axis=1,
    skipna=True,
)

# Larger value = stronger WRN dependency.
dependency = (
    -mapping["WRNGeneEffect"].astype(float)
)

# ------------------------------------------------------------
# Independent historical checksum:
# raw within-lineage burden / WRN correlations.
# ------------------------------------------------------------

lineage_targets = {
    "Bowel": 0.433,
    "Uterus": 0.397,
    "Lymphoid": -0.259,
}

lineage_results = {}

for lineage, target in lineage_targets.items():
    ids = mapping.index[
        mapping["Lineage"] == lineage
    ]

    r = pearsonr(
        burden.loc[ids],
        dependency.loc[ids],
    ).statistic

    lineage_results[lineage] = r

    assert round(r, 3) == target, (
        lineage,
        r,
        target,
    )

# ------------------------------------------------------------
# MSI provenance checksum.
#
# For these 49 models, the default-model MSI score and the
# score attached to the selected repeat SequencingID are
# identical.
# ------------------------------------------------------------

signatures["MSIScore"] = pd.to_numeric(
    signatures["MSIScore"],
    errors="coerce",
)

default_flag = (
    signatures["IsDefaultEntryForModel"]
    .astype(str)
    .str.strip()
    .str.lower()
)

default_msi = (
    signatures[
        default_flag.isin(["yes", "true", "1"])
    ]
    .dropna(subset=["ModelID", "MSIScore"])
    .groupby("ModelID")["MSIScore"]
    .first()
    .reindex(mapping.index)
)

seq_lookup = (
    signatures
    .dropna(subset=["SequencingID", "MSIScore"])
    .drop_duplicates("SequencingID")
    .set_index("SequencingID")["MSIScore"]
)

same_seq_msi = pd.Series(
    index=mapping.index,
    dtype=float,
)

for model_id, row in mapping.iterrows():
    same_seq_msi.loc[model_id] = seq_lookup.loc[
        row["SequencingID"]
    ]

assert default_msi.notna().all()
assert same_seq_msi.notna().all()

assert np.allclose(
    default_msi.to_numpy(float),
    same_seq_msi.to_numpy(float),
)

repeat_msi_p = pearsonr(
    burden,
    default_msi,
).statistic

repeat_msi_s = spearmanr(
    burden,
    default_msi,
).statistic

msi_wrn_p = pearsonr(
    default_msi,
    dependency,
).statistic

msi_wrn_s = spearmanr(
    default_msi,
    dependency,
).statistic

assert round(repeat_msi_p, 3) == 0.946
assert round(repeat_msi_s, 3) == 0.990
assert round(msi_wrn_p, 3) == 0.621
assert round(msi_wrn_s, 3) == 0.724

# ------------------------------------------------------------
# Report
# ------------------------------------------------------------

print("Repeat-WGS reconstruction validated")
print("=" * 64)
print(f"Models: {len(mapping)}")
print(f"Repeat loci: {len(repeats)}")

print()
print("Lineages")
for lineage in ("Bowel", "Uterus", "Lymphoid"):
    print(
        f"  {lineage:<10} "
        f"{actual_lineages[lineage]:>2}"
    )

print()
print("Historical burden / WRN checksum")
for lineage in ("Bowel", "Uterus", "Lymphoid"):
    print(
        f"  {lineage:<10} "
        f"{lineage_results[lineage]:+.3f}"
    )

print()
print("Historical MSI checksum")
print(
    "  Repeat burden vs MSI: "
    f"Pearson {repeat_msi_p:.3f}, "
    f"Spearman {repeat_msi_s:.3f}"
)
print(
    "  MSI vs WRN dependency: "
    f"Pearson {msi_wrn_p:.3f}, "
    f"Spearman {msi_wrn_s:.3f}"
)

print()
print(
    "Default-model MSI == selected-repeat-profile MSI: yes"
)
