from pathlib import Path

import pandas as pd


SOURCE = Path("data/source")
DERIVED = Path("data/derived")

TARGET_LINEAGES = {
    "Bowel",
    "Uterus",
    "Lymphoid",
}

profiles = pd.read_csv(SOURCE / "OmicsProfiles.csv")
repeats = pd.read_csv(SOURCE / "OmicsMicrosatelliteRepeats.csv")
cohort = pd.read_csv(DERIVED / "wrn_msi_26q3.csv")

repeat_ids = set(repeats.columns[6:])

target = cohort[
    cohort["MSI"].astype(str).str.lower().eq("true")
    & cohort["Lineage"].isin(TARGET_LINEAGES)
].copy()

assert len(target) == 49

wgs = profiles[
    profiles["DataType"].astype(str).str.lower().eq("wgs")
    & profiles["SequencingID"].isin(repeat_ids)
    & profiles["ModelID"].isin(target["ModelID"])
].copy()

selected = []

for model_id in target["ModelID"]:
    candidates = wgs[
        wgs["ModelID"] == model_id
    ].copy()

    if len(candidates) == 0:
        raise RuntimeError(
            f"No repeat-matrix WGS profile for {model_id}"
        )

    if len(candidates) == 1:
        chosen = candidates.iloc[0]

    else:
        default = candidates[
            candidates["IsDefaultEntryForModel"]
            .astype(str)
            .str.lower()
            .eq("yes")
            &
            candidates["IsDefaultEntryForMC"]
            .astype(str)
            .str.lower()
            .eq("yes")
        ]

        if len(default) != 1:
            raise RuntimeError(
                f"Ambiguous default profile for {model_id}"
            )

        chosen = default.iloc[0]

    selected.append(
        {
            "ModelID": model_id,
            "SequencingID": chosen["SequencingID"],
        }
    )

mapping = (
    target[
        [
            "ModelID",
            "WRNGeneEffect",
            "Lineage",
            "Disease",
        ]
    ]
    .merge(
        pd.DataFrame(selected),
        on="ModelID",
        how="left",
        validate="one_to_one",
    )
)

mapping_path = DERIVED / "repeat_wgs_49_models.csv"
mapping.to_csv(mapping_path, index=False)

seq_ids = mapping["SequencingID"].tolist()

matrix = repeats[seq_ids].apply(
    pd.to_numeric,
    errors="coerce",
)

observed_per_locus = matrix.notna().sum(axis=1)

print("CANONICAL 49-MODEL REPEAT COHORT")
print("=" * 70)
print("Models:", len(mapping))
print()
print(mapping["Lineage"].value_counts().to_string())

print()
print("Duplicate resolutions")
print("=" * 70)

for model_id in ["ACH-001669", "ACH-001709"]:
    row = mapping[
        mapping["ModelID"] == model_id
    ].iloc[0]

    print(
        f"{model_id} -> {row['SequencingID']}"
    )

print()
print("LOCUS COMPLETENESS")
print("=" * 70)
print("Total repeat loci:", len(repeats))
print(
    "Complete across all 49 models:",
    int((observed_per_locus == 49).sum()),
)
print(
    "At least 48/49:",
    int((observed_per_locus >= 48).sum()),
)
print(
    "At least 47/49:",
    int((observed_per_locus >= 47).sum()),
)
print(
    "Incomplete loci:",
    int((observed_per_locus < 49).sum()),
)

print()
print("Observed-model count distribution")
print("=" * 70)

print(
    observed_per_locus
    .value_counts()
    .sort_index(ascending=False)
    .head(15)
    .to_string()
)

print()
print("Wrote:", mapping_path)
