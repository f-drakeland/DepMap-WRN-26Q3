from pathlib import Path

import pandas as pd


SOURCE = Path("data/source")
DERIVED = Path("data/derived")

WORKBOOK = SOURCE / "WRN_2020_SourceData_Fig4.xlsx"
GENE_EFFECT = SOURCE / "CRISPRGeneEffect.csv"
INVENTORY = SOURCE / "public_wgs_inventory.csv"
PROFILES = SOURCE / "OmicsProfiles.csv"

OUTPUT = DERIVED / "current_wgs_26_manifest.csv"

WRN_COLUMN = "WRN (7486)"


def s3_to_https(uri):
    prefix = "s3://depmap-omics-ccle/"

    if not isinstance(uri, str) or not uri.startswith(prefix):
        raise ValueError(f"Unexpected S3 URI: {uri!r}")

    relative = uri[len(prefix):]

    return (
        "https://depmap-omics-ccle.s3.amazonaws.com/"
        + relative
    )


# ------------------------------------------------------------
# Published MSI cohort
# ------------------------------------------------------------

published = pd.read_excel(
    WORKBOOK,
    sheet_name="4d",
    engine="openpyxl",
)

published = (
    published[
        published["CCLE_MSI"] == "MSI"
    ]
    .rename(
        columns={
            "DepMap_ID": "ModelID",
            "FPBM_mean": "PublishedFPBM",
        }
    )
    [
        [
            "ModelID",
            "PublishedFPBM",
            "Disease",
        ]
    ]
)

assert len(published) == 31


# ------------------------------------------------------------
# Current 26Q3 WRN GeneEffect
# ------------------------------------------------------------

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
    "ModelID",
    "WRNGeneEffect",
]

current["WRNGeneEffect"] = pd.to_numeric(
    current["WRNGeneEffect"],
    errors="coerce",
)

current = current.dropna(
    subset=["WRNGeneEffect"]
)

target = published.merge(
    current,
    on="ModelID",
    how="inner",
    validate="one_to_one",
)

assert len(target) == 26


# ------------------------------------------------------------
# Public hg38 WGS inventory
# ------------------------------------------------------------

inventory = pd.read_csv(INVENTORY)

required = {
    "model_id",
    "sample_id",
    "profile_id",
    "cell_line_name",
    "oncotree_lineage",
    "oncotree_primary_disease",
    "datatype",
    "reference_genome",
    "cram_bam_size",
    "cram_bam_url",
    "crai_bai_url",
}

missing = required - set(inventory.columns)

if missing:
    raise RuntimeError(
        f"Inventory missing columns: {sorted(missing)}"
    )

wgs = inventory[
    inventory["datatype"]
    .astype(str)
    .str.lower()
    .eq("wgs")
    &
    inventory["reference_genome"]
    .astype(str)
    .str.lower()
    .eq("hg38")
    &
    inventory["model_id"].isin(
        target["ModelID"]
    )
].copy()

wgs = wgs.dropna(
    subset=[
        "sample_id",
        "cram_bam_url",
        "crai_bai_url",
    ]
)


# ------------------------------------------------------------
# Default WGS profiles for resolving any duplicate inventory rows
# ------------------------------------------------------------

profiles = pd.read_csv(PROFILES)

profile_wgs = profiles[
    profiles["DataType"]
    .astype(str)
    .str.lower()
    .eq("wgs")
].copy()

default_model = (
    profile_wgs["IsDefaultEntryForModel"]
    .astype(str)
    .str.strip()
    .str.lower()
    .eq("yes")
)

default_mc = (
    profile_wgs["IsDefaultEntryForMC"]
    .astype(str)
    .str.strip()
    .str.lower()
    .eq("yes")
)

profile_defaults = profile_wgs[
    default_model & default_mc
]


# ------------------------------------------------------------
# Select exactly one public hg38 WGS row per model
# ------------------------------------------------------------

selected = []

for model_id in target["ModelID"]:
    candidates = wgs[
        wgs["model_id"] == model_id
    ].copy()

    if len(candidates) == 0:
        raise RuntimeError(
            f"No public hg38 WGS BAM/BAI for {model_id}"
        )

    if len(candidates) == 1:
        chosen = candidates.iloc[0]

    else:
        default_ids = set(
            profile_defaults.loc[
                profile_defaults["ModelID"] == model_id,
                "SequencingID",
            ].dropna()
        )

        matched = candidates[
            candidates["sample_id"].isin(
                default_ids
            )
        ]

        if len(matched) != 1:
            print()
            print(
                f"Ambiguous WGS inventory rows for {model_id}"
            )
            print(
                candidates[
                    [
                        "sample_id",
                        "profile_id",
                        "cram_bam_size",
                        "cram_bam_url",
                        "crai_bai_url",
                    ]
                ].to_string(index=False)
            )

            raise RuntimeError(
                f"Could not uniquely resolve {model_id}"
            )

        chosen = matched.iloc[0]

    selected.append(
        {
            "ModelID": model_id,
            "CellLineName": chosen["cell_line_name"],
            "Lineage": chosen["oncotree_lineage"],
            "SequencingID": chosen["sample_id"],
            "ProfileID": chosen["profile_id"],
            "BAMSizeBytes": int(chosen["cram_bam_size"]),
            "BAMSizeGiB": (
                int(chosen["cram_bam_size"])
                / 1024**3
            ),
            "BAMS3": chosen["cram_bam_url"],
            "BAIS3": chosen["crai_bai_url"],
            "BAMUrl": s3_to_https(
                chosen["cram_bam_url"]
            ),
            "BAIUrl": s3_to_https(
                chosen["crai_bai_url"]
            ),
        }
    )


manifest = target.merge(
    pd.DataFrame(selected),
    on="ModelID",
    how="inner",
    validate="one_to_one",
)

assert len(manifest) == 26
assert manifest["ModelID"].is_unique
assert manifest["SequencingID"].notna().all()
assert manifest["BAMUrl"].notna().all()
assert manifest["BAIUrl"].notna().all()

# Known anchor from the completed KM12 reconstruction.
km12 = manifest[
    manifest["ModelID"] == "ACH-000969"
]

assert len(km12) == 1
assert (
    km12.iloc[0]["SequencingID"]
    == "CDS-vHArJs"
)

OUTPUT.parent.mkdir(
    parents=True,
    exist_ok=True,
)

manifest.to_csv(
    OUTPUT,
    index=False,
)


print("CURRENT-WGS 26-MODEL MANIFEST")
print("=" * 72)
print("Published MSI models:", len(published))
print("Current 26Q3 overlap:", len(target))
print("Public hg38 WGS mapped:", len(manifest))

print()
print("KM12 anchor")
print("=" * 72)

row = km12.iloc[0]

print("ModelID:      ", row["ModelID"])
print("Cell line:    ", row["CellLineName"])
print("SequencingID: ", row["SequencingID"])
print(
    "BAM size:     "
    f"{row['BAMSizeGiB']:.1f} GiB"
)

print()
print("Manifest")
print("=" * 72)

print(
    manifest[
        [
            "ModelID",
            "CellLineName",
            "SequencingID",
            "BAMSizeGiB",
        ]
    ].to_string(
        index=False,
        formatters={
            "BAMSizeGiB": lambda x: f"{x:.1f}",
        },
    )
)

print()
print("Wrote:", OUTPUT)
