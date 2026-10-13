from pathlib import Path

import pandas as pd


SOURCE = Path("data/source")
OUTPUT = Path("data/derived/wrn_msi_26q3.csv")

GENE_EFFECT = SOURCE / "CRISPRGeneEffect.csv"
SUBTYPES = SOURCE / "OmicsInferredMolecularSubtypes.csv"
MODELS = SOURCE / "Model.csv"

WRN_COLUMN = "WRN (7486)"


def parse_msi(value):
    value = str(value).strip().lower()

    if value == "true":
        return True
    if value == "false":
        return False

    return pd.NA


# WRN dependency
gene_effect = pd.read_csv(GENE_EFFECT)

id_column = (
    "ModelID"
    if "ModelID" in gene_effect.columns
    else gene_effect.columns[0]
)

wrn = gene_effect[
    [id_column, WRN_COLUMN]
].copy()

wrn = wrn.rename(
    columns={id_column: "ModelID"}
)

wrn["WRNGeneEffect"] = pd.to_numeric(
    wrn[WRN_COLUMN],
    errors="coerce",
)

wrn = wrn[
    ["ModelID", "WRNGeneEffect"]
].dropna(subset=["WRNGeneEffect"])


# MSI status
subtypes = pd.read_csv(SUBTYPES)

msi = subtypes[
    ["ModelID", "MSI"]
].copy()

msi["MSI"] = msi["MSI"].map(parse_msi)
msi = msi.dropna(subset=["MSI"])
msi["MSI"] = msi["MSI"].astype(bool)


# Model metadata
models = pd.read_csv(MODELS)

metadata = models[
    [
        "ModelID",
        "OncotreeLineage",
        "OncotreePrimaryDisease",
    ]
].copy()

metadata = metadata.rename(
    columns={
        "OncotreeLineage": "Lineage",
        "OncotreePrimaryDisease": "Disease",
    }
)


# Preserve CRISPR GeneEffect row ordering.
result = (
    wrn
    .merge(msi, on="ModelID", how="inner", sort=False)
    .merge(metadata, on="ModelID", how="inner", sort=False)
)

result = result[
    [
        "ModelID",
        "WRNGeneEffect",
        "MSI",
        "Lineage",
        "Disease",
    ]
]


print("WRN / MSI 26Q3 TABLE")
print("=" * 60)
print("Models:", len(result))
print("MSI:", int(result["MSI"].sum()))
print("Non-MSI:", int((~result["MSI"]).sum()))


# Verify against the exploratory artifact before replacing anything.
if OUTPUT.exists():
    existing = pd.read_csv(OUTPUT)

    same = existing.reset_index(drop=True).equals(
        result.reset_index(drop=True)
    )

    print()
    print("Existing derived table:")
    print(f"Rows:        {len(existing)}")
    print(f"Exact match: {same}")

    if not same:
        comparison = Path(
            "data/derived/wrn_msi_26q3.rebuilt.csv"
        )
        result.to_csv(comparison, index=False)

        raise RuntimeError(
            "Reconstruction differs from existing artifact. "
            f"Wrote comparison file to {comparison}"
        )
else:
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(OUTPUT, index=False)
    print(f"\nWrote: {OUTPUT}")
