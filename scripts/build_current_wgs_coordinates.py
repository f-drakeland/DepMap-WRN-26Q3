from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess
import tempfile

import openpyxl


SOURCE = Path(
    "data/source/WRN_2020_SourceData_Fig4.xlsx"
)

CHAIN = Path(
    "data/reference/hg19ToHg38.over.chain.gz"
)

LIFTOVER = Path(
    "tools/liftOver"
)

OUT = Path(
    "data/derived/coordinates"
)

PROVENANCE = (
    OUT / "build_provenance.json"
)

SCRIPT = Path(__file__)

EXPECTED = {
    "broken_hg19": 5400,
    "broken_hg38": 5377,
    "broken_unmapped": 23,
    "nonbroken_hg19": 61244,
    "nonbroken_hg38": 60985,
    "nonbroken_unmapped": 259,
    "canonical_hg38": 66362,
}

EXPECTED_CANONICAL_SHA256 = (
    "e0cb25a38b4867c39b2f23f6222dca6"
    "cdf253489287ff6ff7ddf56d3369fc24b"
)


def sha256(path):
    h = hashlib.sha256()

    with path.open("rb") as f:
        for chunk in iter(
            lambda: f.read(1024 * 1024),
            b"",
        ):
            h.update(chunk)

    return h.hexdigest()


def count_bed_records(path):
    count = 0

    with path.open() as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith("#"):
                continue

            count += 1

    return count


def require_file(path):
    if not path.exists():
        raise RuntimeError(
            f"Required file not found: {path}"
        )


def run(cmd, **kwargs):
    subprocess.run(
        [str(x) for x in cmd],
        check=True,
        **kwargs,
    )


for required in (
    SOURCE,
    CHAIN,
    LIFTOVER,
):
    require_file(required)

if not os.access(
    LIFTOVER,
    os.X_OK,
):
    raise RuntimeError(
        f"liftOver is not executable: {LIFTOVER}"
    )

OUT.mkdir(
    parents=True,
    exist_ok=True,
)

with tempfile.TemporaryDirectory(
    prefix="wrn-ta-coordinates-"
) as tmp_name:

    tmp = Path(tmp_name)

    broken_hg19 = (
        tmp / "broken_ta_hg19.bed"
    )

    nonbroken_hg19 = (
        tmp / "nonbroken_ta_hg19.bed"
    )

    broken_hg38 = (
        tmp / "broken_ta_hg38.bed"
    )

    broken_unmapped = (
        tmp / "broken_ta_unmapped.bed"
    )

    nonbroken_hg38 = (
        tmp / "nonbroken_ta_hg38.bed"
    )

    nonbroken_unmapped = (
        tmp / "nonbroken_ta_unmapped.bed"
    )

    combined = (
        tmp / "all_ta_hg38.bed"
    )

    canonical = (
        tmp / "all_ta_hg38.sorted.bed"
    )

    # --------------------------------------------------------
    # Source workbook -> hg19 BED
    # --------------------------------------------------------

    wb = openpyxl.load_workbook(
        SOURCE,
        read_only=True,
        data_only=True,
    )

    ws = wb["4e"]

    broken_count = 0
    nonbroken_count = 0

    with (
        broken_hg19.open("w") as broken_out,
        nonbroken_hg19.open("w") as nonbroken_out,
    ):
        for row in ws.iter_rows(
            min_row=5,
            values_only=True,
        ):

            chrom, start, stop = row[0:3]

            if (
                chrom is not None
                and start is not None
                and stop is not None
            ):
                broken_count += 1

                broken_out.write(
                    f"{chrom}\t"
                    f"{int(start) - 1}\t"
                    f"{int(stop)}\t"
                    f"brokenTA_{broken_count:04d}\n"
                )

            chrom, start, stop = row[9:12]

            if (
                chrom is not None
                and start is not None
                and stop is not None
            ):
                nonbroken_count += 1

                nonbroken_out.write(
                    f"{chrom}\t"
                    f"{int(start) - 1}\t"
                    f"{int(stop)}\t"
                    f"NB_{nonbroken_count:05d}\n"
                )

    wb.close()

    if broken_count != EXPECTED["broken_hg19"]:
        raise RuntimeError(
            f"Broken hg19 count: {broken_count}"
        )

    if (
        nonbroken_count
        != EXPECTED["nonbroken_hg19"]
    ):
        raise RuntimeError(
            f"Non-broken hg19 count: {nonbroken_count}"
        )

    # --------------------------------------------------------
    # hg19 -> hg38
    # --------------------------------------------------------

    run(
        [
            LIFTOVER,
            broken_hg19,
            CHAIN,
            broken_hg38,
            broken_unmapped,
        ]
    )

    run(
        [
            LIFTOVER,
            nonbroken_hg19,
            CHAIN,
            nonbroken_hg38,
            nonbroken_unmapped,
        ]
    )

    counts = {
        "broken_hg19":
            count_bed_records(
                broken_hg19
            ),
        "broken_hg38":
            count_bed_records(
                broken_hg38
            ),
        "broken_unmapped":
            count_bed_records(
                broken_unmapped
            ),
        "nonbroken_hg19":
            count_bed_records(
                nonbroken_hg19
            ),
        "nonbroken_hg38":
            count_bed_records(
                nonbroken_hg38
            ),
        "nonbroken_unmapped":
            count_bed_records(
                nonbroken_unmapped
            ),
    }

    for key, expected in EXPECTED.items():
        if key == "canonical_hg38":
            continue

        actual = counts[key]

        if actual != expected:
            raise RuntimeError(
                f"{key}: expected {expected}, "
                f"found {actual}"
            )

    # --------------------------------------------------------
    # Canonical merge and exact historical sort.
    # --------------------------------------------------------

    with combined.open("wb") as out:
        for source in (
            broken_hg38,
            nonbroken_hg38,
        ):
            with source.open("rb") as src:
                shutil.copyfileobj(
                    src,
                    out,
                )

    sort_env = os.environ.copy()
    sort_env["LC_ALL"] = "C"

    with canonical.open("wb") as out:
        run(
            [
                "sort",
                "-V",
                "-k1,1",
                "-k2,2n",
                combined,
            ],
            stdout=out,
            env=sort_env,
        )

    canonical_count = (
        count_bed_records(
            canonical
        )
    )

    if (
        canonical_count
        != EXPECTED["canonical_hg38"]
    ):
        raise RuntimeError(
            "Canonical hg38 count: "
            f"expected "
            f"{EXPECTED['canonical_hg38']}, "
            f"found {canonical_count}"
        )

    canonical_sha = sha256(
        canonical
    )

    if (
        canonical_sha
        != EXPECTED_CANONICAL_SHA256
    ):
        raise RuntimeError(
            "Canonical BED SHA256 mismatch:\n"
            f"expected: "
            f"{EXPECTED_CANONICAL_SHA256}\n"
            f"actual:   {canonical_sha}"
        )

    # --------------------------------------------------------
    # Validation passed. Publish generated files atomically
    # with respect to scientific validation.
    # --------------------------------------------------------

    generated = [
        broken_hg19,
        broken_hg38,
        broken_unmapped,
        nonbroken_hg19,
        nonbroken_hg38,
        nonbroken_unmapped,
        combined,
        canonical,
    ]

    for source in generated:
        shutil.copy2(
            source,
            OUT / source.name,
        )

# ------------------------------------------------------------
# Durable tracked provenance.
# ------------------------------------------------------------

provenance = {
    "schema_version": 1,
    "source_workbook": str(SOURCE),
    "source_sheet": "4e",
    "source_workbook_sha256":
        sha256(SOURCE),
    "liftover_chain": str(CHAIN),
    "liftover_chain_sha256":
        sha256(CHAIN),
    "liftover_binary": str(LIFTOVER),
    "liftover_binary_sha256":
        sha256(LIFTOVER),
    "builder_script": str(SCRIPT),
    "builder_script_sha256":
        sha256(SCRIPT),
    "coordinate_convention": (
        "source workbook interpreted as "
        "1-based inclusive; BED written as "
        "0-based half-open"
    ),
    "sort_command": (
        "LC_ALL=C sort -V -k1,1 -k2,2n"
    ),
    "counts": EXPECTED,
    "canonical_bed": str(
        OUT / "all_ta_hg38.sorted.bed"
    ),
    "canonical_bed_sha256":
        EXPECTED_CANONICAL_SHA256,
}

PROVENANCE.write_text(
    json.dumps(
        provenance,
        indent=2,
        sort_keys=True,
    )
    + "\n"
)

print("CURRENT-WGS TA COORDINATE BUILD")
print("=" * 72)
print(
    f"Broken hg19:           "
    f"{EXPECTED['broken_hg19']}"
)
print(
    f"Broken hg38 mapped:    "
    f"{EXPECTED['broken_hg38']}"
)
print(
    f"Broken unmapped:       "
    f"{EXPECTED['broken_unmapped']}"
)
print(
    f"Non-broken hg19:       "
    f"{EXPECTED['nonbroken_hg19']}"
)
print(
    f"Non-broken hg38 mapped:"
    f" {EXPECTED['nonbroken_hg38']}"
)
print(
    f"Non-broken unmapped:   "
    f"{EXPECTED['nonbroken_unmapped']}"
)
print(
    f"Canonical hg38 loci:   "
    f"{EXPECTED['canonical_hg38']}"
)
print()
print(
    "Canonical SHA256: "
    f"{EXPECTED_CANONICAL_SHA256}"
)
print(
    f"Provenance: {PROVENANCE}"
)
print()
print(
    "PASS: canonical coordinate build "
    "reproduced exactly."
)
