from pathlib import Path
import csv
import hashlib
import sys


MANIFEST = Path("data/SOURCE_MANIFEST.tsv")
SOURCE = Path("data/source")


def sha256(path):
    h = hashlib.sha256()

    with path.open("rb") as f:
        for chunk in iter(
            lambda: f.read(1024 * 1024),
            b"",
        ):
            h.update(chunk)

    return h.hexdigest()


with MANIFEST.open(newline="") as f:
    rows = list(
        csv.DictReader(
            f,
            delimiter="\t",
        )
    )

if not rows:
    raise RuntimeError(
        "Source manifest is empty"
    )

failures = []

print("SOURCE MANIFEST VERIFICATION")
print("=" * 72)

for row in rows:
    path = SOURCE / row["filename"]

    if not path.exists():
        failures.append(
            f"MISSING: {path}"
        )
        print(
            f"FAIL  {row['filename']}: missing"
        )
        continue

    expected_size = int(
        row["size_bytes"]
    )

    actual_size = path.stat().st_size

    if actual_size != expected_size:
        failures.append(
            f"SIZE: {path}: "
            f"expected {expected_size}, "
            f"found {actual_size}"
        )
        print(
            f"FAIL  {row['filename']}: "
            "size mismatch"
        )
        continue

    actual_sha = sha256(path)
    expected_sha = row["sha256"]

    if actual_sha != expected_sha:
        failures.append(
            f"SHA256: {path}: "
            f"expected {expected_sha}, "
            f"found {actual_sha}"
        )
        print(
            f"FAIL  {row['filename']}: "
            "SHA256 mismatch"
        )
        continue

    print(
        f"PASS  {row['filename']}"
    )

print()

if failures:
    print("VERIFICATION FAILED")
    print("=" * 72)

    for failure in failures:
        print(failure)

    sys.exit(1)

print(
    f"PASS: {len(rows)} source files "
    "match the tracked manifest."
)
