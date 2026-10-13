from pathlib import Path
import openpyxl

SOURCE = Path("data/source/WRN_2020_SourceData_Fig4.xlsx")
OUT = Path("data/derived/coordinates")

OUT.mkdir(parents=True, exist_ok=True)

wb = openpyxl.load_workbook(
    SOURCE,
    read_only=True,
    data_only=True,
)
ws = wb["4e"]

broken_path = OUT / "broken_ta_hg19.bed"
nonbroken_path = OUT / "nonbroken_ta_hg19.bed"

broken = 0
nonbroken = 0

with broken_path.open("w") as b, nonbroken_path.open("w") as nb:
    for row in ws.iter_rows(min_row=5, values_only=True):

        # Broken TA loci: columns A-C
        chrom, start, stop = row[0:3]

        if chrom is not None and start is not None and stop is not None:
            broken += 1

            # Spreadsheet coordinates are 1-based inclusive.
            # BED uses 0-based, half-open coordinates.
            b.write(
                f"{chrom}\t{int(start)-1}\t{int(stop)}"
                f"\tbrokenTA_{broken:04d}\n"
            )

        # Non-broken TA loci: columns J-L
        chrom, start, stop = row[9:12]

        if chrom is not None and start is not None and stop is not None:
            nonbroken += 1
            nb.write(
                f"{chrom}\t{int(start)-1}\t{int(stop)}"
                f"\tNB_{nonbroken:05d}\n"
            )

print(f"Broken hg19 loci:     {broken}")
print(f"Non-broken hg19 loci: {nonbroken}")
print(f"Wrote: {broken_path}")
print(f"Wrote: {nonbroken_path}")
