from pathlib import Path

import pandas as pd


BED = Path(
    "data/derived/coordinates/all_ta_hg38.sorted.bed"
)

THRESHOLDS = [
    0,
    100,
    500,
    1_000,
    5_000,
    10_000,
    25_000,
    50_000,
    100_000,
    250_000,
    500_000,
    1_000_000,
]


bed = pd.read_csv(
    BED,
    sep="\t",
    header=None,
    names=["chrom", "start", "end", "locus"],
)

assert len(bed) == 66362


def merge_windows(df, max_gap):
    windows = []

    for chrom, group in df.groupby(
        "chrom",
        sort=False,
    ):
        group = group.sort_values(
            ["start", "end"]
        )

        current_start = None
        current_end = None
        loci = 0

        for row in group.itertuples():
            if current_start is None:
                current_start = row.start
                current_end = row.end
                loci = 1
                continue

            gap = row.start - current_end

            if gap <= max_gap:
                current_end = max(
                    current_end,
                    row.end,
                )
                loci += 1
            else:
                windows.append(
                    (
                        chrom,
                        current_start,
                        current_end,
                        loci,
                    )
                )

                current_start = row.start
                current_end = row.end
                loci = 1

        windows.append(
            (
                chrom,
                current_start,
                current_end,
                loci,
            )
        )

    return pd.DataFrame(
        windows,
        columns=[
            "chrom",
            "start",
            "end",
            "loci",
        ],
    )


original_bases = (
    bed["end"] - bed["start"]
).sum()

print("TA FETCH-WINDOW PLANNER")
print("=" * 86)
print(f"TA intervals: {len(bed):,}")
print(f"Original TA bases: {original_bases:,}")
print()

print(
    f"{'Gap':>10} "
    f"{'Windows':>10} "
    f"{'Reduction':>11} "
    f"{'Fetch Mb':>12} "
    f"{'Expansion':>11} "
    f"{'Median kb':>11} "
    f"{'P95 kb':>10}"
)

print("-" * 86)

for gap in THRESHOLDS:
    windows = merge_windows(
        bed,
        gap,
    )

    lengths = (
        windows["end"]
        - windows["start"]
    )

    total_bases = lengths.sum()

    reduction = (
        1
        - len(windows) / len(bed)
    ) * 100

    expansion = (
        total_bases / original_bases
    )

    print(
        f"{gap:>10,} "
        f"{len(windows):>10,} "
        f"{reduction:>10.1f}% "
        f"{total_bases / 1e6:>11.1f} "
        f"{expansion:>10.1f}x "
        f"{lengths.median() / 1000:>10.1f} "
        f"{lengths.quantile(.95) / 1000:>9.1f}"
    )
