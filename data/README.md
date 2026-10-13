# Data

This repository does not version large downloaded DepMap datasets,
public WGS alignments, or third-party source workbooks.

## source/

Downloaded source datasets used by the analyses, including:

- DepMap 26Q3 CRISPR GeneEffect
- DepMap molecular/omics tables
- public CCLE WGS inventory
- 2020 WRN study Figure 4 source-data workbook

These files are intentionally excluded from Git.

## reference/

External reference assets such as the UCSC hg19-to-hg38 liftOver chain.
These are also excluded from Git.

## derived/

Small analysis-derived datasets and coordinate definitions.

The published broken-TA and non-broken-TA coordinates are generated
from the Figure 4 source workbook by:

    scripts/build_ta_coordinates.py

Coordinates are converted from hg19 to hg38 with UCSC liftOver.
Generated BED files are excluded from Git because they are reproducible.
