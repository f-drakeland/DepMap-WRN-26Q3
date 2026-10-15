# Data

This repository does not version large downloaded DepMap datasets, public WGS
alignments, third-party source workbooks, external reference chains, or
generated BED coordinate files.

## source/

Downloaded source datasets used by the analyses include:

- DepMap 26Q3 CRISPR GeneEffect
- DepMap molecular/signature tables
- DepMap omics profile metadata
- public CCLE WGS inventory
- 2020 WRN study Figure 4 source-data workbook

These files are intentionally excluded from Git.

## reference/

External reference assets include:

- UCSC `hg19ToHg38.over.chain.gz`

Reference assets are intentionally excluded from Git.

## derived/

Small derived manifests and provenance records may be committed.

Generated BED coordinate files are excluded because they are reproducibly
generated from the source workbook and liftOver chain.

The canonical coordinate build is:

```bash
python3 scripts/build_current_wgs_coordinates.py
```

That script:

1. reads sheet `4e` of the 2020 Figure 4 source workbook;
2. converts source coordinates from 1-based inclusive to 0-based half-open
   BED coordinates;
3. runs UCSC liftOver from hg19 to hg38 for broken and non-broken TA loci;
4. verifies mapped and unmapped locus counts;
5. concatenates the mapped BED files;
6. reproduces the historical natural/version sort with
   `LC_ALL=C sort -V -k1,1 -k2,2n`;
7. verifies the final 66,362-locus canonical BED SHA256.

Expected coordinate counts:

```text
Broken hg19              5,400
Broken hg38 mapped       5,377
Broken unmapped             23

Non-broken hg19         61,244
Non-broken hg38 mapped  60,985
Non-broken unmapped        259

Canonical hg38 total    66,362
```

Expected canonical SHA256:

```text
e0cb25a38b4867c39b2f23f6222dca6cdf253489287ff6ff7ddf56d3369fc24b
```

Tracked build provenance:

```text
data/derived/coordinates/build_provenance.json
```

The external UCSC liftOver executable is expected at:

```text
tools/liftOver
```

The liftOver binary, chain file, source workbook, and generated BED files are
not committed.
