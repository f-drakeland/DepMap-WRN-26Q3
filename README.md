# DepMap WRN 26Q3

Reproducibility and extension study of the broken-(TA)n repeat phenotype
associated with WRN dependency in microsatellite-instable cancer models.

## Research question

Can the broken-(TA)n repeat phenotype associated with WRN dependency be
reconstructed from current public CCLE whole-genome sequencing, and does it
explain variation in current WRN dependency among MSI models beyond
conventional MSI measurement?

## Current conclusion

The broken-TA phenotype is highly reproducible from current public WGS and
retains its expected relationship with current WRN dependency.

Across 26 MSI models:

- published vs reconstructed FPBM:
  - Pearson `r = 0.969`
  - Spearman `r = 0.973`
- reconstructed FPBM vs current WRN GeneEffect:
  - Pearson `r = 0.643`
  - Spearman `r = 0.621`
- after adjustment for current DepMap MSI score:
  - partial Pearson `r = 0.402`
  - parametric `p = 0.0466`
  - residual permutation `p = 0.0430`
  - residual Spearman `r = 0.301`, `p = 0.135`

The MSI-adjusted association is modest and uncertain. Its direction remains
positive after removal of every individual model, but the bootstrap 95%
interval spans zero and KM12 has substantial leverage.

Out-of-sample prediction does not support a robust incremental prediction
claim. Adding reconstructed FPBM to MSI improves all-26 LOOCV RMSE by 3.7%,
but the improvement reverses when KM12 is removed.

The current evidence therefore supports FPBM as a reproducible and
mechanistically informative phenotype, not as a demonstrated robust predictor
beyond MSI.

## Coordinate reconstruction

The 2020 Figure 4 source workbook contains:

- 5,400 broken-TA hg19 loci
- 61,244 non-broken-TA hg19 loci

UCSC liftOver to hg38 retains:

- 5,377 broken loci; 23 unmapped
- 60,985 non-broken loci; 259 unmapped
- 66,362 total mapped loci

The canonical coordinate file reproduces the recovered historical sort:

```text
LC_ALL=C sort -V -k1,1 -k2,2n
```

Canonical BED SHA256:

```text
e0cb25a38b4867c39b2f23f6222dca6cdf253489287ff6ff7ddf56d3369fc24b
```

The complete coordinate build is executable with:

```bash
python3 scripts/build_current_wgs_coordinates.py
```

## Current-WGS reconstruction

Public hg38 WGS is available for all 26 MSI models in the current comparison
cohort.

The extractor uses indexed remote BAM access rather than downloading complete
WGS BAMs, which are roughly 60–106 GiB each. Chromosomes are fetched in
parallel, coverage is calculated locally with `samtools bedcov`, completed
chromosome work is resumable, transient remote failures are retried, and final
coverage files are checksum-recorded.

The reconstructed phenotype strongly reproduces the published FPBM:

```text
Pearson r   0.968521
Spearman r  0.973333
R²          0.938033
```

A calibration audit shows that 25 of 26 models follow an especially tight
approximately multiplicative scale shift:

```text
Published/current Pearson   0.996253
R²                          0.992521
Median current/published    1.086431
Through-origin scale        1.085697
24/25 within ±2% of median scale
25/25 within ±3% of median scale
```

JHUEM-1 is the single prominent reconstruction discrepancy. Its
current/published ratio is `1.225260`, compared with `1.066891–1.109085`
among the other 25 models.

The discrepancy is not explained by an alternate public WGS candidate, BAM
size, remote-fetch retries, coordinate version, or extractor version. It is
retained as an unresolved sample/model-specific limitation rather than
corrected post hoc.

## Canonical analysis chain

```bash
python3 scripts/analyze_current_wgs_reconstruction.py
python3 scripts/analyze_current_wgs_wrn_association.py
python3 scripts/analyze_current_wgs_beyond_msi.py
python3 scripts/analyze_current_wgs_beyond_msi_robustness.py
python3 scripts/analyze_current_wgs_loocv_value.py
python3 scripts/audit_current_wgs_calibration.py
```

Generated reports are under:

```text
results/current-wgs/
```

The detailed scientific history, including exploratory analyses whose exact
downstream implementations were not retained, is documented in
`docs/ANALYSIS_LOG.md`.

## Reconstructing the WGS cohort

Build the current WGS manifest:

```bash
python3 scripts/build_current_wgs_manifest.py
```

Reconstruct one model:

```bash
python3 scripts/extract_current_wgs_ta.py ACH-000969
```

Run the remaining cohort:

```bash
python3 scripts/run_current_wgs_cohort.py
```

The extractor requires `samtools` and network access to the public CCLE WGS
BAM/BAI objects. Full WGS BAM files are not downloaded.

## External data

Large and third-party source datasets are intentionally not committed.

Required sources include:

- 2020 WRN Figure 4 source-data workbook
- DepMap 26Q3 CRISPR GeneEffect
- DepMap molecular/signature tables
- public CCLE WGS inventory and profile metadata
- UCSC hg19-to-hg38 liftOver chain

See `data/README.md`.

Python dependencies are listed in `requirements.txt`.

External command-line requirements:

- `samtools`
- UCSC `liftOver`

## Limitations

This is a 26-model MSI comparison, not an independent prospective validation
cohort.

The published-FPBM/current-WRN comparison uses overlapping cell-line models
across different data releases and is not independent biological replication.

The reconstructed phenotype retains a modest MSI-independent linear
association, but its uncertainty is substantial and its apparent incremental
predictive value is not robust to removal of KM12.

JHUEM-1 remains an unexplained measurement discrepancy.

The project therefore distinguishes phenotype reproducibility, mechanistic
association, and predictive utility rather than treating them as the same
claim.
