# Analysis Log

## 2026-10-12

Research workflow reconstructed from exploratory analysis.

### Published repeat-locus reproducibility

2,800 loci tested across discovery and validation models.

- Discovery-validation Pearson: 0.080
- Discovery-validation Spearman: 0.078
- Same correlation direction: 44.4%

Conclusion:
Individual repeat-locus correlations did not reproduce.

### Generic repeat burden beyond MSI

Mean repeat shift:
- Pearson: -0.209
- Spearman: -0.157
- Permutation p: 0.2152

Fraction shifted >=2:
- Pearson: -0.164
- Spearman: -0.160
- Permutation p: 0.3132

Conclusion:
Generic repeat burden did not explain WRN dependency beyond MSI.

### Published broken-TA phenotype

2020 source-data workbook recovered and Figure 4 source data inspected.

Published FPBM versus WRN dependency:

MSI:
- N = 31
- Pearson = 0.710
- Spearman = 0.795

MSS:
- N = 207
- Pearson = 0.104
- Spearman = 0.071

Conclusion:
The published broken-TA phenotype is strongly WRN-associated
specifically within MSI models.

### Comparison with DepMap 26Q3 WRN dependency

204 published models overlap current WRN GeneEffect data.

MSI:
- N = 26
- FPBM vs 26Q3 WRN Pearson = 0.651
- Spearman = 0.687

MSS:
- N = 178
- Pearson = 0.064
- Spearman = 0.045

Conclusion:
The published broken-TA phenotype remains associated with
current WRN dependency.

### Broken-TA phenotype beyond MSI

Among 26 MSI models:

- Partial Pearson = 0.395
- p = 0.04593
- Permutation p = 0.04582
- Residual Spearman = 0.305
- p = 0.1294

Robustness:
- Leave-one-out Pearson range = 0.187 to 0.514
- Bootstrap 95% CI = -0.184 to 0.734
- ACH-000969 materially influences the effect.

LOOCV:
- MSI-only RMSE = 0.7691
- MSI + FPBM RMSE = 0.7505
- Improvement = 2.43%
- Paired permutation p = 0.38300

Without ACH-000969:
- MSI-only RMSE = 0.7041
- MSI + FPBM RMSE = 0.7320

Conclusion:
FPBM is mechanistically informative but does not yet demonstrate
robust incremental predictive value beyond MSI.

### Current-WGS reconstruction

- Published broken-TA loci: 5,400 hg19 intervals
- Successfully lifted to hg38: 5,377
- Unmapped: 23
- Public hg38 WGS available: 26/26 MSI models
- Remote indexed BAM access verified with samtools
- KM12 (ACH-000969) selected as first full reconstruction target

Current operation:
Reconstructing broken/non-broken TA coverage for KM12 directly
from current public CCLE WGS without downloading the full BAM.

### KM12 current-WGS reconstruction

Current public hg38 WGS was queried across 5,377 lifted broken-TA
loci and 60,985 lifted non-broken TA loci.

KM12 (ACH-000969):

- Current broken-TA mean FPBM: 5.1926
- Published FPBM_mean: 4.8482
- Current / published ratio: 1.071

The current-WGS reconstruction is therefore within approximately
7.1% of the published KM12 phenotype despite differences in WGS
dataset, alignment/reference build, and coordinate conversion.

Conclusion:
The broken-TA FPBM reconstruction pipeline is sufficiently
validated in KM12 to justify scaling, pending access-performance
optimization.

### Executable provenance

The published-FPBM analytical branch is reproduced by:

- `scripts/reproduce_published_fpbm_wrn.py`
- `scripts/compare_published_fpbm_to_26q3.py`
- `scripts/test_fpbm_beyond_msi.py`
- `scripts/test_fpbm_beyond_msi_robustness.py`
- `scripts/test_fpbm_loocv_value.py`

All scripts reproduce the exploratory numerical results recorded above.
