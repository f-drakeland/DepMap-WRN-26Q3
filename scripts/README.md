# Scripts

The repository contains both the canonical current-WGS analysis pipeline and
earlier exploratory or benchmarking code preserved for provenance.

For the final scientific conclusions, use the canonical scripts below.

## Canonical current-WGS pipeline

### Source and coordinate verification

- `verify_source_manifest.py`
  - verifies local source files against `data/SOURCE_MANIFEST.tsv`
- `build_current_wgs_coordinates.py`
  - reconstructs the published broken/non-broken TA coordinates in hg38
  - verifies expected liftOver counts and the canonical BED SHA256
- `build_current_wgs_manifest.py`
  - identifies the 26-model current public-WGS comparison cohort

### WGS reconstruction

- `extract_current_wgs_ta.py`
  - performs indexed remote BAM extraction
  - computes locus coverage
  - writes per-model provenance and FPBM results
- `run_current_wgs_cohort.py`
  - runs/resumes reconstruction across the cohort
  - writes the cohort summary

### Final scientific analyses

- `analyze_current_wgs_reconstruction.py`
  - published FPBM versus reconstructed current-WGS FPBM
- `analyze_current_wgs_wrn_association.py`
  - reconstructed FPBM versus current WRN GeneEffect
- `analyze_current_wgs_beyond_msi.py`
  - MSI-adjusted FPBM/WRN association using the canonical DepMap default-profile
    selection rule
- `analyze_current_wgs_beyond_msi_robustness.py`
  - leave-one-out and bootstrap robustness
- `analyze_current_wgs_loocv_value.py`
  - out-of-sample incremental predictive-value test
- `audit_current_wgs_calibration.py`
  - cohort-wide FPBM calibration and JHUEM-1 discrepancy audit
- `make_current_wgs_figures.py`
  - deterministically regenerates the three presentation figures

## Historical and recovered analyses

These scripts record earlier project stages and remain useful as provenance,
but they are not the canonical final current-WGS analysis path:

- `build_wrn_msi_table.py`
- `compare_published_fpbm_to_26q3.py`
- `reproduce_published_fpbm_wrn.py`
- `test_fpbm_beyond_msi.py`
- `test_fpbm_beyond_msi_robustness.py`
- `test_fpbm_loocv_value.py`
- `build_repeat_wgs_cohort.py`
- `validate_repeat_wgs_reconstruction.py`
- `build_ta_coordinates.py`

In particular, the older FPBM-beyond-MSI scripts predate recovery of the
correct `IsDefaultEntryForModel == "Yes"` MSI-profile selection rule. Their
historical results are retained but do not supersede the canonical
`analyze_current_wgs_*` results.

## WGS acquisition development and benchmarks

These scripts document how the final remote-WGS extraction strategy was
developed:

- `plan_ta_fetch_windows.py`
- `analyze_km12_fpbm.py`
- `benchmark_km12_multiregion.py`
- `benchmark_km12_chrom_parallel.py`

They established the indexed, chromosome-parallel extraction method later
implemented in `extract_current_wgs_ta.py`.

No additional scientific conclusion should be inferred from benchmark scripts
unless it is also represented in the canonical analysis outputs.
