# Reproducibility Environment

Environment recorded at the scientific closeout checkpoint on 2026-10-15.

## Python

```text
Python 3.13.5
```

Exact project package versions are recorded in:

```text
requirements.lock.txt
```

The less restrictive `requirements.txt` is retained as the normal project
dependency list.

## External tools

```text
samtools 1.24
```

The current-WGS extractor uses `samtools` for indexed remote BAM access,
alignment validation/indexing, and coverage calculation.

The coordinate builder also requires the UCSC `liftOver` executable and the
`hg19ToHg38.over.chain.gz` reference chain.

Exact SHA256 values for the liftOver binary and chain used by the canonical
coordinate build are recorded in:

```text
data/derived/coordinates/build_provenance.json
```

## Source data

Third-party and large source files are not versioned in Git.

Exact filenames, byte sizes, and SHA256 values for the local source corpus are
recorded in:

```text
data/SOURCE_MANIFEST.tsv
```

Verify an assembled local source directory with:

```bash
python3 scripts/verify_source_manifest.py
```

The verifier fails if a required file is missing or differs in size or SHA256.
