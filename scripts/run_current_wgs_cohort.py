from pathlib import Path
import argparse
import json
import subprocess
import sys

import pandas as pd


MANIFEST = Path(
    "data/derived/current_wgs_26_manifest.csv"
)

RESULT_ROOT = Path(
    "results/current-wgs"
)

EXTRACTOR = Path(
    "scripts/extract_current_wgs_ta.py"
)

SUMMARY = (
    RESULT_ROOT
    / "cohort_summary.csv"
)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Run current-WGS broken-TA reconstruction "
            "across the manifest cohort."
        )
    )

    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help=(
            "Chromosome workers passed to extractor "
            "(default: 4)"
        ),
    )

    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help=(
            "Maximum number of unfinished models to run "
            "during this invocation"
        ),
    )

    parser.add_argument(
        "--continue-on-error",
        action="store_true",
        help=(
            "Continue to the next model if one model fails"
        ),
    )

    return parser.parse_args()


def result_path(model_id):
    return (
        RESULT_ROOT
        / model_id
        / "result.json"
    )


def load_result(model_id):
    path = result_path(model_id)

    if not path.exists():
        return None

    return json.loads(
        path.read_text()
    )


def write_summary(manifest):
    rows = []

    for model in manifest.itertuples(
        index=False
    ):
        result = load_result(
            model.ModelID
        )

        row = {
            "ModelID": model.ModelID,
            "CellLineName": model.CellLineName,
            "Lineage": model.Lineage,
            "PublishedFPBM": model.PublishedFPBM,
            "WRNGeneEffect": model.WRNGeneEffect,
            "Status": (
                "complete"
                if result is not None
                else "pending"
            ),
            "CurrentBrokenTAFPBM": None,
            "CurrentNonbrokenTAFPBM": None,
            "CurrentToPublishedRatio": None,
            "WallMinutes": None,
            "FetchRetries": None,
            "ResumedJobs": None,
            "CoverageSHA256": None,
            "ResultSchema": None,
            "ExtractorGitCommit": None,
        }

        if result is not None:
            row.update(
                {
                    "CurrentBrokenTAFPBM":
                        result.get(
                            "current_broken_ta_mean_fpbm"
                        ),
                    "CurrentNonbrokenTAFPBM":
                        result.get(
                            "current_nonbroken_ta_mean_fpbm"
                        ),
                    "CurrentToPublishedRatio":
                        result.get(
                            "current_to_published_ratio"
                        ),
                    "WallMinutes":
                        (
                            result.get(
                                "wall_seconds"
                            ) / 60
                            if result.get(
                                "wall_seconds"
                            ) is not None
                            else None
                        ),
                    "FetchRetries":
                        result.get(
                            "fetch_retry_count"
                        ),
                    "ResumedJobs":
                        result.get(
                            "resumed_job_count"
                        ),
                    "CoverageSHA256":
                        result.get(
                            "coverage_sha256"
                        ),
                    "ResultSchema":
                        result.get(
                            "schema_version"
                        ),
                    "ExtractorGitCommit":
                        result.get(
                            "extractor_git_commit"
                        ),
                }
            )

        rows.append(row)

    summary = pd.DataFrame(rows)

    RESULT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    summary.to_csv(
        SUMMARY,
        index=False,
    )

    complete = (
        summary["Status"]
        .eq("complete")
        .sum()
    )

    print(
        f"Cohort summary: "
        f"{complete}/{len(summary)} complete"
    )

    print(
        f"Summary: {SUMMARY}"
    )


def main():
    args = parse_args()

    if args.workers < 1:
        raise RuntimeError(
            "--workers must be at least 1"
        )

    if (
        args.limit is not None
        and args.limit < 1
    ):
        raise RuntimeError(
            "--limit must be at least 1"
        )

    manifest = pd.read_csv(
        MANIFEST
    )

    write_summary(manifest)

    pending = [
        row
        for row in manifest.itertuples(
            index=False
        )
        if not result_path(
            row.ModelID
        ).exists()
    ]

    if args.limit is not None:
        pending = pending[
            :args.limit
        ]

    print()
    print(
        "CURRENT-WGS COHORT RUNNER"
    )
    print("=" * 72)
    print(
        f"Manifest models: {len(manifest)}"
    )
    print(
        f"Already complete: "
        f"{len(manifest) - len([
            row for row in manifest.itertuples(index=False)
            if not result_path(row.ModelID).exists()
        ])}"
    )
    print(
        f"Models this run: {len(pending)}"
    )
    print(
        f"Workers/model: {args.workers}"
    )
    print()

    if not pending:
        print(
            "Nothing to run."
        )
        return

    failures = []

    for number, model in enumerate(
        pending,
        start=1,
    ):
        print()
        print("=" * 72)
        print(
            f"MODEL {number}/{len(pending)}"
        )
        print(
            f"{model.ModelID} "
            f"{model.CellLineName}"
        )
        print("=" * 72)

        cmd = [
            sys.executable,
            "-u",
            str(EXTRACTOR),
            model.ModelID,
            "--workers",
            str(args.workers),
        ]

        try:
            subprocess.run(
                cmd,
                check=True,
            )

        except subprocess.CalledProcessError as exc:
            failures.append(
                model.ModelID
            )

            print(
                f"FAILED: {model.ModelID} "
                f"(exit {exc.returncode})",
                flush=True,
            )

            if not args.continue_on_error:
                write_summary(
                    manifest
                )
                raise

        write_summary(
            manifest
        )

    print()
    print("=" * 72)
    print("COHORT RUN COMPLETE")
    print("=" * 72)

    if failures:
        print(
            "Failed models:",
            ", ".join(failures),
        )
    else:
        print(
            "All requested models completed."
        )


if __name__ == "__main__":
    main()
