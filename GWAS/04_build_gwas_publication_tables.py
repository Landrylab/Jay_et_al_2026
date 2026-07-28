#!/usr/bin/env python3
"""Export compact publication-facing GWAS tables.

This wrapper runs the submission exporter without rebuilding the four large
compressed point-level Source Data files.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-root", required=True, type=Path)
    parser.add_argument("--figure-root", required=True, type=Path)
    parser.add_argument("--submission-root", required=True, type=Path)
    args = parser.parse_args()

    exporter = Path(__file__).resolve().parent / "07_export_submission_data.py"
    subprocess.check_call(
        [
            sys.executable,
            str(exporter),
            "--results-root",
            str(args.results_root),
            "--figure-root",
            str(args.figure_root),
            "--submission-root",
            str(args.submission_root),
            "--skip-source-data",
        ]
    )


if __name__ == "__main__":
    main()
