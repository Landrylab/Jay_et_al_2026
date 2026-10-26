#!/usr/bin/env python3
"""Export final ploidy, chromosome-CN, gene-CN, and 10 kb source data."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import shutil
from pathlib import Path

import numpy as np
import pandas as pd


BIN_SUFFIXES = (".tsv", ".binCN.tsv", ".bins.cn.tsv", ".bins.tsv", ".cn.tsv")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bin-cn-dir", required=True)
    parser.add_argument("--ploidy", required=True)
    parser.add_argument("--chromosome-cn", required=True)
    parser.add_argument("--gene-cnv", required=True)
    parser.add_argument("--loci", required=True)
    parser.add_argument("--manual-review", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--source-data", required=True)
    parser.add_argument(
        "--gene-copy-number-out",
        default=None,
        help="Optional gzip-compressed copy of the final gene-CN matrix",
    )
    parser.add_argument("--expected-samples", type=int, default=None)
    return parser.parse_args()


def resolve_column(frame: pd.DataFrame, names: tuple[str, ...]) -> str:
    lowered = {column.lower(): column for column in frame.columns}
    for name in names:
        if name.lower() in lowered:
            return lowered[name.lower()]
    raise ValueError(f"Missing required column; accepted names: {names}")


def find_bin_file(directory: Path, sample: str) -> Path:
    for suffix in BIN_SUFFIXES:
        candidate = directory / f"{sample}{suffix}"
        if candidate.exists() and candidate.stat().st_size:
            return candidate
    raise FileNotFoundError(f"No 10 kb copy-number file found for {sample}")


def sha256sum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_final_calls(args: argparse.Namespace):
    ploidy = pd.read_csv(args.ploidy, sep="\t")
    sample_column = resolve_column(ploidy, ("sample_id", "sample", "isolate"))
    ploidy_column = resolve_column(
        ploidy, ("ploidy_hat", "baseline_ploidy", "ploidy")
    )
    ploidy = ploidy.set_index(sample_column)
    baseline = pd.to_numeric(ploidy[ploidy_column], errors="raise")
    if baseline.index.has_duplicates or (args.expected_samples is not None and len(baseline) != args.expected_samples):
        raise ValueError("Final ploidy table has duplicate samples or does not match --expected-samples.")
    if baseline.empty or not baseline.isin([2, 3, 4, 6]).all():
        raise ValueError("Final baseline assignments must be 2N, 3N, 4N or 6N.")

    chromosome_cn = pd.read_csv(args.chromosome_cn, sep="\t")
    chromosome_sample = resolve_column(
        chromosome_cn, ("sample_id", "sample", "isolate")
    )
    chromosome_cn = chromosome_cn.set_index(chromosome_sample)
    if set(chromosome_cn.index) != set(baseline.index):
        raise ValueError("Ploidy and chromosome-CN sample sets differ.")

    manual_review = pd.read_csv(args.manual_review, sep="\t")
    manual_sample_column = resolve_column(
        manual_review, ("sample_id", "sample", "isolate")
    )
    manual_review = manual_review.set_index(manual_sample_column)
    if manual_review.index.has_duplicates:
        raise ValueError("Manual-review table has duplicate samples.")
    if not set(manual_review.index).issubset(baseline.index):
        raise ValueError("Manual-review table contains unknown isolates.")
    if not manual_review.empty:
        review_ploidy_column = resolve_column(manual_review, ("baseline_ploidy", "ploidy_hat", "ploidy"))
        reviewed = pd.to_numeric(manual_review[review_ploidy_column], errors="raise")
        if not reviewed.eq(baseline.loc[reviewed.index]).all():
            raise ValueError("Final baseline differs from the accepted review; rerun sample calling with the reviewed baseline.")
    return ploidy, baseline, chromosome_cn, manual_review


def clean_ploidy_table(
    baseline: pd.Series,
    manual_review: pd.DataFrame,
) -> pd.DataFrame:
    output = baseline.rename("baseline_ploidy").to_frame().reset_index()
    output = output.rename(columns={output.columns[0]: "sample_id"})
    manual_samples = set(manual_review.index)
    output["assignment_method"] = output["sample_id"].map(
        lambda sample: (
            "manual_review"
            if sample in manual_samples
            else "candidate_allele_fraction_baseline"
        )
    )
    return output


def clean_chromosome_table(
    chromosome_cn: pd.DataFrame,
    baseline: pd.Series,
) -> pd.DataFrame:
    chromosome_columns = list(chromosome_cn.columns)
    if not chromosome_columns:
        raise ValueError("The chromosome-CN table has no chromosome columns.")
    output = chromosome_cn.loc[baseline.index, chromosome_columns].copy()
    output.insert(0, "baseline_ploidy", baseline)
    output.index.name = "sample_id"
    return output.reset_index()


def clean_locus_table(loci_path: str) -> pd.DataFrame:
    loci = pd.read_csv(loci_path, sep=None, engine="python")
    required = [
        "locus_id",
        "gene_ids",
        "chromosome",
        "chromosome_label",
        "start_1based_inclusive",
        "end_1based_inclusive",
    ]
    missing = set(required).difference(loci.columns)
    if missing:
        raise ValueError(f"Locus table is missing: {sorted(missing)}")
    output = loci.loc[:, required].dropna(
        subset=[
            "locus_id",
            "gene_ids",
            "chromosome",
            "start_1based_inclusive",
            "end_1based_inclusive",
        ]
    ).copy()
    output["start_1based_inclusive"] = pd.to_numeric(
        output["start_1based_inclusive"], errors="raise"
    ).astype(int)
    output["end_1based_inclusive"] = pd.to_numeric(
        output["end_1based_inclusive"], errors="raise"
    ).astype(int)
    if output["locus_id"].duplicated().any():
        raise ValueError("Locus identifiers are not unique.")
    return output.sort_values(
        ["chromosome_label", "start_1based_inclusive", "end_1based_inclusive"]
    ).reset_index(drop=True)


def normalize_bin_table(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    lowered = {column.lower(): column for column in frame.columns}
    chromosome_column = next(
        (
            lowered[name]
            for name in ("contig", "chr", "chrom", "chromosome")
            if name in lowered
        ),
        None,
    )
    start_column = lowered.get("start")
    end_column = lowered.get("end")
    copy_number_column = next(
        (
            lowered[name]
            for name in ("cn", "copy_number", "cn_bin", "cn_abs_round", "cn_abs")
            if name in lowered
        ),
        None,
    )
    if None in {
        chromosome_column,
        start_column,
        end_column,
        copy_number_column,
    }:
        raise ValueError(f"Unrecognized 10 kb copy-number columns: {path}")
    output = pd.DataFrame(
        {
            "chromosome": frame[chromosome_column].astype(str),
            "window_start_1based_inclusive": pd.to_numeric(
                frame[start_column], errors="raise"
            ).astype(int),
            "window_end_1based_inclusive": pd.to_numeric(
                frame[end_column], errors="raise"
            ).astype(int),
            "copy_number": pd.to_numeric(
                frame[copy_number_column], errors="raise"
            ).astype(float),
        }
    )
    for source, target in (
        ("n_bins", "number_1kb_bins"),
        ("n_above", "number_1kb_bins_above_baseline"),
        ("n_below", "number_1kb_bins_below_baseline"),
        ("thresh", "directional_support_threshold"),
    ):
        if source in lowered:
            output[target] = pd.to_numeric(frame[lowered[source]], errors="coerce")
        else:
            output[target] = np.nan
    if "baseline_ploidy" in lowered:
        output["track_baseline_ploidy"] = pd.to_numeric(frame[lowered["baseline_ploidy"]], errors="raise")
    if output["copy_number"].isna().any() or not np.isfinite(output["copy_number"]).all():
        raise ValueError(f"Uncalled copy-number rows in {path}; omit unsupported windows.")
    return output


def export_10kb_source_data(
    bin_directory: Path,
    source_path: Path,
    baseline: pd.Series,
    chromosome_cn: pd.DataFrame,
    manual_review: pd.DataFrame,
) -> tuple[int, int]:
    source_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "sample_id",
        "chromosome",
        "window_start_1based_inclusive",
        "window_end_1based_inclusive",
        "copy_number",
        "baseline_ploidy",
        "number_1kb_bins",
        "number_1kb_bins_above_baseline",
        "number_1kb_bins_below_baseline",
        "directional_support_threshold",
        "manual_review_applied",
    ]
    manual_samples = set(manual_review.index)
    row_count = 0
    with gzip.open(source_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for sample in baseline.index:
            frame = normalize_bin_table(find_bin_file(bin_directory, sample))
            sample_baseline = float(baseline.loc[sample])
            if sample in manual_samples and "track_baseline_ploidy" not in frame:
                raise ValueError(f"{sample}: reviewed sample requires a track baseline; aggregate the reviewed segment calls before exporting.")
            if "track_baseline_ploidy" in frame and not np.all(frame["track_baseline_ploidy"] == sample_baseline):
                raise ValueError(f"{sample}: 10-kb track and final baseline differ; aggregate reviewed segment calls.")
            for row in frame.itertuples(index=False):
                writer.writerow(
                    {
                        "sample_id": sample,
                        "chromosome": row.chromosome,
                        "window_start_1based_inclusive":
                            row.window_start_1based_inclusive,
                        "window_end_1based_inclusive":
                            row.window_end_1based_inclusive,
                        "copy_number": row.copy_number,
                        "baseline_ploidy": sample_baseline,
                        "number_1kb_bins": row.number_1kb_bins,
                        "number_1kb_bins_above_baseline":
                            row.number_1kb_bins_above_baseline,
                        "number_1kb_bins_below_baseline":
                            row.number_1kb_bins_below_baseline,
                        "directional_support_threshold":
                            row.directional_support_threshold,
                        "manual_review_applied": sample in manual_samples,
                    }
                )
                row_count += 1
    return row_count, len(baseline)


def gzip_copy(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with source.open("rb") as input_handle, gzip.open(
        destination, "wb", compresslevel=9
    ) as output_handle:
        shutil.copyfileobj(input_handle, output_handle)


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    source_path = Path(args.source_data)

    _, baseline, chromosome_cn, manual_review = load_final_calls(args)
    ploidy_output = clean_ploidy_table(baseline, manual_review)
    chromosome_output = clean_chromosome_table(chromosome_cn, baseline)
    locus_output = clean_locus_table(args.loci)

    ploidy_path = output_dir / "cnv_ploidy_calls.tsv"
    chromosome_path = output_dir / "cnv_chromosome_copy_number.tsv"
    locus_path = output_dir / "cnv_gene_loci.tsv"
    manual_path = output_dir / "cnv_manual_review_corrections.tsv"
    ploidy_output.to_csv(ploidy_path, sep="\t", index=False)
    chromosome_output.to_csv(chromosome_path, sep="\t", index=False)
    locus_output.to_csv(locus_path, sep="\t", index=False)
    manual_review.reset_index().rename(
        columns={manual_review.index.name: "sample_id"}
    ).to_csv(manual_path, sep="\t", index=False)

    source_rows, source_samples = export_10kb_source_data(
        Path(args.bin_cn_dir),
        source_path,
        baseline,
        chromosome_cn,
        manual_review,
    )

    gene_output = (
        Path(args.gene_copy_number_out)
        if args.gene_copy_number_out
        else source_path.parent / "CNV_gene_copy_number_matrix.tsv.gz"
    )
    gzip_copy(Path(args.gene_cnv), gene_output)

    manifest = pd.DataFrame(
        [
            {
                "file": source_path.name,
                "content": "Main Figure 1G 10 kb copy-number states",
                "rows": source_rows,
                "samples": source_samples,
                "sha256": sha256sum(source_path),
            },
            {
                "file": gene_output.name,
                "content": "Final manually reviewed gene-by-isolate copy-number matrix",
                "rows": sum(1 for _ in gzip.open(gene_output, "rt")) - 1,
                "samples": len(baseline),
                "sha256": sha256sum(gene_output),
            },
        ]
    )
    manifest.to_csv(
        source_path.parent / "cnv_source_data_manifest.tsv",
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    main()
