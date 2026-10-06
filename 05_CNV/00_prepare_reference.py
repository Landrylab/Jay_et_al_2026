#!/usr/bin/env python3
"""Prepare nuclear gene intervals, merged exclusion masks and locus metadata."""

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote


def write_table(path, columns, rows, header=True):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        if header:
            writer.writerow(columns)
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fai", required=True)
    parser.add_argument("--gff", required=True)
    parser.add_argument("--gene-types", default="gene,ORF,pseudogene")
    parser.add_argument("--mask-types", default="centromere,telomere,repeat_region,rRNA_gene")
    parser.add_argument("--mask-bed", action="append", default=[], help="Additional BED mask; repeat as needed")
    parser.add_argument("--exclude-contig-regex", default="(?i)chrM|mitochond")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    excluded = re.compile(args.exclude_contig_regex)
    lengths = {}
    with open(args.fai) as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if not excluded.search(fields[0]):
                lengths[fields[0]] = int(fields[1])
    genes, masks = {}, defaultdict(list)
    gene_types, mask_types = set(args.gene_types.split(",")), set(args.mask_types.split(","))
    with open(args.gff) as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip().split("\t")
            chrom, feature = fields[0], fields[2]
            if chrom not in lengths:
                continue
            start, end = int(fields[3]) - 1, int(fields[4])
            if not 0 <= start < end <= lengths[chrom]:
                raise ValueError(f"GFF interval outside reference: {chrom}:{start}-{end}")
            if feature in mask_types:
                masks[chrom].append((start, end))
            if feature in gene_types:
                attributes = dict(item.split("=", 1) for item in fields[8].split(";") if "=" in item)
                identifier = unquote(attributes.get("ID", ""))
                if not identifier or identifier in genes:
                    raise ValueError(f"Missing or duplicate gene ID: {identifier}")
                genes[identifier] = (chrom, start, end, identifier)
    for path in args.mask_bed:
        with open(path) as handle:
            for line in handle:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                chrom, first, last = line.split()[:3]
                if chrom not in lengths:
                    continue
                start, end = int(first), int(last)
                if not 0 <= start < end <= lengths[chrom]:
                    raise ValueError(f"BED interval outside reference: {chrom}:{start}-{end}")
                masks[chrom].append((start, end))
    if not genes:
        raise ValueError("No genes match --gene-types on the selected reference contigs.")
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    order = {chrom: i for i, chrom in enumerate(lengths)}
    gene_rows = sorted(genes.values(), key=lambda row: (order[row[0]], row[1], row[2], row[3]))
    write_table(out / "genes.tsv", ["chr", "start", "end", "gene_id"], gene_rows)
    grouped = defaultdict(list)
    for chrom, start, end, identifier in gene_rows:
        grouped[(chrom, start, end)].append(identifier)
    locus_rows = []
    for (chrom, start, end), identifiers in grouped.items():
        match = re.search(r"chr([1-7R])A", chrom)
        label = f"Chr{match.group(1)}" if match else chrom
        locus_rows.append((";".join(identifiers), ";".join(identifiers), chrom, label, start + 1, end))
    write_table(out / "gene_loci.tsv", ["locus_id", "gene_ids", "chromosome", "chromosome_label", "start_1based_inclusive", "end_1based_inclusive"], locus_rows)
    merged = []
    for chrom in lengths:
        for start, end in sorted(masks[chrom]):
            if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
                merged[-1][2] = max(merged[-1][2], end)
            else:
                merged.append([chrom, start, end])
    write_table(out / "mask.bed", [], merged, header=False)
    write_table(out / "nuclear.genome.tsv", [], lengths.items(), header=False)


if __name__ == "__main__":
    main()
