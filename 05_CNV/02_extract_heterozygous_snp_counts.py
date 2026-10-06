#!/usr/bin/env python3
"""Extract per-sample AD counts at undecomposed biallelic heterozygous SNPs."""

import argparse
import csv
import json
import re
import subprocess
from collections import Counter
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True, help="Undecomposed VCF/BCF with GT and FORMAT/AD")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--minimum-allele-depth", type=int, default=1)
    parser.add_argument("--minimum-total-depth", type=int, default=1)
    parser.add_argument("--exclude-contig-regex", default="(?i)chrM|mitochond")
    parser.add_argument("--bcftools", default="bcftools")
    args = parser.parse_args()
    if min(args.minimum_allele_depth, args.minimum_total_depth) < 0:
        parser.error("Depth thresholds must be non-negative.")
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    excluded, counts, seen = re.compile(args.exclude_contig_regex), Counter(), set()
    # Allele cardinality is inspected before any decomposition or sample trimming.
    command = [args.bcftools, "query", "-s", args.sample, "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%AD]\n", args.vcf]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, text=True)
    temporary = out.with_name(out.name + ".partial")
    try:
        with temporary.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["chr", "pos", "ref_count", "alt_count", "sample_id", "ref", "alt", "GT"])
            for line in process.stdout:
                counts["records_examined"] += 1
                chrom, pos, ref, alt, gt, ad = line.rstrip("\n").split("\t")
                if excluded.search(chrom):
                    counts["non_nuclear"] += 1
                    continue
                if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT" or ref == alt:
                    counts["not_biallelic_snp"] += 1
                    continue
                if re.split(r"[/|]", gt) not in (["0", "1"], ["1", "0"]):
                    counts["not_diploid_heterozygote"] += 1
                    continue
                fields = ad.split(",")
                if len(fields) != 2 or not all(value.isdigit() for value in fields):
                    counts["missing_or_invalid_AD"] += 1
                    continue
                ref_count, alt_count = map(int, fields)
                if min(ref_count, alt_count) < args.minimum_allele_depth or ref_count + alt_count < args.minimum_total_depth or ref_count + alt_count == 0:
                    counts["insufficient_allele_depth"] += 1
                    continue
                key = (chrom, int(pos))
                if key in seen:
                    raise ValueError(f"Repeated qualifying position {chrom}:{pos}; use one undecomposed record per position.")
                seen.add(key)
                writer.writerow([chrom, pos, ref_count, alt_count, args.sample, ref, alt, gt])
                counts["retained"] += 1
        if process.wait() != 0:
            raise RuntimeError("bcftools query failed; verify the sample ID and FORMAT/AD field.")
        temporary.replace(out)
    except BaseException:
        process.kill()
        process.wait()
        temporary.unlink(missing_ok=True)
        raise
    with out.with_suffix(out.suffix + ".qc.json").open("w") as handle:
        json.dump({"sample_id": args.sample, "minimum_allele_depth": args.minimum_allele_depth, "minimum_total_depth": args.minimum_total_depth, "counts": dict(counts)}, handle, indent=2)
        handle.write("\n")


if __name__ == "__main__":
    main()
