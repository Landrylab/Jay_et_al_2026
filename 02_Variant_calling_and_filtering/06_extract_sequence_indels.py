#!/usr/bin/env python3
r"""Represent sequence indel ALTs as diploid target-allele dosages in PLINK BED.

Multiallelic records are split. Each allele retains its original coordinate and parent-record
identity. Other called ALTs represent zero copies of the target allele;
partially missing and non-diploid calls are missing.

Example:
  python 06_extract_sequence_indels.py --vcf path/to/cohort.ann.vcf.gz \
    --out-prefix path/to/output/sequence_indels

All VCF samples are retained by default. Use --samples to select and order
sample IDs, or --sample-fam to use an existing six-column FAM unchanged.

Optional --exclude-bim omits exact sequence alleles represented in another
marker panel. Its variant identifiers must have CHROM:POS:REF:ALT syntax,
where ALT may contain the complete original comma-separated ALT list.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile


SEQUENCE = re.compile(r"[ACGTNacgtn]+")
GT_CODE = {
    "0/0": 3, "0|0": 3,
    "0/1": 2, "1/0": 2, "0|1": 2, "1|0": 2,
    "1/1": 0, "1|1": 0,
}
SUFFIXES = (".bed", ".bim", ".fam", ".alleles.tsv.gz", ".audit.json")
MANIFEST_FIELDS = (
    "chromosome", "marker_id", "position", "ref", "target_alt",
    "original_parent_id", "original_alt_index", "original_chromosome",
    "original_position", "original_ref", "original_alt_list", "record_split",
    "marker_index", "n_homalt", "n_heterozygous", "n_target_absent", "n_missing",
)


def read_fam(path: Path) -> list[str]:
    samples = []
    for number, line in enumerate(path.read_text().splitlines(), 1):
        fields = line.split()
        if len(fields) != 6:
            raise ValueError(f"FAM row {number} must contain six fields")
        samples.append(fields[1])
    if not samples or len(set(samples)) != len(samples):
        raise ValueError("FAM must contain at least one sample and unique individual IDs")
    return samples


def read_samples(path: Path) -> list[str]:
    samples = []
    for number, line in enumerate(path.read_text().splitlines(), 1):
        if not line.strip():
            continue
        fields = line.split()
        if len(fields) != 1:
            raise ValueError(f"Sample-list row {number} must contain exactly one sample ID")
        samples.append(fields[0])
    if not samples or len(set(samples)) != len(samples):
        raise ValueError("Sample list must contain at least one sample and unique IDs")
    return samples


def read_chromosome_map(path: Path | None) -> dict[str, str] | None:
    if path is None:
        return None
    mapping = {}
    for number, line in enumerate(path.read_text().splitlines(), 1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.split()
        if len(fields) != 2 or fields[0] in mapping:
            raise ValueError(f"Invalid or repeated chromosome mapping at row {number}")
        mapping[fields[0]] = fields[1]
    if not mapping or len(set(mapping.values())) != len(mapping):
        raise ValueError("Chromosome mapping must be nonempty and one-to-one")
    return mapping


def represented_keys(path: Path | None) -> set[tuple[str, int, str, str]]:
    if path is None:
        return set()
    keys = set()
    with path.open() as handle:
        for number, line in enumerate(handle, 1):
            fields = line.split()
            if len(fields) != 6:
                raise ValueError(f"BIM row {number} must contain six fields")
            _chrom, marker, _cm, pos, a1, a2 = fields
            if not SEQUENCE.fullmatch(a1) or not SEQUENCE.fullmatch(a2) or len(a1) == len(a2):
                continue
            identity = marker.split(":")
            if len(identity) != 4 or identity[1] != pos:
                raise ValueError(f"BIM indel at row {number} needs a CHROM:POS:REF:ALT identifier")
            chrom, _position, ref, alts = identity
            if a2 == ref:
                alt = a1
            elif a1 == ref:
                alt = a2
            else:
                raise ValueError(f"BIM alleles do not contain the identifier REF at row {number}")
            if alt not in alts.split(","):
                raise ValueError(f"BIM target allele is absent from the identifier ALT list at row {number}")
            keys.add((chrom, int(pos), ref, alt))
    return keys


def pack_calls(calls: list[str]) -> tuple[bytes, Counter]:
    packed = bytearray((len(calls) + 3) // 4)
    counts = Counter()
    for index, call in enumerate(calls):
        code = GT_CODE.get(call, 1)
        counts[code] += 1
        packed[index // 4] |= code << (2 * (index % 4))
        if code == 1 and "." not in call:
            counts["non_diploid_or_unexpected"] += 1
    return bytes(packed), counts


def parent_identity(chrom: str, pos: int, ref: str, alt: str, original: str) -> tuple:
    if original == ".":
        return chrom, pos, ref, alt, 1, False
    fields = original.split("|")
    if len(fields) != 5:
        raise ValueError(f"Unexpected original-record tag: {original}")
    parent_chrom, parent_pos, parent_ref, parent_alts, alt_index = fields
    parent_pos, alt_index = int(parent_pos), int(alt_index)
    alleles = parent_alts.split(",")
    if (
        (parent_chrom, parent_pos, parent_ref) != (chrom, pos, ref)
        or not 1 <= alt_index <= len(alleles)
        or alleles[alt_index - 1] != alt
    ):
        raise ValueError("Split allele differs from its original coordinate, REF, or ALT")
    return parent_chrom, parent_pos, parent_ref, parent_alts, alt_index, True


def extract(args: argparse.Namespace) -> dict:
    executable = shutil.which(args.bcftools)
    if executable is None:
        raise ValueError(f"bcftools executable not found: {args.bcftools}")
    for path in (args.vcf, args.sample_fam, args.samples, args.exclude_bim, args.chromosome_map):
        if path is not None and not path.is_file():
            raise ValueError(f"Input file not found: {path}")
    mapping = read_chromosome_map(args.chromosome_map)
    exclusions = represented_keys(args.exclude_bim)
    output_paths = [Path(str(args.out_prefix) + suffix) for suffix in SUFFIXES]
    existing = [str(path) for path in output_paths if path.exists()]
    if existing:
        raise ValueError("Output files already exist: " + ", ".join(existing))
    sample_query = subprocess.run(
        [executable, "query", "-l", str(args.vcf)], check=True, capture_output=True, text=True
    )
    vcf_samples = sample_query.stdout.splitlines()
    if not vcf_samples or len(set(vcf_samples)) != len(vcf_samples):
        raise ValueError("VCF must contain at least one sample and unique sample IDs")
    if args.sample_fam:
        samples = read_fam(args.sample_fam)
        sample_order = "Individual IDs and order from the supplied FAM; FAM copied unchanged."
    elif args.samples:
        samples = read_samples(args.samples)
        sample_order = "IDs and order from the supplied sample list; generated FAM uses FID=IID, parents=0, sex=0, phenotype=-9."
    else:
        samples = vcf_samples
        sample_order = "All samples in VCF order; generated FAM uses FID=IID, parents=0, sex=0, phenotype=-9."
    if any(any(character.isspace() for character in sample) for sample in samples):
        raise ValueError("Sample IDs must not contain whitespace")
    lookup = {sample: index for index, sample in enumerate(vcf_samples)}
    absent = [sample for sample in samples if sample not in lookup]
    if absent:
        raise ValueError("Selected sample IDs absent from VCF: " + ", ".join(absent[:10]))
    indices = [lookup[sample] for sample in samples]
    tag = "INDEL_PARENT_RECORD"
    header = subprocess.run(
        [executable, "view", "-h", str(args.vcf)], check=True, capture_output=True, text=True
    ).stdout
    if f"##INFO=<ID={tag}," in header:
        raise ValueError(f"VCF already defines INFO/{tag}; use an undecomposed source VCF")
    commands = [
        [executable, "view", "-v", "indels", "-Ou", str(args.vcf)],
        [executable, "norm", "-m", "-any", "-N", "--multi-overlaps", "0",
         "--old-rec-tag", tag, "-Ou"],
        [executable, "query", "-f", f"%CHROM\t%POS\t%REF\t%ALT\t%INFO/{tag}[\t%GT]\n"],
    ]
    counts = Counter({"source_samples": len(vcf_samples), "selected_samples": len(samples),
                      "excluded_exact_allele_keys": len(exclusions), "alleles_written": 0})
    distinct_positions = set()
    seen: dict[tuple[str, int, str, str], str] = {}
    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".sequence-indels-", dir=args.out_prefix.parent) as temporary:
        stage = Path(temporary)
        staged = [stage / ("output" + suffix) for suffix in SUFFIXES]
        logs = [open(stage / f"command_{index}.stderr", "w+") for index in range(3)]
        processes = []
        try:
            upstream = None
            for index, command in enumerate(commands):
                process = subprocess.Popen(command, stdin=upstream, stdout=subprocess.PIPE,
                                           stderr=logs[index], text=(index == 2))
                if upstream is not None:
                    upstream.close()
                processes.append(process)
                upstream = process.stdout
            with staged[0].open("wb") as bed, staged[1].open("w") as bim, gzip.open(staged[3], "wt") as manifest:
                bed.write(bytes((108, 27, 1)))
                writer = csv.DictWriter(manifest, fieldnames=MANIFEST_FIELDS, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                for line in processes[-1].stdout:
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) != 5 + len(vcf_samples):
                        raise ValueError("VCF query returned an unexpected number of sample fields")
                    chrom, position, ref, alt, original = fields[:5]
                    pos = int(position)
                    counts["split_rows_seen"] += 1
                    if not SEQUENCE.fullmatch(ref) or not SEQUENCE.fullmatch(alt):
                        counts["symbolic_or_nonsequence_rows_skipped"] += 1
                        continue
                    if len(ref) == len(alt):
                        counts["equal_length_rows_skipped"] += 1
                        continue
                    counts["sequence_indel_rows"] += 1
                    key = (chrom, pos, ref, alt)
                    if key in exclusions:
                        counts["represented_exact_alleles_skipped"] += 1
                        continue
                    calls = [fields[5 + index] for index in indices]
                    packed, genotype_counts = pack_calls(calls)
                    digest = hashlib.sha256(packed).hexdigest()
                    if key in seen:
                        if seen[key] != digest:
                            raise ValueError(f"Duplicate sequence allele has discordant genotypes: {key}")
                        counts["duplicate_exact_alleles_skipped"] += 1
                        continue
                    seen[key] = digest
                    parent_chrom, parent_pos, parent_ref, parent_alts, alt_index, split = parent_identity(
                        chrom, pos, ref, alt, original
                    )
                    if mapping is not None and chrom not in mapping:
                        raise ValueError(f"Chromosome is absent from --chromosome-map: {chrom}")
                    output_chrom = chrom if mapping is None else mapping[chrom]
                    marker = f"{chrom}:{pos}:{ref}:{alt}"
                    parent = f"{parent_chrom}:{parent_pos}:{parent_ref}:{parent_alts}"
                    marker_index = counts["alleles_written"]
                    bed.write(packed)
                    bim.write(f"{output_chrom}\t{marker}\t0\t{pos}\t{alt}\t{ref}\n")
                    writer.writerow(dict(zip(MANIFEST_FIELDS, (
                        output_chrom, marker, pos, ref, alt, parent, alt_index, parent_chrom,
                        parent_pos, parent_ref, parent_alts, split, marker_index,
                        genotype_counts[0], genotype_counts[2], genotype_counts[3], genotype_counts[1],
                    ))))
                    counts["alleles_written"] += 1
                    counts["insertions_by_length" if len(alt) > len(ref) else "deletions_by_length"] += 1
                    counts["alleles_with_multiallelic_parent"] += "," in parent_alts
                    counts["non_diploid_or_unexpected_calls_made_missing"] += genotype_counts["non_diploid_or_unexpected"]
                    distinct_positions.add((chrom, pos))
            statuses = [process.wait() for process in processes]
            errors = []
            for command, status, log in zip(commands, statuses, logs):
                if status:
                    log.seek(0)
                    errors.append(f"Command failed ({status}): {command}\n{log.read()[-4000:]}")
            if errors:
                raise ValueError("\n".join(errors))
            logs[1].seek(0)
            normalization_log = logs[1].read()
            normalization_counts = re.search(r"Lines\s+total/[^:\n]*:\s*(\d+)/(\d+)", normalization_log)
            if normalization_counts:
                counts["source_indel_containing_records"] = int(normalization_counts[1])
                counts["source_multiallelic_records_split"] = int(normalization_counts[2])
        finally:
            for process in processes:
                if process.poll() is None:
                    process.terminate()
            for process in processes:
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
            for log in logs:
                log.close()
        if args.sample_fam:
            shutil.copyfile(args.sample_fam, staged[2])
        else:
            staged[2].write_text("".join(f"{sample}\t{sample}\t0\t0\t0\t-9\n" for sample in samples))
        counts["distinct_positions_written"] = len(distinct_positions)
        expected_bytes = 3 + counts["alleles_written"] * ((len(samples) + 3) // 4)
        if staged[0].stat().st_size != expected_bytes:
            raise ValueError("BED size does not match the sample and allele counts")
        report = {
            "counts": dict(counts), "source_vcf": str(args.vcf),
            "sample_fam": str(args.sample_fam) if args.sample_fam else None,
            "samples_file": str(args.samples) if args.samples else None,
            "exclude_bim": str(args.exclude_bim) if args.exclude_bim else None,
            "chromosome_map": mapping, "commands": commands,
            "normalization": "Multiallelic splitting only.",
            "genotype_coding": "Target ALT dosage 0/1/2; other called ALTs non-target; partial/non-diploid calls missing.",
            "sample_order": sample_order,
            "allele_definition": "Sequence REF/ALT containing A,C,G,T,N; unequal lengths; exact representation keys.",
            "position_definition": "Number of distinct original CHROM/POS pairs among written alleles.",
            "outputs": [str(path) for path in output_paths],
            "normalization_log": normalization_log.strip(),
        }
        staged[4].write_text(json.dumps(report, indent=2) + "\n")
        published = []
        try:
            for source, target in zip(staged, output_paths):
                # Exclusive links prevent an existing output from being overwritten.
                os.link(source, target)
                published.append(target)
        except BaseException:
            for target in published:
                target.unlink()
            raise
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--vcf", type=Path, required=True, help="Joint VCF containing sequence indels and GT fields")
    sample_selection = parser.add_mutually_exclusive_group()
    sample_selection.add_argument("--sample-fam", type=Path, help="Six-column FAM defining individual IDs and order")
    sample_selection.add_argument("--samples", type=Path, help="One sample ID per line, in desired order; default retains all VCF samples")
    parser.add_argument("--out-prefix", type=Path, required=True, help="Prefix for BED, BIM, FAM, allele table, and audit")
    parser.add_argument("--exclude-bim", type=Path, help="Marker panel whose exact sequence indels should be excluded")
    parser.add_argument("--chromosome-map", type=Path, help="Headerless source-to-output chromosome mapping, two columns")
    parser.add_argument("--bcftools", default="bcftools", help="bcftools executable with --multi-overlaps support")
    args = parser.parse_args()
    try:
        report = extract(args)
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        parser.exit(1, f"Error: {error}\n")
    print(json.dumps(report["counts"], indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
