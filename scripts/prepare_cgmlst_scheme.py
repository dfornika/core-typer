#!/usr/bin/env python3
"""
Combine a directory of per-locus allele FASTA files (as downloaded from
cgmlst.org, one file per locus named "<locus_id>.fasta", allele headers just
the allele number) into a single multi-FASTA with headers rewritten as
"<locus_id>_<allele_id>_<md5_hash>", ready for `kma_index`.

The md5 hash is a content-derived identifier for the allele's own sequence
(the same convention CGE cgMLSTFinder uses for its "hypothetical new allele"
output). core-typer's parsers recognize and use this hash suffix when
present, but tolerate plain "<locus_id>_<allele_id>" schemes too.
"""

import argparse
import glob
import os

from core_typer import allele_calling


def iter_locus_fasta_files(input_dir):
    return sorted(glob.glob(os.path.join(input_dir, "*.fasta")))


def iter_fasta_records(fasta_path):
    header = None
    seq_lines = []
    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_lines)
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        if header is not None:
            yield header, "".join(seq_lines)


def rewrite_headers(input_dir, output_fasta):
    """
    :return: (num_loci, num_alleles) written
    :rtype: tuple[int, int]
    """
    num_loci = 0
    num_alleles = 0
    with open(output_fasta, 'w') as out_f:
        for locus_fasta in iter_locus_fasta_files(input_dir):
            locus_id = os.path.splitext(os.path.basename(locus_fasta))[0]
            num_loci += 1
            for allele_id, sequence in iter_fasta_records(locus_fasta):
                allele_hash = allele_calling.hash_sequence(sequence)
                out_f.write(f">{locus_id}_{allele_id}_{allele_hash}\n")
                out_f.write(sequence + "\n")
                num_alleles += 1

    return num_loci, num_alleles


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', required=True, help='Directory of per-locus FASTA files (e.g. extracted from a cgmlst.org alleles/ download)')
    parser.add_argument('--output-fasta', required=True)
    args = parser.parse_args()

    num_loci, num_alleles = rewrite_headers(args.input_dir, args.output_fasta)
    print(f"Wrote {num_alleles} alleles across {num_loci} loci: {args.output_fasta}")


if __name__ == '__main__':
    main()
