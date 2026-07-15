#!/usr/bin/env python3
"""
Pick one random allele per locus from a directory of per-locus FASTA files
(e.g. downloaded from cgmlst.org, before header rewriting) and stitch them
into a synthetic "true genome" with a known ground-truth allele profile.
Useful for validating core-typer against a real scheme without needing real
sequencing reads.
"""

import argparse
import glob
import json
import os
import random


def read_fasta_records(path):
    records = []
    header = None
    seq_lines = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_lines)))
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            records.append((header, "".join(seq_lines)))
    return records


def random_spacer(length, rng):
    return "".join(rng.choice("ACGT") for _ in range(length))


def build_true_genome(input_dir, spacer_length=50, seed=1234):
    """
    :return: (true_genome_sequence, true_profile) where true_profile maps
        locus_id -> allele_id for the alleles used in true_genome_sequence
    :rtype: tuple[str, dict[str, str]]
    """
    rng = random.Random(seed)
    locus_fastas = sorted(glob.glob(os.path.join(input_dir, "*.fasta")))
    true_profile = {}
    pieces = [random_spacer(spacer_length, rng)]
    for locus_fasta in locus_fastas:
        locus_id = os.path.splitext(os.path.basename(locus_fasta))[0]
        allele_id, allele_seq = rng.choice(read_fasta_records(locus_fasta))
        true_profile[locus_id] = allele_id
        pieces.append(allele_seq)
        pieces.append(random_spacer(spacer_length, rng))
    return "".join(pieces), true_profile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', required=True, help='Directory of per-locus FASTA files (cgmlst.org format, before header rewriting)')
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--spacer-length', type=int, default=50)
    parser.add_argument('--seed', type=int, default=1234)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    true_genome_sequence, true_profile = build_true_genome(args.input_dir, args.spacer_length, args.seed)

    genome_fasta = os.path.join(args.outdir, "true_genome.fasta")
    profile_path = os.path.join(args.outdir, "true_profile.json")
    with open(genome_fasta, 'w') as f:
        f.write(">true_genome\n")
        for i in range(0, len(true_genome_sequence), 70):
            f.write(true_genome_sequence[i:i + 70] + "\n")
    with open(profile_path, 'w') as f:
        json.dump(true_profile, f, indent=2)

    print(f"Wrote true genome ({len(true_genome_sequence)} bp): {genome_fasta}")
    print(f"Wrote true profile: {profile_path}")


if __name__ == '__main__':
    main()
