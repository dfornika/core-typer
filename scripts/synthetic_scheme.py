#!/usr/bin/env python3
"""
Generate a tiny synthetic cgMLST scheme for testing/validating core-typer,
without needing any real scheme data.

For each locus, a random reference sequence is generated, then a handful of
other alleles are derived from it by introducing a few SNPs. One allele per
locus is chosen as "true" and stitched together (with random spacer sequence
between loci) into a synthetic genome, giving a known ground-truth allele
profile to validate core-typer's calls against.
"""

import argparse
import json
import os
import random

BASES = "ACGT"


def random_sequence(length, rng):
    return "".join(rng.choice(BASES) for _ in range(length))


def mutate_sequence(sequence, num_snps, rng):
    sequence = list(sequence)
    positions = rng.sample(range(len(sequence)), min(num_snps, len(sequence)))
    for pos in positions:
        original = sequence[pos]
        sequence[pos] = rng.choice([base for base in BASES if base != original])
    return "".join(sequence)


def generate_scheme(num_loci=10, alleles_per_locus=3, locus_length=300, spacer_length=50, seed=1234):
    """
    :return: (scheme_records, true_genome_sequence, true_profile)
        scheme_records: list of (header, sequence) for the scheme fasta, header format "<locus_id>_<allele_id>"
        true_genome_sequence: single contig string containing one allele per locus, joined by random spacers
        true_profile: dict of locus_id -> allele_id, matching the alleles used in true_genome_sequence
    :rtype: tuple[list[tuple[str, str]], str, dict[str, str]]
    """
    rng = random.Random(seed)
    scheme_records = []
    true_profile = {}
    # Flanking spacer at each end so the first/last locus aren't sitting at
    # the exact edge of the (linear) genome - read simulators like wgsim
    # can't simulate reads wrapping past a linear contig's ends, so a locus
    # placed right at the edge sees thinner coverage there than a real
    # genome would give it.
    true_genome_pieces = [random_sequence(spacer_length, rng)]

    for locus_index in range(1, num_loci + 1):
        locus_id = f"locus{locus_index:03d}"
        reference = random_sequence(locus_length, rng)
        true_allele_index = rng.randint(1, alleles_per_locus)
        for allele_index in range(1, alleles_per_locus + 1):
            allele_seq = reference if allele_index == 1 else mutate_sequence(reference, rng.randint(3, 8), rng)
            scheme_records.append((f"{locus_id}_{allele_index}", allele_seq))
            if allele_index == true_allele_index:
                true_profile[locus_id] = str(allele_index)
                true_genome_pieces.append(allele_seq)
        true_genome_pieces.append(random_sequence(spacer_length, rng))

    return scheme_records, "".join(true_genome_pieces), true_profile


def write_fasta(records, path, line_length=70):
    with open(path, 'w') as f:
        for header, sequence in records:
            f.write(f">{header}\n")
            for i in range(0, len(sequence), line_length):
                f.write(sequence[i:i + line_length] + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--num-loci', type=int, default=10)
    parser.add_argument('--alleles-per-locus', type=int, default=3)
    parser.add_argument('--locus-length', type=int, default=300)
    parser.add_argument('--spacer-length', type=int, default=50)
    parser.add_argument('--seed', type=int, default=1234)
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    scheme_records, true_genome_sequence, true_profile = generate_scheme(
        num_loci=args.num_loci,
        alleles_per_locus=args.alleles_per_locus,
        locus_length=args.locus_length,
        spacer_length=args.spacer_length,
        seed=args.seed,
    )

    scheme_fasta = os.path.join(args.outdir, "scheme.fasta")
    genome_fasta = os.path.join(args.outdir, "true_genome.fasta")
    profile_path = os.path.join(args.outdir, "true_profile.json")

    write_fasta(scheme_records, scheme_fasta)
    write_fasta([("synthetic_genome", true_genome_sequence)], genome_fasta)
    with open(profile_path, 'w') as f:
        json.dump(true_profile, f, indent=2)

    print(f"Wrote scheme FASTA ({len(scheme_records)} alleles across {args.num_loci} loci): {scheme_fasta}")
    print(f"Wrote true genome FASTA ({len(true_genome_sequence)} bp): {genome_fasta}")
    print(f"Wrote true allele profile: {profile_path}")


if __name__ == '__main__':
    main()
