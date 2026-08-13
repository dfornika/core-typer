#!/usr/bin/env python3
"""
Simulate paired-end reads from a genome FASTA using wgsim, at a target mean
depth, for validating core-typer's behavior at different coverage/error
levels. Requires wgsim to be on PATH.
"""

import argparse
import os
import subprocess


def read_fasta_length(fasta_path):
    length = 0
    with open(fasta_path) as f:
        for line in f:
            if not line.startswith(">"):
                length += len(line.strip())
    return length


def simulate(genome_fasta, outdir, depth=30.0, read_length=125, error_rate=0.02,
             mutation_rate=0.0, indel_fraction=0.15, seed=-1):
    """
    Run wgsim against genome_fasta to produce paired-end reads at approximately
    the requested mean depth.

    :return: (r1_path, r2_path)
    :rtype: tuple[str, str]
    """
    os.makedirs(outdir, exist_ok=True)
    genome_length = read_fasta_length(genome_fasta)
    num_read_pairs = max(1, round(depth * genome_length / (2 * read_length)))

    r1_path = os.path.join(outdir, "R1.fastq")
    r2_path = os.path.join(outdir, "R2.fastq")

    wgsim_command = [
        "wgsim",
        "-e", str(error_rate),
        "-N", str(num_read_pairs),
        "-1", str(read_length),
        "-2", str(read_length),
        "-r", str(mutation_rate),
        "-R", str(indel_fraction),
        "-S", str(seed),
        genome_fasta, r1_path, r2_path,
    ]
    subprocess.run(wgsim_command, check=True, capture_output=True, text=True)

    return r1_path, r2_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--genome-fasta', required=True)
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--depth', type=float, default=30.0, help='Target mean depth of coverage (default: 30.0)')
    parser.add_argument('--read-length', type=int, default=125)
    parser.add_argument('--error-rate', type=float, default=0.02, help='wgsim base error rate (default: 0.02)')
    parser.add_argument('--mutation-rate', type=float, default=0.0, help='wgsim rate of mutations relative to the reference (default: 0.0)')
    parser.add_argument('--indel-fraction', type=float, default=0.15, help='wgsim fraction of mutations that are indels (default: 0.15)')
    parser.add_argument('--seed', type=int, default=-1)
    args = parser.parse_args()

    r1_path, r2_path = simulate(
        genome_fasta=args.genome_fasta,
        outdir=args.outdir,
        depth=args.depth,
        read_length=args.read_length,
        error_rate=args.error_rate,
        mutation_rate=args.mutation_rate,
        indel_fraction=args.indel_fraction,
        seed=args.seed,
    )
    print(f"Wrote simulated reads: {r1_path}, {r2_path}")


if __name__ == '__main__':
    main()
