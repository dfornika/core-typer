"""
End-to-end tests against a small synthetic cgMLST scheme (see
scripts/synthetic_scheme.py). These require kma, kma_index, and wgsim on
PATH, and are skipped otherwise - they're meant to activate once a full
core-typer environment (not just the unit-test dependencies) is available.
"""
import csv
import os
import shutil
import subprocess
import sys

import pytest

import simulate_reads
import synthetic_scheme

KMA_INDEX = shutil.which("kma_index")
KMA = shutil.which("kma")
WGSIM = shutil.which("wgsim")

pytestmark = pytest.mark.skipif(
    not (KMA_INDEX and KMA and WGSIM),
    reason="kma, kma_index, and wgsim must all be on PATH for end-to-end integration tests",
)


def _run_core_typer(r1, r2, scheme_db, outdir, tmpdir, min_identity=100.0, min_coverage=100.0):
    subprocess.run(
        [
            sys.executable, "-m", "core_typer",
            "--R1", r1,
            "--R2", r2,
            "--scheme", scheme_db,
            "--outdir", outdir,
            "--tmpdir", tmpdir,
            "--min-identity", str(min_identity),
            "--min-coverage", str(min_coverage),
            "--no-cleanup",
        ],
        check=True, capture_output=True, text=True,
    )


def _read_profile(outdir):
    with open(os.path.join(outdir, "allele_profile.csv")) as f:
        rows = list(csv.reader(f))
    return dict(zip(rows[0], rows[1]))


@pytest.fixture
def synthetic_scheme_db(tmp_path):
    scheme_dir = tmp_path / "scheme"
    scheme_records, true_genome_sequence, true_profile = synthetic_scheme.generate_scheme(
        num_loci=10, alleles_per_locus=3, locus_length=300, seed=42,
    )
    scheme_dir.mkdir()
    scheme_fasta = scheme_dir / "scheme.fasta"
    genome_fasta = scheme_dir / "true_genome.fasta"
    synthetic_scheme.write_fasta(scheme_records, str(scheme_fasta))
    synthetic_scheme.write_fasta([("synthetic_genome", true_genome_sequence)], str(genome_fasta))

    scheme_db = str(scheme_dir / "scheme_db")
    subprocess.run(["kma_index", "-i", str(scheme_fasta), "-o", scheme_db], check=True, capture_output=True, text=True)

    return scheme_db, str(genome_fasta), true_profile


def test_exact_high_depth_reads_call_true_profile(tmp_path, synthetic_scheme_db):
    scheme_db, genome_fasta, true_profile = synthetic_scheme_db

    reads_dir = tmp_path / "reads"
    r1, r2 = simulate_reads.simulate(
        genome_fasta=genome_fasta, outdir=str(reads_dir),
        depth=50, read_length=125, error_rate=0.0, mutation_rate=0.0, seed=42,
    )

    outdir = tmp_path / "out"
    _run_core_typer(r1, r2, scheme_db, str(outdir), str(tmp_path / "tmp"))

    assert _read_profile(str(outdir)) == true_profile


def test_low_depth_reads_reduce_percent_called(tmp_path, synthetic_scheme_db):
    scheme_db, genome_fasta, true_profile = synthetic_scheme_db

    reads_dir = tmp_path / "reads"
    r1, r2 = simulate_reads.simulate(
        genome_fasta=genome_fasta, outdir=str(reads_dir),
        depth=1, read_length=125, error_rate=0.02, mutation_rate=0.0, seed=42,
    )

    outdir = tmp_path / "out"
    _run_core_typer(r1, r2, scheme_db, str(outdir), str(tmp_path / "tmp"))

    with open(outdir / "qc.csv") as f:
        qc_row = list(csv.DictReader(f))[0]

    assert float(qc_row["percent_called"]) < 100.0
