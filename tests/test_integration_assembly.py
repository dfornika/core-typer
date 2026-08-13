"""
End-to-end tests for blast-based CDS input, against the same small synthetic
cgMLST scheme used by test_integration.py (see scripts/synthetic_scheme.py).
These require blastn and makeblastdb on PATH, and are skipped otherwise.

These drive the blast path via --cds rather than --assembly: the synthetic
"genome" is a bare concatenation of allele sequences, not a genome with real
ORFs, so it exercises the blast/parse/call path directly without depending on
pyrodigal gene prediction (which --assembly would run first).
"""
import csv
import os
import shutil
import subprocess
import sys

import pytest

import synthetic_scheme

BLASTN = shutil.which("blastn")
MAKEBLASTDB = shutil.which("makeblastdb")

pytestmark = pytest.mark.skipif(
    not (BLASTN and MAKEBLASTDB),
    reason="blastn and makeblastdb must both be on PATH for assembly-input integration tests",
)


def _run_core_typer(cds, scheme_db, outdir, tmpdir, min_identity=100.0, min_coverage=100.0):
    subprocess.run(
        [
            sys.executable, "-m", "core_typer",
            "--cds", cds,
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
def synthetic_scheme_blastdb(tmp_path):
    scheme_dir = tmp_path / "scheme"
    scheme_records, true_genome_sequence, true_profile = synthetic_scheme.generate_scheme(
        num_loci=10, alleles_per_locus=3, locus_length=300, seed=42,
    )
    scheme_dir.mkdir()
    scheme_fasta = scheme_dir / "scheme_db.fasta"
    genome_fasta = scheme_dir / "true_genome.fasta"
    synthetic_scheme.write_fasta(scheme_records, str(scheme_fasta))
    synthetic_scheme.write_fasta([("synthetic_genome", true_genome_sequence)], str(genome_fasta))

    scheme_db = str(scheme_dir / "scheme_db")
    subprocess.run(
        ["makeblastdb", "-in", str(scheme_fasta), "-dbtype", "nucl", "-out", scheme_db],
        check=True, capture_output=True, text=True,
    )

    return scheme_db, str(genome_fasta), true_profile


def test_assembly_exact_match_calls_true_profile(tmp_path, synthetic_scheme_blastdb):
    scheme_db, genome_fasta, true_profile = synthetic_scheme_blastdb

    outdir = tmp_path / "out"
    _run_core_typer(genome_fasta, scheme_db, str(outdir), str(tmp_path / "tmp"))

    assert _read_profile(str(outdir)) == true_profile


def test_assembly_depth_is_none(tmp_path, synthetic_scheme_blastdb):
    scheme_db, genome_fasta, true_profile = synthetic_scheme_blastdb

    outdir = tmp_path / "out"
    _run_core_typer(genome_fasta, scheme_db, str(outdir), str(tmp_path / "tmp"))

    with open(outdir / "allele_calls.csv") as f:
        rows = list(csv.DictReader(f))

    assert all(row["depth"] == "" for row in rows)
