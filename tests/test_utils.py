import argparse

import pytest

from core_typer import utils


class DummyParser:
    def print_help(self):
        pass


def make_args(**overrides):
    defaults = {'R1': None, 'R2': None, 'assembly': None, 'cds': None, 'scheme': None, 'outdir': None}
    defaults.update(overrides)
    return argparse.Namespace(**defaults)


@pytest.fixture
def scheme_prefix(tmp_path, monkeypatch):
    monkeypatch.setattr(utils.shutil, 'which', lambda name: f'/usr/bin/{name}')
    prefix = tmp_path / "scheme"
    (tmp_path / "scheme.name").write_text("")
    (tmp_path / "scheme.comp.b").write_text("")
    return str(prefix)


@pytest.fixture
def blast_scheme_prefix(tmp_path, monkeypatch):
    monkeypatch.setattr(utils.shutil, 'which', lambda name: f'/usr/bin/{name}')
    prefix = tmp_path / "scheme"
    (tmp_path / "scheme.nin").write_text("")
    (tmp_path / "scheme.nsq").write_text("")
    (tmp_path / "scheme.fasta").write_text("")
    return str(prefix)


def test_validate_args_requires_some_input_mode(scheme_prefix, tmp_path):
    args = make_args(scheme=scheme_prefix, outdir=str(tmp_path))

    with pytest.raises(SystemExit):
        utils.validate_args(args, DummyParser())


def test_validate_args_rejects_both_reads_and_assembly(scheme_prefix, tmp_path):
    r1 = tmp_path / "R1.fastq"
    r2 = tmp_path / "R2.fastq"
    assembly = tmp_path / "assembly.fasta"
    for f in (r1, r2, assembly):
        f.write_text("")

    args = make_args(R1=str(r1), R2=str(r2), assembly=str(assembly), scheme=scheme_prefix, outdir=str(tmp_path))

    with pytest.raises(SystemExit):
        utils.validate_args(args, DummyParser())


def test_validate_args_rejects_partial_reads(scheme_prefix, tmp_path):
    r1 = tmp_path / "R1.fastq"
    r1.write_text("")

    args = make_args(R1=str(r1), scheme=scheme_prefix, outdir=str(tmp_path))

    with pytest.raises(SystemExit):
        utils.validate_args(args, DummyParser())


def test_validate_args_accepts_reads_only(scheme_prefix, tmp_path):
    r1 = tmp_path / "R1.fastq"
    r2 = tmp_path / "R2.fastq"
    r1.write_text("")
    r2.write_text("")

    args = make_args(R1=str(r1), R2=str(r2), scheme=scheme_prefix, outdir=str(tmp_path))

    assert utils.validate_args(args, DummyParser()) is args


def test_validate_args_accepts_assembly_only(blast_scheme_prefix, tmp_path):
    assembly = tmp_path / "assembly.fasta"
    assembly.write_text("")

    args = make_args(assembly=str(assembly), scheme=blast_scheme_prefix, outdir=str(tmp_path))

    assert utils.validate_args(args, DummyParser()) is args


def test_validate_args_rejects_missing_assembly_file(blast_scheme_prefix, tmp_path):
    args = make_args(assembly=str(tmp_path / "does_not_exist.fasta"), scheme=blast_scheme_prefix, outdir=str(tmp_path))

    with pytest.raises(SystemExit):
        utils.validate_args(args, DummyParser())


def test_validate_args_accepts_cds_only(blast_scheme_prefix, tmp_path):
    cds = tmp_path / "cds.fasta"
    cds.write_text("")

    args = make_args(cds=str(cds), scheme=blast_scheme_prefix, outdir=str(tmp_path))

    assert utils.validate_args(args, DummyParser()) is args


def test_validate_args_rejects_missing_cds_file(blast_scheme_prefix, tmp_path):
    args = make_args(cds=str(tmp_path / "does_not_exist.fasta"), scheme=blast_scheme_prefix, outdir=str(tmp_path))

    with pytest.raises(SystemExit):
        utils.validate_args(args, DummyParser())


def test_validate_args_rejects_assembly_and_cds(blast_scheme_prefix, tmp_path):
    assembly = tmp_path / "assembly.fasta"
    cds = tmp_path / "cds.fasta"
    assembly.write_text("")
    cds.write_text("")

    args = make_args(assembly=str(assembly), cds=str(cds), scheme=blast_scheme_prefix, outdir=str(tmp_path))

    with pytest.raises(SystemExit):
        utils.validate_args(args, DummyParser())


def test_validate_args_accepts_assembly_with_volume_split_blast_db(tmp_path, monkeypatch):
    monkeypatch.setattr(utils.shutil, 'which', lambda name: f'/usr/bin/{name}')
    assembly = tmp_path / "assembly.fasta"
    assembly.write_text("")
    # A large scheme's blast db is split into volumes (scheme.00.nin ...) plus
    # a scheme.nal alias, with no bare scheme.nin/.nsq. This must be accepted.
    (tmp_path / "scheme.00.nin").write_text("")
    (tmp_path / "scheme.00.nsq").write_text("")
    (tmp_path / "scheme.nal").write_text("")
    (tmp_path / "scheme.fasta").write_text("")

    args = make_args(assembly=str(assembly), scheme=str(tmp_path / "scheme"), outdir=str(tmp_path))

    assert utils.validate_args(args, DummyParser()) is args


def test_validate_args_rejects_assembly_without_blast_db(tmp_path, monkeypatch):
    monkeypatch.setattr(utils.shutil, 'which', lambda name: f'/usr/bin/{name}')
    assembly = tmp_path / "assembly.fasta"
    assembly.write_text("")
    # scheme prefix has no blast db files (.nin/.nsq) or companion .fasta
    args = make_args(assembly=str(assembly), scheme=str(tmp_path / "scheme"), outdir=str(tmp_path))

    with pytest.raises(SystemExit):
        utils.validate_args(args, DummyParser())


def test_validate_args_rejects_assembly_missing_scheme_fasta(tmp_path, monkeypatch):
    monkeypatch.setattr(utils.shutil, 'which', lambda name: f'/usr/bin/{name}')
    assembly = tmp_path / "assembly.fasta"
    assembly.write_text("")
    (tmp_path / "scheme.nin").write_text("")
    (tmp_path / "scheme.nsq").write_text("")
    # deliberately no scheme.fasta

    args = make_args(assembly=str(assembly), scheme=str(tmp_path / "scheme"), outdir=str(tmp_path))

    with pytest.raises(SystemExit):
        utils.validate_args(args, DummyParser())
