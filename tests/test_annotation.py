import importlib.util

import pytest

from core_typer import annotation


def test_read_fasta_parses_records(tmp_path):
    fasta = tmp_path / "seqs.fasta"
    fasta.write_text(">contig1 some description\nACGT\nACGT\n>contig2\nTTTT\n")

    records = list(annotation.read_fasta(str(fasta)))

    # record_id is only the first whitespace-delimited token; multi-line
    # sequence is concatenated.
    assert records == [("contig1", "ACGTACGT"), ("contig2", "TTTT")]


def test_read_fasta_empty_file(tmp_path):
    fasta = tmp_path / "empty.fasta"
    fasta.write_text("")

    assert list(annotation.read_fasta(str(fasta))) == []


def test_extract_cds_raises_on_empty_assembly(tmp_path):
    empty = tmp_path / "empty.fasta"
    empty.write_text("")

    with pytest.raises(ValueError):
        annotation.extract_cds(str(empty), str(tmp_path / "cds.fasta"))


@pytest.mark.skipif(
    importlib.util.find_spec("pyrodigal") is None,
    reason="pyrodigal must be installed for CDS extraction",
)
def test_extract_cds_writes_genes(tmp_path):
    # A short synthetic sequence with a clear ORF: ATG ... TAA in frame.
    orf = "ATG" + "AAAGCACTG" * 30 + "TAA"
    fasta = tmp_path / "genome.fasta"
    fasta.write_text(f">contig1\n{orf}\n")
    out = tmp_path / "cds.fasta"

    num_cds = annotation.extract_cds(str(fasta), str(out))

    records = list(annotation.read_fasta(str(out)))
    assert num_cds == len(records)
    assert num_cds >= 1
    # CDS records are named after the source contig.
    assert all(rid.startswith("contig1") for rid, _ in records)
