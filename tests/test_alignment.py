from core_typer import alignment
from core_typer import parsers


def test_prepend_blast_header_adds_header_row(tmp_path):
    blast_out = tmp_path / "blast-out.tsv"
    body = "contig1\tlocusA_1\t100.0\n"
    blast_out.write_text(body)

    alignment._prepend_blast_header(str(blast_out))

    lines = blast_out.read_text().splitlines()
    assert lines[0] == "\t".join(parsers.BLAST_OUTFMT_FIELDS)
    assert lines[1] == "contig1\tlocusA_1\t100.0"


def test_build_kma_command_reads_mode():
    params = {'threads': 4, 'scheme': '/path/to/scheme_db', 'tmpdir': '/tmp/work', 'R1': 'R1.fastq', 'R2': 'R2.fastq'}

    command = alignment.build_kma_command(params)

    assert command[0] == 'kma'
    assert '-ipe' in command
    assert 'R1.fastq' in command
    assert 'R2.fastq' in command
    assert '-boot' in command


def test_build_blastn_command_cds_mode():
    params = {'scheme': '/path/to/scheme_db', 'tmpdir': '/tmp/work', 'query': 'cds.fasta', 'threads': 4}

    command = alignment.build_blastn_command(params)

    assert command[0] == 'blastn'
    assert '-query' in command
    assert 'cds.fasta' in command
    assert '-db' in command
    assert '/path/to/scheme_db' in command
    assert '-max_target_seqs' in command
    assert str(alignment.BLAST_MAX_TARGET_SEQS) in command
    assert '-evalue' in command
    assert str(alignment.BLAST_MIN_EVALUE) in command
    assert '-num_threads' in command


def test_run_alignment_dispatches_to_blastn_for_query(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(alignment, 'build_blastn_command', lambda params: calls.append('blastn') or ['blastn', '-query', params['query']])
    monkeypatch.setattr(alignment, 'build_kma_command', lambda params: calls.append('kma') or ['kma'])

    def fake_run_command(command):
        # blast-out.tsv is the file run_alignment checks for afterward
        (tmp_path / "blast-out.tsv").write_text("")
        return None

    monkeypatch.setattr(alignment.utils, 'run_command', fake_run_command)

    alignment.run_alignment({'query': 'cds.fasta', 'scheme': 'scheme_db', 'tmpdir': str(tmp_path)})

    assert calls == ['blastn']


def test_run_alignment_dispatches_to_kma_for_reads(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(alignment, 'build_blastn_command', lambda params: calls.append('blastn') or ['blastn'])
    monkeypatch.setattr(alignment, 'build_kma_command', lambda params: calls.append('kma') or ['kma'])

    def fake_run_command(command):
        (tmp_path / "kma-out.res").write_text("")
        (tmp_path / "kma-out.mapstat").write_text("")
        (tmp_path / "kma-out.fsa").write_text("")
        return None

    monkeypatch.setattr(alignment.utils, 'run_command', fake_run_command)

    alignment.run_alignment({'R1': 'R1.fastq', 'R2': 'R2.fastq', 'scheme': 'scheme_db', 'tmpdir': str(tmp_path)})

    assert calls == ['kma']
