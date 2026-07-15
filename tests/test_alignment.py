from core_typer import alignment


def test_build_kma_command_reads_mode():
    params = {'threads': 4, 'scheme': '/path/to/scheme_db', 'tmpdir': '/tmp/work', 'R1': 'R1.fastq', 'R2': 'R2.fastq'}

    command = alignment.build_kma_command(params)

    assert command[0] == 'kma'
    assert '-ipe' in command
    assert 'R1.fastq' in command
    assert 'R2.fastq' in command
    assert '-boot' in command


def test_build_blastn_command_assembly_mode():
    params = {'scheme': '/path/to/scheme_db', 'tmpdir': '/tmp/work', 'assembly': 'assembly.fasta'}

    command = alignment.build_blastn_command(params)

    assert command[0] == 'blastn'
    assert '-query' in command
    assert 'assembly.fasta' in command
    assert '-db' in command
    assert '/path/to/scheme_db' in command


def test_run_alignment_dispatches_to_blastn_for_assembly(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(alignment, 'build_blastn_command', lambda params: calls.append('blastn') or ['blastn', '-query', params['assembly']])
    monkeypatch.setattr(alignment, 'build_kma_command', lambda params: calls.append('kma') or ['kma'])

    def fake_run_command(command):
        # blast-out.tsv is the file run_alignment checks for afterward
        (tmp_path / "blast-out.tsv").write_text("")
        return None

    monkeypatch.setattr(alignment.utils, 'run_command', fake_run_command)

    alignment.run_alignment({'assembly': 'assembly.fasta', 'scheme': 'scheme_db', 'tmpdir': str(tmp_path)})

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
