import logging


def read_fasta(fasta_path):
    """
    Minimal FASTA reader yielding (record_id, sequence) tuples. record_id is
    the first whitespace-delimited token of the header line (the '>' stripped).

    :param fasta_path: Path to a FASTA file
    :type fasta_path: str
    :return: Generator of (record_id, sequence) tuples
    :rtype: Iterator[tuple[str, str]]
    """
    record_id = None
    seq_chunks = []
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if record_id is not None:
                    yield record_id, ''.join(seq_chunks)
                record_id = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
    if record_id is not None:
        yield record_id, ''.join(seq_chunks)


def extract_cds(assembly_path, output_cds_path):
    """
    Predict coding sequences (CDS) in an assembly with pyrodigal and write the
    nucleotide CDS to output_cds_path as FASTA.

    core-typer's blast-based typing aligns query records against the scheme's
    per-allele database. Because a single cgMLST locus can have thousands of
    near-identical alleles in that database, aligning whole contigs (each
    spanning hundreds of genes) lets one locus's alleles crowd out every other
    locus on the contig under blastn's -max_target_seqs cap. Extracting CDS
    first - one query record per gene, matching how tools like locidex operate
    - keeps each query's hits confined to its own locus's alleles.

    A single-genome training pass is used when there is enough sequence for
    prodigal to train (its usual >20kb heuristic); otherwise, or if training
    fails, pyrodigal's metagenomic mode is used as a fallback.

    :param assembly_path: Path to the assembly/contigs FASTA
    :type assembly_path: str
    :param output_cds_path: Path to write the extracted nucleotide CDS FASTA
    :type output_cds_path: str
    :return: The number of CDS written
    :rtype: int
    """
    import pyrodigal

    contigs = list(read_fasta(assembly_path))
    if not contigs:
        raise ValueError(f"No sequences found in assembly: {assembly_path}")

    total_length = sum(len(seq) for _, seq in contigs)
    gene_finder = pyrodigal.GeneFinder()
    trained = False
    if total_length >= 20000:
        try:
            gene_finder.train(*[seq for _, seq in contigs])
            trained = True
        except (ValueError, RuntimeError) as e:
            logging.warning(f"pyrodigal single-genome training failed ({e}); falling back to metagenomic mode.")
    if not trained:
        logging.info("Using pyrodigal metagenomic mode for CDS prediction.")
        gene_finder = pyrodigal.GeneFinder(meta=True)

    num_cds = 0
    with open(output_cds_path, 'w') as out:
        for contig_id, seq in contigs:
            genes = gene_finder.find_genes(seq)
            genes.write_genes(out, sequence_id=contig_id)
            num_cds += len(genes)

    return num_cds
