import datetime
import os
import shutil
import sys
import logging

from . import parsers
from . import utils


def _prepend_blast_header(blast_out_path):
    """
    Prepend a tab-separated column header (BLAST_OUTFMT_FIELDS) to blastn's
    tabular output. blastn's outfmt 6 writes no header, which makes the file
    hard to read; the header is streamed on without loading the whole file
    into memory. parse_blast_result skips this header row when reading.
    """
    header = "\t".join(parsers.BLAST_OUTFMT_FIELDS) + "\n"
    tmp_path = blast_out_path + ".tmp"
    with open(tmp_path, 'w') as out, open(blast_out_path, 'r') as body:
        out.write(header)
        shutil.copyfileobj(body, out)
    os.replace(tmp_path, blast_out_path)


def build_kma_command(params):
    """
    Build the kma alignment command line for paired-end read input.

    :param params: Dictionary of parameters. Keys are 'threads', 'scheme',
        'tmpdir', 'R1', 'R2'.
    :type params: dict
    :return: Alignment command line
    :rtype: list
    """
    return [
        "kma",
        "-t", str(params['threads']),
        "-ef",
        "-cge",
        "-1t1",
        "-mem_mode",
        "-and",
        "-boot",
        "-t_db", params['scheme'],
        "-ipe", params['R1'], params['R2'],
        # kma (>=1.6) rejects -tmp with "Invalid output directory specified"
        # unless the path has a trailing separator.
        "-tmp", os.path.join(params['tmpdir'], "kma-tmp") + os.sep,
        "-o", os.path.join(params['tmpdir'], "kma-out"),
    ]


# blastn's -max_target_seqs caps how many subject sequences it keeps per query,
# and -evalue filters which alignments qualify in the first place (only hits
# with evalue <= threshold are reported). They compose: -evalue decides what is
# significant, -max_target_seqs decides how many of the significant hits to keep
# (best-scoring first).
#
# For per-CDS typing we only need each query's best-matching allele (plus a
# little headroom to resolve near-ties), not every allele of the locus. A cap
# of 50 was chosen empirically against a real S. enterica assembly + the
# senterica scheme: 50 recovers the same set of loci as an effectively
# unlimited cap (50000) while reliably returning the true best allele, whereas
# a cap of 10 (as some tools use) occasionally dropped the best allele for a
# locus with thousands of near-identical alleles. Raising the cap far higher
# only adds a handful of borderline paralogous hits at a large output-size cost
# (3M+ rows vs ~150k). -evalue 1e-4 removes low-significance chance alignments.
BLAST_MAX_TARGET_SEQS = 50
BLAST_MIN_EVALUE = 0.0001


def build_blastn_command(params):
    """
    Build the blastn alignment command line for CDS input.

    The query is a FASTA of coding sequences - one record per gene - extracted
    from the assembly (see annotation.extract_cds) or supplied directly via
    --cds. Whole contigs are deliberately not queried: see BLAST_MAX_TARGET_SEQS
    and annotation.extract_cds for why per-CDS querying is required.

    kma's own assembly-genefinding preset (-asm) was tried and rejected:
    empirically, against a synthetic scheme it missed 10-15% of loci that
    were exact matches to a catalogued allele, for reasons that couldn't be
    resolved by relaxing -p or -mrs. blastn found every exact match
    reliably in the same tests.

    :param params: Dictionary of parameters. Keys are 'scheme', 'tmpdir', 'query', 'threads'.
    :type params: dict
    :return: Alignment command line
    :rtype: list
    """
    return [
        "blastn",
        "-query", params['query'],
        "-db", params['scheme'],
        "-outfmt", parsers.BLAST_OUTFMT,
        "-max_target_seqs", str(BLAST_MAX_TARGET_SEQS),
        "-evalue", str(BLAST_MIN_EVALUE),
        "-num_threads", str(params.get('threads', 1)),
        "-out", os.path.join(params['tmpdir'], "blast-out.tsv"),
    ]


def run_alignment(alignment_params):
    """
    Run the alignment: kma for paired-end reads, blastn for an assembly.

    :param alignment_params: Dictionary of parameters (see build_kma_command/build_blastn_command)
    :type alignment_params: dict
    :return: None
    :rtype: None
    """
    is_blast = bool(alignment_params.get('query'))

    if is_blast:
        alignment_command = build_blastn_command(alignment_params)
        expected_alignment_result_files = {
            "blast_out": os.path.abspath(os.path.join(alignment_params['tmpdir'], "blast-out.tsv")),
        }
    else:
        os.makedirs(os.path.join(alignment_params['tmpdir'], "kma-tmp"), exist_ok=True)
        alignment_command = build_kma_command(alignment_params)
        expected_alignment_result_files = {
            "kma_res": os.path.abspath(os.path.join(alignment_params['tmpdir'], "kma-out.res")),
            "kma_mapstat": os.path.abspath(os.path.join(alignment_params['tmpdir'], "kma-out.mapstat")),
            "kma_fsa": os.path.abspath(os.path.join(alignment_params['tmpdir'], "kma-out.fsa")),
        }

    logging.info(f"Alignment started with command: {' '.join(alignment_command)}")
    alignment_start_timestamp = datetime.datetime.now()
    utils.run_command(alignment_command)
    alignment_end_timestamp = datetime.datetime.now()
    alignment_elapsed_time = alignment_end_timestamp - alignment_start_timestamp

    for alignment_result_file_type, alignment_result_file in expected_alignment_result_files.items():
        if not os.path.exists(alignment_result_file):
            logging.error(f"Alignment failed. Missing output file: {alignment_result_file}")
            sys.exit(-1)

    if is_blast:
        _prepend_blast_header(expected_alignment_result_files["blast_out"])

    alignment_elapsed_time_seconds_str = str(round(alignment_elapsed_time.total_seconds(), 2))
    logging.info(f"Alignment completed. Elapsed time: {alignment_elapsed_time_seconds_str} seconds.")
