import datetime
import os
import sys
import logging

from . import parsers
from . import utils


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


def build_blastn_command(params):
    """
    Build the blastn alignment command line for assembly/contigs input.

    kma's own assembly-genefinding preset (-asm) was tried and rejected:
    empirically, against a synthetic scheme it missed 10-15% of loci that
    were exact matches to a catalogued allele, for reasons that couldn't be
    resolved by relaxing -p or -mrs. blastn found every exact match
    reliably in the same tests.

    :param params: Dictionary of parameters. Keys are 'scheme', 'tmpdir', 'assembly'.
    :type params: dict
    :return: Alignment command line
    :rtype: list
    """
    return [
        "blastn",
        "-query", params['assembly'],
        "-db", params['scheme'],
        "-outfmt", parsers.BLAST_OUTFMT,
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
    is_assembly = bool(alignment_params.get('assembly'))

    if is_assembly:
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
    alignment_elapsed_time_seconds_str = str(round(alignment_elapsed_time.total_seconds(), 2))
    logging.info(f"Alignment completed. Elapsed time: {alignment_elapsed_time_seconds_str} seconds.")
