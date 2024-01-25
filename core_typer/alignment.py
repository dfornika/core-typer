import datetime
import json
import os
import logging

from . import utils

def build_alignment_command(params):
    """
    Build the kma alignment command line.

    :param params: Dictionary of parameters. Keys of params are: 'threads', 'scheme', 'R1', 'R2', 'tmpdir'
    :type params: dict
    :return: Alignment command line
    :rtype: list
    """
    kma_command = [
        "kma",
        "-t", str(params['threads']),
        "-ef",
        "-cge",
        "-boot",
        "-1t1",
        "-mem_mode",
        "-and",
        "-t_db", params['scheme'],
        "-ipe", params['R1'], params['R2'],
        "-tmp", os.path.join(params['tmpdir'], "kma-tmp"),
        "-o", os.path.join(params['tmpdir'], "kma-out"),
    ]

    return kma_command


def run_alignment(alignment_params):
    """
    Run the alignment.

    :param alignment_command: The alignment command
    :type alignment_command: list
    :return: None
    :rtype: None
    """
    alignment_command = build_alignment_command(alignment_params)
    logging.info(json.dumps({
        "event_type": "alignment_started",
        "alignment_command": " ".join(alignment_command),
    }))
    alignment_start_timestamp = datetime.datetime.now()
    alignment_result = utils.run_command(" ".join(alignment_command))
    alignment_end_timestamp = datetime.datetime.now()
    alignment_elapsed_time = alignment_end_timestamp - alignment_start_timestamp
    expected_alignment_result_files = {
        "kma_res": os.path.abspath(os.path.join(alignment_params['tmpdir'], "kma-out.res")),
        "kma_mapstat": os.path.abspath(os.path.join(alignment_params['tmpdir'], "kma-out.mapstat")),
    }
    for alignment_result_file_type, alignment_result_file in expected_alignment_result_files.items():
        if not os.path.exists(alignment_result_file):
            logging.error(json.dumps({
                "event_type": "alignment_failed",
                "alignment_result_file_type": alignment_result_file_type,
                "alignment_result_file": alignment_result_file,
            }))
            sys.exit(-1)
    logging.info(json.dumps({
        "event_type": "alignment_completed",
        "elapsed_time_seconds": str(alignment_elapsed_time.total_seconds()),
    }))
    
