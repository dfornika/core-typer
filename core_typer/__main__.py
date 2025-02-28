#!/usr/bin/env python3

import argparse
import csv
import datetime
import json
import logging
import os
import subprocess
import sys

from . import __version__
from . import alignment
from . import allele_calling
from . import config
from . import parsers
from . import utils


def main():
    parser = argparse.ArgumentParser(description='A cgMLST Typing Tool')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of CPU threads to use (default: 1)')
    parser.add_argument('-p', '--prefix', help='Prefix for output files (default: taken from R1 fastq file name)')
    parser.add_argument('--min-identity', type=float, default=100.0, help='Minimum percent identity (default: 100.0)')
    parser.add_argument('--min-coverage', type=float, default=100.0, help='Minimum percent coverage (default: 100.0)')
    parser.add_argument('--R1', help='Read 1')
    parser.add_argument('--R2', help='Read 2')
    parser.add_argument('--scheme', help='cgMLST scheme')
    parser.add_argument('--tmpdir', default='/tmp', help='Temporary directory (default: /tmp)')
    parser.add_argument('--log-level', default='info', help='Log level (default: info)')
    parser.add_argument('--outdir', help='Output directory')
    args = parser.parse_args()

    args = utils.validate_args(args, parser)

    config.configure_logging({'log_level': args.log_level})

    if not os.path.exists(args.tmpdir):
        os.makedirs(args.tmpdir)

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    alignment_params = {
        'R1': args.R1,
        'R2': args.R2,
        'threads': args.threads,
        'scheme': args.scheme,
        'tmpdir': args.tmpdir,
    }

    alignment.run_alignment(alignment_params)

    kma_result_file = os.path.join(args.tmpdir, "kma-out.res")
    logging.debug(f"Parsing kma result file: {kma_result_file}")
    parsed_kma_result = parsers.parse_kma_result(kma_result_file)
    logging.debug(f"Parsing kma result file completed: {kma_result_file}")

    allele_calls = []
    for locus_id, kma_results in parsed_kma_result.items():
        best_allele = allele_calling.choose_best_allele(kma_results, min_identity=args.min_identity, min_coverage=args.min_coverage)
        allele_calls.append(best_allele)

    kma_mapstat_file = os.path.join(args.tmpdir, "kma-out.mapstat")
    logging.info(f"Parsing kma mapstat file: {kma_mapstat_file}")
    parsed_kma_mapstat = parsers.parse_kma_mapstat(kma_mapstat_file)
    logging.debug(f"Parsing kma mapstat file completed: {kma_mapstat_file}")

    allele_calls_file = os.path.abspath(os.path.join(args.outdir, "allele_calls.csv"))
    logging.info(f"Writing allele calls: {allele_calls_file}")
    allele_calling.write_allele_calls(allele_calls_file, allele_calls)
    logging.debug(f"Writing allele calls completed: {allele_calls_file}")
