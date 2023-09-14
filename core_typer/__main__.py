#!/usr/bin/env python3

import argparse
import csv
import json
import logging
import os
import subprocess
import sys

from . import __version__

def run_command(command):
    try:
        result = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print("Error running command: " + command)
        print(e)
        sys.exit(1)

    return result

def parse_kma_result(kma_result_file):
    """
    """
    kma_result_by_locus_id = {}
    int_fields = [
        "score",
        "expected",
        "template_length",
    ]
    float_fields = [
        "template_identity",
        "template_coverage",
        "query_identity",
        "query_coverage",
        "depth",
        "q_value",
        "p_value",
    ]
    with open(kma_result_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            record = {}
            for k, v in row.items():
                key = k.lower().replace("#", "")
                if key in int_fields:
                    try:
                        record[key] = int(v.strip())
                    except ValueError as e:
                        record[key] = None
                elif key in float_fields:
                    try:
                        record[key] = float(v.strip())
                    except ValueError as e:
                        record[key] = None
                else:
                    record[key] = v.strip()
                locus_id = record["template"].split("_")[0]
                allele_id = record["template"].split("_")[1]
                record["locus_id"] = locus_id
                record["allele_id"] = allele_id
                if locus_id not in kma_result_by_locus_id:
                    kma_result_by_locus_id[locus_id] = []
                kma_result_by_locus_id[locus_id].append(record)

    for locus_id, kma_results in kma_result_by_locus_id.items():
        kma_result_by_locus_id[locus_id] = sorted(kma_results, key=lambda k: k["score"], reverse=True)

    return kma_result_by_locus_id


def choose_best_allele(kma_results, min_identity=100.0, min_coverage=100.0):
    """
    """
    best_allele = None
    for kma_result in kma_results:
        if best_allele is None:
            best_allele = kma_result
        else:
            if kma_result["score"] > best_allele["score"]:
                if kma_result["template_identity"] < min_identity:
                    kma_result["allele_id"] = "-"
                if kma_result["template_coverage"] < min_coverage:
                    kma_result["allele_id"] = "-"
                best_allele = kma_result

    return best_allele


def parse_kma_mapstat(kma_mapstat_file):
    """
    """
    parsed_kma_mapstat_by_locus_id = {}
    header = [
        "ref_sequence",
        "read_count",
        "fragment_count",
        "map_score_sum",
        "ref_covered_positions",
        "ref_consensus_sum",
        "bp_total",
        "depth_variance",
        "nuc_high_depth_variance",
        "depth_max",
        "snp_sum",
        "insert_sum",
        "deletion_sum",
        "read_count_aln",
        "fragment_count_aln",
    ]
    int_fields = [
        "read_count",
        "fragment_count",
        "map_score_sum",
        "ref_covered_positions",
        "ref_consensus_sum",
        "bp_total",
        "nuc_high_depth_variance",
        "depth_max",
        "snp_sum",
        "insert_sum",
        "deletion_sum",
        "read_count_aln",
        "fragment_count_aln",
    ]
    float_fields = [
        "depth_variance",
    ]
    with open(kma_mapstat_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            else:
                record = dict(zip(header, line.strip().split("\t")))
                locus_id = record['ref_sequence'].split("_")[0]
                allele_id = record['ref_sequence'].split("_")[1]
                record['locus_id'] = locus_id
                record['allele_id'] = allele_id
                for field in int_fields:
                    try:
                        record[field] = int(record[field])
                    except ValueError as e:
                        record[field] = None
                for field in float_fields:
                    try:
                        record[field] = float(record[field])
                    except ValueError as e:
                        record[field] = None
                if locus_id not in parsed_kma_mapstat_by_locus_id:
                    parsed_kma_mapstat_by_locus_id[locus_id] = []
                parsed_kma_mapstat_by_locus_id[locus_id].append(record)

    for locus_id, kma_mapstat in parsed_kma_mapstat_by_locus_id.items():
        parsed_kma_mapstat_by_locus_id[locus_id] = sorted(kma_mapstat, key=lambda k: k["map_score_sum"], reverse=True)

    return parsed_kma_mapstat_by_locus_id


def write_allele_calls(allele_calls_file, allele_calls):
    """
    """
    output_fieldnames = [
        "locus_id",
        "allele_id",
        "score",
        "template_length",
        "percent_identity",
        "percent_coverage",
        "depth",
    ]
    with open(allele_calls_file, 'w') as f:
        writer = csv.DictWriter(f, delimiter=',', fieldnames=output_fieldnames)
        writer.writeheader()
        for allele_call in allele_calls:
            output_record = {}
            for fieldname in output_fieldnames:
                if fieldname in allele_call:
                    output_record[fieldname] = allele_call[fieldname]
                elif fieldname == "percent_identity":
                    output_record[fieldname] = allele_call["template_identity"]
                elif fieldname == "percent_coverage":
                    output_record[fieldname] = allele_call["template_coverage"]
            writer.writerow(output_record)


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

    if args.R1 is None or args.R2 is None or args.scheme is None or args.outdir is None:
        parser.print_help()
        sys.exit(1)

    try:
        log_level = getattr(logging, args.log_level.upper())
    except AttributeError as e:
        log_level = logging.INFO

    logging.basicConfig(
        format='{"timestamp": "%(asctime)s.%(msecs)03d", "level": "%(levelname)s", "message": %(message)s}',
        datefmt='%Y-%m-%dT%H:%M:%S',
        encoding='utf-8',
        level=log_level,
    )
    logging.debug(json.dumps({"event_type": "debug_logging_enabled"}))
    
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    # Run the pipeline
    kma_command = [
        "kma",
        "-t", str(args.threads),
        "-ef",
        "-cge",
        "-boot",
        "-1t1",
        "-mem_mode",
        "-and",
        "-t_db", args.scheme,
        "-ipe", args.R1, args.R2,
        "-tmp", os.path.join(args.tmpdir, "kma-tmp"),
        "-o", os.path.join(args.tmpdir, "kma-out"),
    ]

    logging.info(json.dumps({
        "event_type": "kma_analysis_started",
        "command": " ".join(kma_command),
    }))
    kma_result = run_command(" ".join(kma_command))
    logging.info(json.dumps({
        "event_type": "kma_analysis_completed",
    }))

    kma_result_file = os.path.join(args.tmpdir, "kma-out.res")
    logging.info(json.dumps({
        "event_type": "parse_kma_result_started",
        "kma_result_file": kma_result_file,
    }))
    parsed_kma_result = parse_kma_result(kma_result_file)
    logging.info(json.dumps({
        "event_type": "parse_kma_result_completed",
        "kma_result_file": kma_result_file,
    }))

    allele_calls = []
    for locus_id, kma_results in parsed_kma_result.items():
        best_allele = choose_best_allele(kma_results, min_identity=args.min_identity, min_coverage=args.min_coverage)
        allele_calls.append(best_allele)

    logging.info(json.dumps({
        "event_type": "parse_kma_mapstat_started",
    }))
    parsed_kma_mapstat = parse_kma_mapstat(os.path.join(args.tmpdir, "kma-out.mapstat"))
    logging.info(json.dumps({
        "event_type": "parse_kma_mapstat_completed",
    }))

    allele_calls_file = os.path.abspath(os.path.join(args.outdir, "allele_calls.csv"))
    logging.info(json.dumps({
        "event_type": "write_allele_calls_started",
        "allele_calls_file": allele_calls_file,
    }))
    write_allele_calls(allele_calls_file, allele_calls)
    logging.info(json.dumps({
        "event_type": "write_allele_calls_completed",
        "allele_calls_file": allele_calls_file,
    }))
