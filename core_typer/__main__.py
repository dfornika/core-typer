#!/usr/bin/env python3

import argparse
import datetime
import logging
import os
import shutil

from . import __version__
from . import alignment
from . import annotation
from . import allele_calling
from . import qc
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
    parser.add_argument('--multicopy-depth-ratio', type=float, default=allele_calling.MULTICOPY_MIN_DEPTH_RATIO,
                        help='For read input, minimum depth of a secondary hit relative to a locus\'s best hit for the locus to be flagged as possible multi-copy (default: %(default)s)')
    parser.add_argument('--R1', help='Read 1')
    parser.add_argument('--R2', help='Read 2')
    parser.add_argument('--assembly', help='Assembly/contigs fasta; CDS are extracted with pyrodigal before typing (alternative to --R1/--R2)')
    parser.add_argument('--cds', help='Pre-extracted CDS fasta, one record per gene, e.g. from prodigal/prokka/bakta (alternative to --assembly)')
    parser.add_argument('--scheme', help='cgMLST scheme')
    parser.add_argument('--tmpdir', default='./tmp', help='Temporary directory (default: ./tmp)')
    parser.add_argument('--no-cleanup', action='store_true', help='Do not cleanup temporary directory')
    parser.add_argument('--log-level', default='info', help='Log level (default: info)')
    parser.add_argument('--outdir', help='Output directory')
    args = parser.parse_args()

    config.configure_logging({'log_level': args.log_level})

    args = utils.validate_args(args, parser)

    now_str = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    analysis_tmpdir = os.path.join(args.tmpdir, f"{now_str}-core-typer-tmp")
    if not os.path.exists(analysis_tmpdir):
        os.makedirs(analysis_tmpdir)

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    is_blast = bool(args.assembly or args.cds)

    alignment_params = {
        'threads': args.threads,
        'scheme': args.scheme,
        'tmpdir': analysis_tmpdir,
    }
    if is_blast:
        if args.cds:
            query_cds_file = args.cds
        else:
            query_cds_file = os.path.join(analysis_tmpdir, "cds.fasta")
            logging.info(f"Extracting CDS from assembly with pyrodigal: {args.assembly}")
            num_cds = annotation.extract_cds(args.assembly, query_cds_file)
            logging.info(f"Extracted {num_cds} CDS: {query_cds_file}")
        alignment_params['query'] = query_cds_file
    else:
        alignment_params['R1'] = args.R1
        alignment_params['R2'] = args.R2

    alignment.run_alignment(alignment_params)

    if is_blast:
        blast_result_file = os.path.join(analysis_tmpdir, "blast-out.tsv")
        logging.debug(f"Parsing blast result file: {blast_result_file}")
        parsed_alignment_result = parsers.parse_blast_result(blast_result_file)
        logging.debug(f"Parsing blast result file completed: {blast_result_file}")
        # Persist the raw blast hits (header included) to the output dir so they
        # remain reviewable after the tmp directory is cleaned up.
        blast_hits_file = os.path.join(args.outdir, "blast_hits.tsv")
        logging.info(f"Writing blast hits: {blast_hits_file}")
        shutil.copyfile(blast_result_file, blast_hits_file)
        locus_ids = parsers.parse_locus_names_from_fasta(f"{args.scheme}.fasta")
        # Novel-allele extraction from an assembly hit's own matched region
        # isn't implemented yet - only reads-derived (kma -ef) consensus is
        # currently supported.
        consensus_by_template = {}
    else:
        kma_result_file = os.path.join(analysis_tmpdir, "kma-out.res")
        logging.debug(f"Parsing kma result file: {kma_result_file}")
        parsed_alignment_result = parsers.parse_kma_result(kma_result_file)
        logging.debug(f"Parsing kma result file completed: {kma_result_file}")
        locus_ids = parsers.parse_locus_names(f"{args.scheme}.name")
        kma_fsa_file = os.path.join(analysis_tmpdir, "kma-out.fsa")
        logging.debug(f"Parsing kma consensus fasta file: {kma_fsa_file}")
        consensus_by_template = parsers.parse_kma_consensus_fasta(kma_fsa_file)
        logging.debug(f"Parsing kma consensus fasta file completed: {kma_fsa_file}")

    allele_calls = allele_calling.build_complete_allele_calls(
        locus_ids,
        parsed_alignment_result,
        min_identity=args.min_identity,
        min_coverage=args.min_coverage,
    )

    novel_alleles = allele_calling.extract_novel_alleles(allele_calls, consensus_by_template)
    novel_alleles_file = os.path.join(args.outdir, "novel_alleles.fsa")
    logging.info(f"Writing {len(novel_alleles)} candidate novel allele(s): {novel_alleles_file}")
    allele_calling.write_novel_alleles(novel_alleles, novel_alleles_file)
    logging.debug(f"Writing candidate novel alleles completed: {novel_alleles_file}")

    possible_multicopy_loci = allele_calling.find_possible_multicopy_loci(parsed_alignment_result, min_depth_ratio=args.multicopy_depth_ratio)
    multicopy_loci_file = os.path.join(args.outdir, "possible_multicopy_loci.csv")
    num_multicopy_loci = len({record["locus_id"] for record in possible_multicopy_loci})
    logging.info(f"Writing {num_multicopy_loci} possible multicopy locus/loci: {multicopy_loci_file}")
    allele_calling.write_possible_multicopy_loci(possible_multicopy_loci, multicopy_loci_file)
    logging.debug(f"Writing possible multicopy loci completed: {multicopy_loci_file}")

    qc_stats = qc.calculate_qc_stats(allele_calls, possible_multicopy_loci)

    qc_stats_file = os.path.join(args.outdir, 'qc.csv')
    logging.info(f"Writing QC stats: {qc_stats_file}")
    qc.write_qc_stats(qc_stats, qc_stats_file)
    logging.debug(f"Writing QC stats completed: {qc_stats_file}")

    allele_calls_file = os.path.join(args.outdir, "allele_calls.csv")
    logging.info(f"Writing allele calls: {allele_calls_file}")
    allele_calling.write_allele_calls(allele_calls_file, allele_calls)
    logging.debug(f"Writing allele calls completed: {allele_calls_file}")

    allele_profile_file = os.path.join(args.outdir, "allele_profile.csv")
    logging.info(f"Writing allele profile: {allele_profile_file}")
    allele_calling.write_allele_profile(allele_calls, allele_profile_file)
    logging.debug(f"Writing allele profile completed: {allele_profile_file}")
    
    if not args.no_cleanup:
        shutil.rmtree(analysis_tmpdir)
        logging.info(f"Deleted tmp directory: {analysis_tmpdir}")
    else:
        logging.info(f"Skipped deleting tmp directory: {analysis_tmpdir}")


if __name__ == '__main__':
    main()
