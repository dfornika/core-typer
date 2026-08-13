import glob
import importlib.util
import json
import logging
import os
import shutil
import subprocess
import sys

def validate_args(args, parser):
    """
    Validate the arguments, including pre-flight checks that the kma
    executable is available and that input files exist, so failures are
    reported clearly instead of surfacing as a cryptic downstream error.
    """
    has_reads = args.R1 is not None or args.R2 is not None
    has_assembly = args.assembly is not None
    has_cds = args.cds is not None
    # Assembly and CDS input both take the blast path; assembly just has CDS
    # extracted first.
    has_blast_input = has_assembly or has_cds

    if args.scheme is None or args.outdir is None or (not has_reads and not has_blast_input):
        parser.print_help()
        sys.exit(1)

    errors = []

    if sum([has_reads, has_assembly, has_cds]) > 1:
        errors.append("Provide exactly one input type: --R1/--R2 (reads), --assembly, or --cds.")
    elif has_reads and (args.R1 is None or args.R2 is None):
        errors.append("Both --R1 and --R2 are required when typing from reads.")
    elif has_reads:
        for label, path in [('R1', args.R1), ('R2', args.R2)]:
            if not os.path.isfile(path):
                errors.append(f"{label} file not found: {path}")
    elif has_assembly:
        if not os.path.isfile(args.assembly):
            errors.append(f"assembly file not found: {args.assembly}")
    elif has_cds:
        if not os.path.isfile(args.cds):
            errors.append(f"cds file not found: {args.cds}")

    if has_assembly and importlib.util.find_spec('pyrodigal') is None:
        errors.append("pyrodigal is required for --assembly (CDS extraction). Install it (pip install pyrodigal) or supply pre-extracted CDS with --cds.")

    if has_blast_input:
        if shutil.which('blastn') is None:
            errors.append("blastn executable not found on PATH. Install blast (e.g. conda install -c bioconda blast) and ensure it is on PATH.")
        # A blast nucleotide db can be single-volume (<scheme>.nin/.nsq),
        # or - for any db over ~1GB, which includes essentially every real
        # cgMLST scheme - multi-volume: individual volumes (<scheme>.00.nin,
        # <scheme>.01.nin, ...) plus a <scheme>.nal alias file that blastn
        # resolves from the same -db <scheme> argument. Accept any of these.
        single_volume = os.path.isfile(f"{args.scheme}.nin") and os.path.isfile(f"{args.scheme}.nsq")
        multi_volume = os.path.isfile(f"{args.scheme}.nal") or bool(glob.glob(f"{args.scheme}.[0-9][0-9].nin"))
        if not (single_volume or multi_volume):
            errors.append(
                f"Missing blast-indexed scheme files for {args.scheme} (expected "
                f"{args.scheme}.nin/.nsq, a {args.scheme}.nal alias, or {args.scheme}.NN.nin "
                f"volumes). Run makeblastdb -in <scheme.fasta> -dbtype nucl -out {args.scheme}"
            )
        if not os.path.isfile(f"{args.scheme}.fasta"):
            errors.append(
                f"Missing {args.scheme}.fasta - for --assembly runs, the combined scheme FASTA passed to "
                f"makeblastdb must also be kept alongside the blast db, named <scheme>.fasta, so the full "
                f"locus list can be read directly from it."
            )
    else:
        if shutil.which('kma') is None:
            errors.append("kma executable not found on PATH. Install kma (https://bitbucket.org/genomicepidemiology/kma) and ensure it is on PATH.")
        for scheme_ext in ['.name', '.comp.b']:
            scheme_file = f"{args.scheme}{scheme_ext}"
            if not os.path.isfile(scheme_file):
                errors.append(f"Missing expected kma-indexed scheme file: {scheme_file}")

    if errors:
        for error in errors:
            logging.error(error)
        sys.exit(1)

    return args


def run_command(command):
    """
    Runs a command as a subprocess, and returns the result.

    :param command: The command to run, as a list of args (no shell involved)
    :type command: list
    :return: The result of the command
    :rtype: subprocess.CompletedProcess
    """
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        logging.error(json.dumps({
            "event_type": "command_failed",
            "command": command,
            "returncode": e.returncode,
            "stdout": e.stdout,
            "stderr": e.stderr,
        }))
        sys.exit(-1)
    except FileNotFoundError as e:
        logging.error(json.dumps({
            "event_type": "command_not_found",
            "command": command,
            "error": str(e),
        }))
        sys.exit(-1)

    return result
