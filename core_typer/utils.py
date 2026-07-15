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

    if args.scheme is None or args.outdir is None or (not has_reads and not has_assembly):
        parser.print_help()
        sys.exit(1)

    errors = []

    if has_reads and has_assembly:
        errors.append("Provide either --R1/--R2 (reads) or --assembly, not both.")
    elif has_reads and (args.R1 is None or args.R2 is None):
        errors.append("Both --R1 and --R2 are required when typing from reads.")
    elif has_reads:
        for label, path in [('R1', args.R1), ('R2', args.R2)]:
            if not os.path.isfile(path):
                errors.append(f"{label} file not found: {path}")
    elif has_assembly:
        if not os.path.isfile(args.assembly):
            errors.append(f"assembly file not found: {args.assembly}")

    if has_assembly:
        if shutil.which('blastn') is None:
            errors.append("blastn executable not found on PATH. Install blast (e.g. conda install -c bioconda blast) and ensure it is on PATH.")
        for scheme_ext in ['.nin', '.nsq']:
            scheme_file = f"{args.scheme}{scheme_ext}"
            if not os.path.isfile(scheme_file):
                errors.append(f"Missing expected blast-indexed scheme file: {scheme_file}. Run makeblastdb -in <scheme.fasta> -dbtype nucl -out {args.scheme}")
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
