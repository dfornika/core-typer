import json
import logging
import subprocess
import sys

def validate_args(args, parser):
    """
    Validate the arguments.
    
    """
    conditions = [
        args.R1 is not None,
        args.R2 is not None,
        args.scheme is not None,
        args.outdir is not None,
    ]
    all_conditions_met = all(conditions)
    if not all_conditions_met:
        parser.print_help()
        sys.exit(1)
    else:
        return args


def run_command(command):
    """
    Runs a command as a subprocess, and returns the result.

    :param command: The command to run
    :type command: str
    :return: The result of the command
    :rtype: subprocess.CompletedProcess
    """
    try:
        result = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        logging.error(json.dumps({
            "event_type": "command_failed",
            "command": command,
            "returncode": e.returncode,
            "stdout": e.stdout,
            "stderr": e.stderr,
        }))
        sys.exit(-1)

    return result
