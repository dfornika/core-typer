import json
import logging


def configure_logging(params):
    """
    Configure logging.
    """
    try:
        log_level = getattr(logging, params['log_level'].upper())
    except AttributeError as e:
        log_level = logging.INFO

    logging.basicConfig(
        format='%(asctime)s.%(msecs)03d\t%(levelname)s:\t%(message)s',
        datefmt='%Y-%m-%dT%H:%M:%S',
        encoding='utf-8',
        level=log_level,
    )
    logging.debug(f"Debug logging enabled.")
