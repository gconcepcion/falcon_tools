#!/usr/bin/env python
"""Check fasta for 0 length sequences / blank lines and clean it up"""
import os
import sys
import logging
import argparse

from falcon_tools import utils

log = logging.getLogger(__name__)
log.setLevel(logging.INFO)


def main():
    """Clean a fasta file"""
    args = get_parser()
    fastafile = os.path.abspath(args.fastafile)
    debug = args.debug
    logfile = args.log

    if debug:
        utils.setup_log(log, file_name=logfile, level=logging.DEBUG)
    else:
        utils.setup_log(log, file_name=logfile, level=logging.INFO)

    cleaned = utils.clean_fasta(fastafile, log)
    log.info("Your cleaned fasta file can be found here: %s", cleaned)

    return


def get_parser():
    """Return an argparse instance"""

    __version__ = 0.1
    parser = argparse.ArgumentParser(version=__version__)
    parser.add_argument("fastafile", type=str, help='path to a Fasta File')
    parser.add_argument("--log", type=str, default=None)
    parser.add_argument('--debug', action='store_true',
                        help="Print debug logging to stdout")

    return parser.parse_args()


if __name__ == "__main__":
    sys.exit(main())
