"""Misc utilities"""
import os
import sys
import logging
import subprocess

import pbcore.io.FastaIO as fi

def setup_log(alog, level=logging.INFO, file_name=None, log_filter=None,
              str_formatter='[%(levelname)s] %(asctime)-15s ' \
                            '[%(name)s %(funcName)s %(lineno)d] ' \
                            '%(message)s'):
    """Core Util to setup log handler"""
    alog.setLevel(logging.DEBUG)
    if file_name is None:
        handler = logging.StreamHandler(sys.stdout)
    else:
        handler = logging.FileHandler(file_name)
    formatter = logging.Formatter(str_formatter)
    handler.setFormatter(formatter)
    handler.setLevel(level)
    if log_filter:
        handler.addFilter(log_filter)
    alog.addHandler(handler)


def run(cmd, cwd, log):
    """Run shell process"""
    log.debug("Running cmd %s", cmd)

    process = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=cwd)
    stdout, stderr = process.communicate()

    if stderr:
        log.debug(stderr)

    return stdout.rstrip(), stderr


def clean_fasta(fastafile, log):
    """Check fasta for 0 length sequences / blank lines and clean it up"""
    log.info("Cleaning fasta: %s", fastafile)
    reads = fi.FastaReader(fastafile)
    output = "{x}_cleaned.fa".format(x=os.path.basename(fastafile).split('.', 1)[0])
#    output = "{x}_cleaned.fa".format(x=os.split('.', os.path.basename(fastafile))[0])
    with open(output, 'w') as outfile:
        for record in reads:
            if record.sequence:

            #            if len(record.sequence) > 0:

                outfile.write('>{r}\n{s}'.format(r=record.header, s=record.sequence))
            else:
                log.info("Dropped!: %s", record.header)

    return output
