"""Misc utilities"""
import os
import sys
import csv
import glob
import logging
import subprocess

import pbcore.io.FastaIO as fi


def setup_log(alog, level=logging.INFO, file_name=None, log_filter=None,
              str_formatter='[%(levelname)s] %(asctime)-15s '
                            '[%(name)s %(funcName)s %(lineno)d] '
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
    output = "{x}_cleaned.fa".format(
        x=os.path.basename(fastafile).split('.', 1)[0])
    with open(output, 'w') as outfile:
        for record in reads:
            if record.sequence:
                outfile.write('>{r}\n{s}\n'.format(
                    r=record.header, s=record.sequence))
            else:
                log.info("Dropped!: %s", record.header)

    return output


def explode_fasta(fasta, log):
    """split input fasta into individuals"""
    in_file = False
    outdir = "fastas"

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    with open(fasta, 'r') as infile:

        for line in infile:
            if line.startswith(">"):
                if in_file:
                    outfile.close()

                fname = line.rstrip().partition(">")[2].rstrip('|arrow')
                fname = "{s}.fasta".format(s=fname)
                fout = os.path.join(outdir, fname)
                outfile = open(fout, "w")
                in_file = True

                outfile.write(line)
            elif in_file:
                outfile.write(line)
            else:
                log.debug("Line %r, but no previous > found ")

    return glob.glob('fastas/*')


def write_qfile(reference, contigs, length_dict, log):
    """Write qfile of homologous IDs for mummerplot"""

    qfiledir = os.path.join(os.getcwd(), 'qfiles')
    if not os.path.exists(qfiledir):
        os.mkdir(qfiledir)
    qfile = os.path.join(qfiledir, "{r}.qfile".format(
        r=os.path.basename(reference.rstrip('.fasta'))))

    with open(qfile, 'w') as fout:

        csv_out = csv.writer(fout, delimiter=' ', lineterminator="\n")

        for contig in contigs:
            log.debug('Writing qfile for contig: %s', contig)
            name = contig[0].rstrip('|arrow')
            length = length_dict[name]
            csv_out.writerow((contig[0], length, '+'))

    return qfile
