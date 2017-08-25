#!/usr/bin/env python
"""
Check raw read, pread and overlap distributions in a completed FALCON job_root
Author: Greg Concepcion gconcepcion@pacificbiosciences.com
"""

import os
import sys
import logging
import argparse
import subprocess

import matplotlib
from matplotlib import pyplot as plt
import pandas

matplotlib.style.use('ggplot')
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)


def run(cmd, cwd):
    """Run shell process"""
    log.debug("Running cmd %s", cmd)

    process = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=cwd)
    stdout, stderr = process.communicate()

    if stderr:
        log.debug(stderr)

    return stdout.rstrip(), stderr


def get_overlaps(jobdir, nproc):
    """Get overlap distributions"""
    log.info("Gathering overlap stats")

    lasfofn = os.path.join(jobdir, '1-preads_ovl/merge-gather/las.fofn')
    if not os.path.exists(lasfofn):
        log.debug("No las.fofn!")

    cmd = ['fc_ovlp_stats', '--n_core', str(nproc), '--fofn', lasfofn]

    cwd = os.path.join(jobdir, '2-asm-falcon')
    stdout, stderr = run(cmd, cwd)
    if stderr:
        log.debug(stderr)

    return stdout.splitlines()


def plot_ovlp_stats(jobdir, nproc):
    """Plot 5' and 3' Overlap distributions"""

    log.info("Generating overlap plots")
    overlaps = get_overlaps(jobdir, nproc)
    ovlp_dict = {}
    for ovlp in overlaps:
        rid, length, fiveprime, threeprime = ovlp.split()
        ovlp_dict[rid] = fiveprime, threeprime

    fiveprime_ovlps = [int(v[0]) for k, v in ovlp_dict.items()]
    threeprime_ovlps = [int(v[1]) for k, v in ovlp_dict.items()]
    fig, axs = plt.subplots(2, 1)

    dataframe = pandas.DataFrame({'five_prime_ovlps': fiveprime_ovlps,
                                  'three_prime_ovlps': threeprime_ovlps})

    binsfive = dataframe['five_prime_ovlps'].max() + 1
    binsthree = dataframe['three_prime_ovlps'].max() + 1

    dataframe['five_prime_ovlps'].plot.hist(
        bins=binsfive, ax=axs[0], figsize=(5, 10))
    axs[0].set_title('5\' overlaps')
    axs[0].set_xlim(0, 100)

    dataframe['three_prime_ovlps'].plot.hist(
        bins=binsthree, ax=axs[1], figsize=(5, 10))
    axs[1].set_title('3\' overlaps')
    axs[1].set_xlim(0, 100)

    outfig = os.path.join('outfigs', 'overlap_distribution.png')
    plt.savefig(outfig)
    return dataframe


def get_length_distribution(dbpath):
    "Get sequence lengths from the DB"

    log.info("Getting lengths from %s ", dbpath)
    cmd = ['DBdump', '-h', dbpath]
    cwd = os.path.dirname(os.path.dirname(dbpath))
    stdout, stderr = run(cmd, cwd)
    length_list = []
    if stderr:
        log.debug(stderr)

    for line in stdout.splitlines():
        if line.startswith('L '):
            _, _, start, end = line.split()
            seqlen = int(end) - int(start)
            length_list.append(seqlen)

    log.info("Entries in DB: %d", len(length_list))
    return length_list


def plot_length_distribution_raw(path):
    "Plot Raw Read length distribution"

    rawdb = os.path.join(path, '0-rawreads', 'raw_reads.db')
    lengths = None
    if os.path.exists(rawdb):
        lengths = get_length_distribution(rawdb)
        plot_length_distribution(lengths, context='Raw_read')
    else:
        log.info("Raw read DB not found, skipping...")
    return lengths


def plot_length_distribution_preads(path):
    "Plot Pread length distribution"

    preaddb = os.path.join(path, '1-preads_ovl', 'preads.db')
    lengths = None
    if os.path.exists(preaddb):

        lengths = get_length_distribution(preaddb)
        plot_length_distribution(lengths, context='Pread')

    else:
        log.info("Pread DB not found, skipping...")

    return lengths


def plot_length_distribution(lengths, context='None'):
    "Plot length distribution"
    title = "{t} Length Distribution".format(t=context)

    length_df1 = pandas.DataFrame({'lengths': lengths})
    ax = length_df1.plot.hist(bins=100, alpha=0.5, legend=False)
    _, ymax = tuple(ax.get_ylim())
    _, xmax = tuple(ax.get_xlim())
    mean = int(length_df1['lengths'].mean())
    plt.xlabel('Read length', fontsize=16)
    plt.ylabel('Count', fontsize=16)
    plt.title(title, fontsize=18)
    maxlen = "Max: {d} Bp".format(d=xmax)
    meanlen = "Mean: {d} Bp".format(d=mean)

    plt.text(xmax / 2, ymax / 1.5, maxlen, fontsize=14)
    plt.text(xmax / 2, ymax / 2, meanlen, fontsize=14)

    outfig = os.path.join('outfigs', "{c}_out.png".format(c=context))
    plt.savefig(outfig)


def plot_dual_lengths(rlens, plens):
    "Generate dual length distribution"
    title = "Read Length Distributions"
    length_df1 = pandas.DataFrame({'raw': rlens})
    length_df2 = pandas.DataFrame({'pread': plens})

    ax1 = length_df1.hist(bins=100, alpha=0.5)
    for ax, (colname, values) in zip(ax1.flat, length_df2.iteritems()):
        values.plot.hist(ax=ax, bins=100, alpha=0.5)
    plt.xlim(0, 60000)
    plt.legend(['Pread', 'Raw'])
    plt.title(title)
    outfig = os.path.join('outfigs', 'length_distributions.png')
    log.info("Saving figure:\n")
    log.info(outfig)
    plt.savefig(outfig)


def setup_log(alog, level=logging.INFO, file_name=None, log_filter=None,
              str_formatter='[%(levelname)s] %(asctime)-15s '
                            '[%(funcName)s %(lineno)d] '
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


def validate_falcon_root(dirpath):
    """ Validate which distributions we can create """
    raw_reads_db = os.path.join(dirpath, '0-rawreads/raw_reads.db')
    preads_db = os.path.join(dirpath, '1-preads_ovl/preads.db')
    overlaps = os.path.join(dirpath, '2-asm-falcon/preads.ovl')
    raw = False
    pread = False
    overlap = False

    if os.path.exists(raw_reads_db):
        raw = True
        log.info("raw_reads.db exists")
    else:
        log.info("No raw_reads.db found, skipping raw_read distribution!")

    if os.path.exists(preads_db):
        pread = True
        log.info("preads.db exists")
    else:
        log.info("No preads.db found, skipping pread distribution!")

    if os.path.exists(overlaps):
        log.info("preads.ovl exists")
        overlap = True
    else:
        log.info("No preads.ovl found, skipping overlap distribution!")
        
    return raw, pread, overlap


def main():
    """Generate Read length Distribution and Overlap stat distribution"""
    args = get_parser()
    jobdir = args.jobdir
    nproc = args.nproc
    debug = args.debug
    jobdir = os.path.abspath(jobdir)
    if debug:
        setup_log(log, file_name='log.out', level=logging.DEBUG)
    else:
        setup_log(log, file_name='log.out', level=logging.INFO)

    outdir = os.path.join(jobdir, 'outfigs')
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    log.info("Validating FALCON_ROOT %s...", jobdir)

    raw, pread, overlaps = validate_falcon_root(jobdir)

    if raw:
        raw_lengths = plot_length_distribution_raw(jobdir)
    if pread:
        pread_lengths = plot_length_distribution_preads(jobdir)
    if raw and pread:
        plot_dual_lengths(raw_lengths, pread_lengths)
    if overlaps:
        plot_ovlp_stats(jobdir, nproc)
    if not raw and not pread and not overlaps:
        log.info("No data found, are you sure %s is a FALCON job_root?", jobdir)
    log.info("Finished! Find your plots here: %s", outdir)

    return


def get_parser():
    """Return an argparse instance"""

    __version__ = 0.1
    parser = argparse.ArgumentParser(version=__version__)
    parser.add_argument("jobdir", type=str, default='./',
                        help='path to a complete FALCON job directory')
    parser.add_argument("--nproc", type=int, default=4)
    parser.add_argument('--debug', action='store_true',
                        help="Print debug logging to stdout")

    return parser.parse_args()


if __name__ == "__main__":
    sys.exit(main())
