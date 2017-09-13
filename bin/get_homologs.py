#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Simple python script to identify homologous contigs in FALCON primary contig output
Greg Concepcion gconcepcion@pacificbiosciences.com
"""

import os
import csv
import sys
import argparse
import logging
import subprocess
from collections import Counter

import pbcore.io.FastaIO as fi
from falcon_tools import utils

log = logging.getLogger(__name__)
log.setLevel(logging.INFO)


def _get_mummer_bins():
    """I install this into an environment where both MUMmer 3.23 & MUMmer 4.0.0 co-exist
       so this is a little messy. My MUMmer 4.0.0 bins are all suffixed with 4 so this just 
       does a quick check.
    """
    mummer_bins =  ["nucmer", "show-coords", "delta-filter"]
    binlist=[]

    for binary in mummer_bins:
        version4check = "{s}4".format(s=binary)
        try:
            subprocess.call([version4check])
            binlist.append(version4check)
        except OSError as e:
            log.debug("{s} Not found. falling back to {r}".format(s=version4check, r=binary))

            try:
                subprocess.call([binary])
                binlist.append(binary)
            except OSError as e:
                log.error("{s} not found. Please ensure MUMmer 4.0.0 binaries are in your $PATH".format(s=binary))

   
    return binlist[0], binlist[1], binlist[2]

NUCMER_BIN, SHOW_COORDS_BIN, DELTA_FILTER_BIN = _get_mummer_bins()

def run_nucmer(reference, queries, threads):
    """run nucmer for each reference against all queries"""
    refname = os.path.basename(reference.rstrip(".fasta"))

    log.info("Searching all queries for alignments with reference %s", refname)

    if not os.path.exists('deltas'):
        os.mkdir('deltas')
    prefix = "deltas/{f}".format(f=refname)
    cmd = [NUCMER_BIN, "--maxmatch", "-l", "100", "-c", "500",
           "-t", str(threads), "-p", prefix, reference, queries]
    log.debug(cmd)
    stdout, stderr = utils.run(cmd, os.getcwd(), log)

    if stderr:
        log.debug(stderr)

    deltafile = "{p}.delta".format(p=prefix)
    return deltafile


def run_show_coords(deltafile):
    """run show-coords for each reference"""
    log.debug("Converting delta to coords")

    cmd = [SHOW_COORDS_BIN, "-HT", deltafile]
    log.debug(cmd)
    stdout, stderr = utils.run(cmd, os.getcwd(), log)
    coords_file = [line for line in stdout.split(os.linesep)]

    if stderr:
        log.debug(stderr)

    return coords_file


def run_delta_filter(deltafile):
    """Generate filtered delta file for mummerplot"""

    log.debug("filtering delta for plotting")
    new_delta = "{d}_filtered.delta".format(d=deltafile.rstrip('.delta'))
    
    cmd = [DELTA_FILTER_BIN, "-g", deltafile]

    stdout, stderr = utils.run(cmd, os.getcwd(), log)

    if stderr:
        log.debug(stderr)

    with open(new_delta, 'w') as output:
        for line in stdout:
            output.write(line)

    return new_delta


def process_coords(coordsfile, delta, length_dict):
    """Process *.coords output for significant matches"""

    log.debug("Processing coords file")
    reference = os.path.basename(delta.rstrip('.delta'))

    coords = [tuple(i.split()) for i in coordsfile]
    count = Counter(i[8] for i in coords)

    filtered = [i for i, c in count.iteritems() if c > 3]

    query_dict = {}
    self = None
    for query in filtered:
        query_hits = []

        for line in coords:
            if line[7] == line[8]:
                self = line[7]
                #pass
            elif line[8] == query:
                query_hits.append(line)
        query_dict[query] = query_hits
        query_dict.pop(self, None)

    new_list = []

    for key, value in query_dict.iteritems():
        startends = []
        total_bp = sum([int(i[4]) for i in value])
        percent_ref = round(total_bp / float(length_dict[reference]), 4)

        for hit in value:
            startends.append((hit[0], hit[1]))

        if len(merge(startends)) / float(len(startends)) > 0.75:
            if percent_ref >= 0.03:
                ratio = len(merge(startends)) / float(len(startends))
                new_list.append((key, total_bp, percent_ref, ratio))

    log.info("%s shares homology with %s", reference, ",".join([i[0] for i in new_list]))
    qfile = write_qfile(reference, new_list, length_dict)

    return qfile


def write_qfile(reference, contigs, length_dict):
    """Write qfile of homologous IDs for mummerplot"""

    qfiledir = os.path.join(os.getcwd(), 'qfiles')
    if not os.path.exists(qfiledir):
        os.mkdir(qfiledir)
    qfile = os.path.join(qfiledir, "{r}.qfile".format(
        r=os.path.basename(reference.rstrip('.fasta'))))

    with open(qfile, 'w') as fout:

        csv_out = csv.writer(fout, delimiter=' ', lineterminator="\n")

        for contig in contigs:
            name = contig[0].rstrip('|arrow')
            length = length_dict[name]
            csv_out.writerow((contig[0], length, '+'))

    return qfile


def get_length_dict(fastas):
    """Generate length dictionary for all contigs"""

    length_dict = {}

    for fasta in fastas:
        fastadata = fi.FastaReader(fasta)
        for record in fastadata:
            length_dict[record.header] = len(record.sequence)

    return length_dict



def merge(qhits):
    """Merge overlapping hits together to get a contiguous interval"""

    intervals = [(int(i[0]), int(i[1])) for i in qhits]

    if not intervals:
        return []
    data = []
    for interval in intervals:
        data.append((interval[0], 0))
        data.append((interval[1], 1))
    data.sort()
    merged = []
    stack = [data[0]]
    for i in xrange(1, len(data)):
        datum = data[i]
        if datum[1] == 0:
            # this is a lower bound, push this onto the stack
            stack.append(datum)
        elif datum[1] == 1:
            if stack:
                start = stack.pop()
            if len(stack) == 0:
                # we have found our merged interval
                merged.append((start[0], datum[0]))
    return merged


def get_mummerplot_cmd(delta, qfile):
    """Generate mummerplot command for each reference"""
    prefix = os.path.basename(delta).rstrip('.delta')
    new_delta = run_delta_filter(delta)

    cmd = "mummerplot -layout -Q {q} -postscript -p {p} {d}".format(q=qfile,
                                                                    d=new_delta,
                                                                    p=prefix)
    return cmd


def write_bed(reference, query, qhits):
    """Write *.bed annotation file"""

    ref = reference.split('|')[0]
    que = query.split('|')[0]
    outdir = os.path.join(os.getcwd(), 'bedfiles')

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    bedname = "{r}_{q}.bed".format(r=ref, q=que)
    bedfile = os.path.join(outdir, bedname)

    with open(bedfile, 'w') as out:
        csv_out = csv.writer(out, delimiter='\t')

        for row in qhits:
            rstart, rend, _, _, _, _, _ = row
            csv_out.writerow((reference, int(rstart), int(rend), query))

    return bedfile


def get_parser():
    """Return an argparse instance"""

    __version__ = 0.1
    parser = argparse.ArgumentParser(version=__version__)
    parser.add_argument("infile", type=str)
    parser.add_argument("--nproc", type=int, default=8)
    parser.add_argument("--log", type=str, default=None)

    parser.add_argument('--debug', action='store_true',
                        help="Print debug logging to stdout")

    return parser.parse_args()


def main():
    """Main run loop"""

    args = get_parser()
    infile = args.infile
    threads = args.nproc
    debug = args.debug
    logfile = args.log

    if debug:
        utils.setup_log(log, file_name=logfile, level=logging.DEBUG)
    else:
        utils.setup_log(log, file_name=logfile, level=logging.INFO)

    if infile.endswith(('.fasta', '.fa')):
        fastas = utils.explode_fasta(infile, log)
    else:
        log.info("Please provide FASTA as your input file")

    length_dict = get_length_dict(fastas)
    total_seqs = len(length_dict.keys())
    length_sum = sum(length_dict.values())
    log.info("Beginning homology search in %s", infile)
    log.info("Total Contigs: %d", total_seqs)
    log.info("Total Bp: %d", length_sum)

    plot_out = 'plots.sh'

    with open(plot_out, 'w') as plot_out:
        for fasta in sorted(fastas):

            delta = run_nucmer(fasta, infile, threads)
            coordsfile = run_show_coords(delta)

            qfile = process_coords(coordsfile, delta, length_dict)
            plot_cmd = get_mummerplot_cmd(delta, qfile)

            plot_out.write(plot_cmd + '\n')


if __name__ == "__main__":
    sys.exit(main())
