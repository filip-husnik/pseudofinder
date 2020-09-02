#!/usr/bin/env python3

from . import reannotate, common

import sys
import os
import re
import shutil
from contextlib import contextmanager

import pandas as pd
import numpy
from plotly.offline import plot
from plotly.graph_objs import Surface, Layout, Scene, Figure


@contextmanager
def suppress_output_to_console():
    """Prevents reannotate.py from writing to stdout, which would be a huge mess."""
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout


def settings_loop(args):
    """This function will run annotate.py with a variety of settings."""

    interval = 1 / args.resolution  # interval will be used to print the progress to stdout

    # Preparing arguments for reannotate.py
    command_line_args = args
    logged_args = common.parse_log(args.logfile)
    args = common.reconcile_args(command_line_args, logged_args)
    basename = args.outprefix

    for length_pseudo in numpy.arange(0.0, 1.01, interval):
        print("%s\tCollecting data: %d%% completed" % (common.current_time(), round(length_pseudo*100, 2)), end='\r')

        for shared_hits in numpy.arange(0.0, 1.01, interval):
            args.length_pseudo = length_pseudo
            args.shared_hits = shared_hits
            args.outprefix = "%s/L%s_S%s" % (basename, length_pseudo, shared_hits)

            with suppress_output_to_console():  # Prevents writing to stdout
                reannotate.reannotate(args)

    args.outprefix = basename  # Have to put this back to its original value
    print('')  # Necessary because the previous print was rolling back on itself


def parse_summary_files(args):
    """This function will parse the summary files for the values needed to generate a 3D plot."""

    outfile = open(args.outprefix + '_matrix.tsv', 'w')
    outfile.write("length_pseudo\tshared_hits\tpseudogenes\n")  # header for the output file

    for root, dirs, files in os.walk(args.outprefix):
        log_files = [file for file in files if "_log.txt" in file]
        for file in log_files:
            path = os.path.join(root, file)
            with open(path, 'r') as infile:
                lines = infile.readlines()
                data = []
                # regex to find 'length_pseudo', 'shared_hits', 'num_pseudos'
                for line in lines:
                    line = line.replace('\n', '')
                    sep = '\t'
                    if re.match('Length_pseudo', line):
                        data.append(float(line.split(sep=sep)[1]))
                    elif re.match('Shared_hits', line):
                        data.append(float(line.split(sep=sep)[1]))
                    elif re.match('Pseudogenes \(total\)', line):
                        data.append(float(line.split(sep=sep)[1]))
            outfile.write('\t'.join(map(str, data))+"\n")
    outfile.close()


def make_plot(args):
    """This function will generate a 3D surface plot."""
    raw_data = pd.read_csv(args.outprefix+'_matrix.tsv', sep="\t", dtype=float,
                           names=['length_pseudo', 'shared_hits', 'vals'], header=0)
    matrix = raw_data.pivot(index='length_pseudo', columns='shared_hits', values='vals')

    data = [Surface(x=matrix.columns,
                    y=matrix.index,
                    z=matrix.values)]

    layout = Layout(
        scene=Scene(
            xaxis=dict(title='shared_hits',
                       autorange=True),
            yaxis=dict(title='length_pseudo',
                       autorange=True),
            zaxis=dict(title=args.title,
                       autorange=True)
        )
    )

    fig = Figure(data=data, layout=layout)
    plot(fig, filename=args.outprefix+".html", auto_open=False)
    common.print_with_time("Figure plotted: %s.html" % args.outprefix)


def main():
    args = common.get_args('visualize')
    args.length_pseudo = None
    args.shared_hits = None
    args.dnds_out = None
    # Reset the folder specified to contain the outputs
    if os.path.exists(args.outprefix):
        shutil.rmtree(args.outprefix)
    os.makedirs(args.outprefix)

    settings_loop(args)  # Generates the data
    parse_summary_files(args)  # Collects and builds the matrix
    make_plot(args)  # Plots the data

    if args.keep_files is False:
        shutil.rmtree(args.outprefix)

if __name__ == '__main__':
    main()
