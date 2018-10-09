#!/usr/bin/env python3

from . import reannotate

import sys
import os
import argparse
import re
import shutil
import subprocess
from time import localtime, strftime
from contextlib import contextmanager

import pandas as pd
import numpy
from plotly.offline import plot
from plotly.graph_objs import Surface, Layout, Scene, Figure


def current_time() -> str:
    """Returns the current time when this function was executed."""
    return str(strftime("%Y-%m-%d %H:%M:%S", localtime()))


def get_args():
    parser = argparse.ArgumentParser(
        usage='\033[1m'+"[pseudofinder.py visualize -g GENOME -op OUTPREFIX -p BLASTP -x BLASTX -log LOGFILE] or "
                        "[pseudofinder.py visualize --help] for more options."+'\033[0m')

    # Always required
    always_required = parser.add_argument_group('\033[1m' + 'Required arguments' + '\033[0m')

    always_required.add_argument('-g', '--genome',
                                 help='Please provide your genome file in the genbank format.',
                                 required=True)
    always_required.add_argument('-op', '--outprefix',
                                 help='Specify an output prefix. A folder will also be created with this name.',
                                 required=True)
    always_required.add_argument('-p', '--blastp', required=True,
                                 help='Specify an input blastp file.')
    always_required.add_argument('-x', '--blastx', required=True,
                                 help='Specify an input blastx file.')
    always_required.add_argument('-log', '--logfile', required=True,
                                 help='Provide the log file from the run that generated the blast files.')

    # Optional
    optional = parser.add_argument_group('\033[1m' + 'Optional parameters' + '\033[0m')

    optional.add_argument('-r', '--resolution',
                          default=20,
                          type=int,
                          help='Specifies the resolution of your 3D plot. '
                               'Lowering the resolution will increase the speed of this process.'
                               '\'-r 100\' means the program will generate a 100x100 data matrix. '
                               'Default is %(default)s.')
    optional.add_argument('-k', '--keep_files',
                          default=False,
                          type=bool,
                          help='Specifies whether to keep all output files. \'-k False\' will remove'
                               ' all files except for the graph and data matrix after running. Default is %(default)s.')
    optional.add_argument('-t', '--title',
                          default=None,
                          type=str,
                          help='Specifies a title for your plot. Default is None.')
    optional.add_argument('-it', '--intergenic_threshold', default=None, type=float,
                          help='Number of BlastX hits needed to annotate an intergenic region as a pseudogene.\n'
                               'Calculated as a percentage of maximum number of allowed hits (--hitcap).\n'
                               'Default is %(default)s.')
    optional.add_argument('-d', '--distance', default=None, type=int,
                          help='Maximum distance between two regions to consider joining them. Default is %(default)s.')

    # "parse_known_args" will create a tuple of known arguments in the first position and unknown in the second.
    # We only care about the known arguments, so we take [0].
    args = parser.parse_known_args()[0]

    return args


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
    logged_args = reannotate.parse_log(args.logfile)
    reannotate.fix_args(command_line_args, logged_args)
    basename = args.outprefix

    for length_pseudo in numpy.arange(0.0, 1.01, interval):
        print("%s\tCollecting data: %d%% completed" % (current_time(), round(length_pseudo*100, 2)), end='\r')

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
                    if re.match('Length_pseudo', line):
                        data.append(float(line.split(sep='\t')[1]))
                    elif re.match('Shared_hits', line):
                        data.append(float(line.split(sep='\t')[1]))
                    elif re.match('Pseudogenes \(total\)', line):
                        data.append(float(line.split(sep='\t')[1]))
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
    print("%s\tFigure plotted: %s.html" % (current_time(), args.outprefix))


def main():
    args = get_args()
    args.length_pseudo = None
    args.shared_hits = None

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
