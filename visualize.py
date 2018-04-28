#!/usr/bin/env python3

from plotly.offline import plot
from plotly.graph_objs import Surface, Layout, Scene, Figure
import pandas as pd
import numpy
from os import system, path
from sys import argv
import argparse
import re


def get_args():
    parser = argparse.ArgumentParser(usage='\033[1m'+"TODO: write usage."+'\033[0m')

    ##########################  Always required #################################
    always_required = parser.add_argument_group('\033[1m' + 'Required arguments' + '\033[0m')

    always_required.add_argument('-g', '--genome',
                          help='Please provide your genome file in the genbank format.',
                          required=True)
    always_required.add_argument('-op', '--outprefix',
                        help='Specify an output prefix.',
                        required=True)
    always_required.add_argument('-p', '--blastp',
                        help='Specify an input blastp file.')
    always_required.add_argument('-x', '--blastx',
                        help='Specify an input blastx file.')
    always_required.add_argument('-hc', '--hitcap',
                          type=int,
                          help='The hitcap from used to generate your blast files is required.\n\n')

    ##############################################################################
    ##############################################################################
    optional = parser.add_argument_group('\033[1m' + 'Adjustable parameters' + '\033[0m')

    optional.add_argument('-r', '--resolution',
                          default=20,
                          type=int,
                          help='Specifies the resolution of your 3D plot. '
                               'Lowering the resolution will increase the speed of this process.'
                               '\'-r 100\' means the program will generate a 100x100 data matrix. Default is %(default)s.')
    optional.add_argument('-k', '--keep_files',
                          default=True,
                          type=bool,
                          help='Specifies whether to keep all output files. \'-k False\' will remove all files except for'
                               ' the graph and data matrix after running. Default is %(default)s.')
    optional.add_argument('-t', '--title',
                          default=None,
                          type=str,
                          help='Specifies a title for your plot. Default is None.')

    #parse_known_args will create a tuple of known arguments in the first position and unknown in the second.
    #We only care about the known arguments, so we take [0].
    args = parser.parse_known_args()[0]

    return args

def settings_loop(args):
    "This function will run annotate.py with a variety of settings."
    interval = 1 / args.resolution
    path_to_annotate = path.dirname(__file__)+"/annotate.py"

    for length_pseudo in numpy.arange(0.0,1.01,interval):
        for shared_hits in numpy.arange(0.0,1.01,interval):

            filename = "%s/L%s_S%s" % (args.outprefix, length_pseudo, shared_hits)

            #Runs annotate.py as if you were using the command line
            system("python %s "
                   "--genome %s "
                   "--outprefix %s "
                   "--blastp %s "
                   "--blastx %s "
                   "--length_pseudo %s "
                   "--shared_hits %s "
                   "--hitcap %s " % (path_to_annotate,
                                     args.genome,
                                     filename,
                                     args.blastp,
                                     args.blastx,
                                     length_pseudo,
                                     shared_hits,
                                     args.hitcap))


def parse_summary_files(args):
    "This function will parse the summary files for the values needed to generate a 3D plot."
    #pseudocode:
    for file in folder_of_summary_files:
        with open(file, 'r') as infile:
            lines = file.readlines()
            #regex to find 'length_pseudo', 'shared_hits', 'num_pseudos'
            for line in lines:
                if re.match('length_pseudo', line):
                    parse_dict['length_pseudo'] = 'blah'

    with open(args.outprefix+'table.tsv', 'w') as outfile:
        outfile.write("Put all the dictionary vals here")


def make_plot(args):
    "This function will generate a 3D surface plot."
    raw_data = pd.read_csv(args.outprefix+'table.tsv', sep="\t", dtype=float, names=['length_pseudo', 'shared_hits', 'vals'])
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
    plot(fig, image_filename=argv[2], image='png')


def main():
    print('\033[1m' + "This feature is currently under construction." + '\033[0m')
    args = get_args()
    settings_loop(args)
    parse_summary_files(args)
    make_plot(args)
