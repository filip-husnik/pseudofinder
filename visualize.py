#!/usr/bin/env python3

from plotly.offline import plot
from plotly.graph_objs import Surface, Layout, Scene, Figure
import pandas as pd
import numpy
import sys, os
import argparse
import re
import shutil

def get_args():
    parser = argparse.ArgumentParser(
        usage='\033[1m'+"[pseudo_finder.py visualize -g GENOME -op OUTPREFIX -p BLASTP -x BLASTX -hc HITCAP] or "
                        "[pseudo_finder.py visualize --help] for more options."+'\033[0m')

    ##########################  Always required #################################
    always_required = parser.add_argument_group('\033[1m' + 'Required arguments' + '\033[0m')

    always_required.add_argument('-g', '--genome',
                          help='Please provide your genome file in the genbank format.',
                          required=True)
    always_required.add_argument('-op', '--outprefix',
                        help='Specify an output prefix. A folder will also be created with this name.',
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
    optional = parser.add_argument_group('\033[1m' + 'Optional parameters' + '\033[0m')

    optional.add_argument('-r', '--resolution',
                          default=20,
                          type=int,
                          help='Specifies the resolution of your 3D plot. '
                               'Lowering the resolution will increase the speed of this process.'
                               '\'-r 100\' means the program will generate a 100x100 data matrix. Default is %(default)s.')
    optional.add_argument('-k', '--keep_files',
                          default=False,
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

    if args.hitcap is None:
        sys.stderr.write("Error: Hitcap is required to correctly call pseudogenes. The value can be found in the log file from your annotation.\n")
        exit()

    return args


def settings_loop(args):
    "This function will run annotate.py with a variety of settings."
    interval = 1 / args.resolution
    path_to_annotate = os.path.dirname(__file__)+"/annotate.py"

    for length_pseudo in numpy.arange(0.0,1.01,interval):
        print("Collecting data: %d%% completed" % round(length_pseudo*100, 2), end='\r')
        for shared_hits in numpy.arange(0.0,1.01,interval):
            filename = "%s/L%s_S%s" % (args.outprefix, length_pseudo, shared_hits)
            # Runs annotate.py as if you were using the command line
            os.system("python %s "
                   "--genome %s "
                   "--outprefix %s "
                   "--blastp %s "
                   "--blastx %s "
                   "--length_pseudo %s "
                   "--shared_hits %s "
                   "--hitcap %s "
                   "--outformat 0"
                   " > /dev/null" % (path_to_annotate,
                                     args.genome,
                                     filename,
                                     args.blastp,
                                     args.blastx,
                                     length_pseudo,
                                     shared_hits,
                                     args.hitcap))
    print('\n') #Necessary because the previous print was rolling back on itself


def parse_summary_files(args):
    "This function will parse the summary files for the values needed to generate a 3D plot."

    outfile = open(args.outprefix + '_matrix.tsv', 'w')
    outfile.write("length_pseudo\tshared_hits\tpseudogenes\n") #header for the output file

    for root, dirs, files in os.walk(args.outprefix):
        for file in files:
            path = os.path.join(root,file)
            with open(path, 'r') as infile:
                lines = infile.readlines()
                data = []
                #regex to find 'length_pseudo', 'shared_hits', 'num_pseudos'
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
    "This function will generate a 3D surface plot."
    raw_data = pd.read_csv(args.outprefix+'_matrix.tsv', sep="\t", dtype=float, names=['length_pseudo', 'shared_hits', 'vals'], header=0)
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
    plot(fig, image_filename=args.outprefix, image='png')


def main():
    args = get_args()

    #Reset the folder specified to contain the outputs
    if os.path.exists(args.outprefix):
        shutil.rmtree(args.outprefix)
    os.makedirs(args.outprefix)

    settings_loop(args) #Generates the data
    parse_summary_files(args) #Collects and builds the matrix
    make_plot(args) #Plots the data

    if args.keep_files is False:
        shutil.rmtree(args.outprefix)