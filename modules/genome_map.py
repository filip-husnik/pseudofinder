#!/usr/bin/env python3

import argparse
import re
from time import localtime, strftime

from reportlab.lib import colors
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature


def current_time() -> str:
    """Returns the current time when this function was executed."""
    return str(strftime("%Y-%m-%d %H:%M:%S", localtime()))


def get_args():
    parser = argparse.ArgumentParser(
        usage='\033[1m'+"[pseudofinder.py map -g GENOME -gff GFF -op OUTPREFIX] or "
                        "[pseudofinder.py map --help] for more options."+'\033[0m')

    # Always required
    always_required = parser.add_argument_group('\033[1m' + 'Required arguments' + '\033[0m')

    always_required.add_argument('-g', '--genome',
                                 help='Provide your genome file in genbank format.',
                                 required=True)
    always_required.add_argument('-gff', '--gff',
                                 help='Provide your pseudogene calls in GFF format.',
                                 required=True)
    always_required.add_argument('-op', '--outprefix',
                                 help='Specify an output prefix.',
                                 required=True)

    # TODO: If at some point it would be good to have the plots have a title, this is ready.
    # optional = parser.add_argument_group('\033[1m' + 'Optional parameters' + '\033[0m')
    #
    # optional.add_argument('-t', '--title',
    #                       default='map',
    #                       type=str,
    #                       help='Specifies a title for your plot. Default is %(default)s.')

    # "parse_known_args" will create a tuple of known arguments in the first position and unknown in the second.
    # We only care about the known arguments, so we take [0].
    args = parser.parse_known_args()[0]

    return args


def read_gbk(args) -> SeqIO.SeqRecord:
    """Reads the input genome file and concatenates all contigs into a single SeqRecord."""
    whole_record = SeqIO.SeqRecord(seq="", id="", name="", features=None)  # Blank SeqRecord that will be added to
    for record in SeqIO.parse(handle=args.genome, format='genbank'):  # Merge all contigs into one large SeqRecord
        whole_record += record
    return whole_record


def read_gff(args) -> SeqIO.SeqRecord:
    """Reads the input GFF file and creates a single SeqRecord with all features."""

    feature_list = []   # Blank list that will be added to the SeqRecord at the end
    contig_dict = {}    # Stores info specific for each contig in the GFF file

    with open(args.gff, 'r') as gff:
        lines = gff.readlines()
        previous_absolute_end = 0   # GFF is written with every contig starting at 1, so we use this to track distance
        for line in lines:
            if re.match("##sequence-region", line):  # These lines just list the lengths of all the contigs
                name, length = line.split(sep=" ")[1], int(line.split(sep=" ")[3])
                contig_dict[name] = {'name': name,
                                     'length': length,
                                     'absolute_start': previous_absolute_end+1}
                # Add the length so that the next contig will start at that value:
                previous_absolute_end += contig_dict[name]['length']

            elif not re.match("#", line):   # These lines will contain Feature information
                contig = line.split(sep="\t")[0]
                start = int(line.split(sep="\t")[3])
                end = int(line.split(sep="\t")[4])
                strand = line.split(sep="\t")[6]
                # Convert strand from symbol to number
                if strand == '+':
                    strand = 1
                elif strand == "-":
                    strand = -1
                feature_list.append(SeqFeature(FeatureLocation(contig_dict[contig]['absolute_start']+start,
                                                               contig_dict[contig]['absolute_start']+end,
                                                               strand=strand), type='gene'))
        # Store all features in this record
        record = SeqIO.SeqRecord(seq="", id="", name="", features=feature_list)

    return record


# Loosely based on tutorial at http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc254
def make_diagram(args, genome_record: SeqIO.SeqRecord, pseudo_record: SeqIO.SeqRecord):
    """Plots the genome with pseudogenes on another track"""
    diagram = GenomeDiagram.Diagram()

    original_features = GenomeDiagram.FeatureSet()  # These features will be from the original genbank file
    for feature in genome_record.features:  # genome_record is the record from the original genbank file
        if feature.type != "gene":
            # Exclude this feature
            continue
        if len(original_features) % 2 == 0:    # Alternate colours
            color = colors.blue
        else:
            color = colors.lightblue
        original_features.add_feature(feature, color=color)

    track_for_original_features = GenomeDiagram.Track(name="Original Features",
                                                      scale_largetick_interval=100000,
                                                      scale_largeticks=5,
                                                      scale_fontangle=180,
                                                      scale_fontsize=10)
    track_for_original_features.add_set(original_features)
    diagram.add_track(track=track_for_original_features, track_level=1)

    pseudo_features = GenomeDiagram.FeatureSet()    # These features will be from the pseudogene annotation
    for feature in pseudo_record.features:
        if len(pseudo_features) % 2 == 0:    # Alternate colours
            color = colors.red
        else:
            color = colors.lightcoral
        pseudo_features.add_feature(feature, color=color)

    track_for_pseudogenes = GenomeDiagram.Track(name="Pseudogenes", scale_largetick_labels=0)
    track_for_pseudogenes.add_set(pseudo_features)
    diagram.add_track(track=track_for_pseudogenes, track_level=2)

    diagram.draw(format="circular", circular=True,
                 start=0, end=len(genome_record), circle_core=0.8)
    diagram.write(filename=args.outprefix+".pdf", output="PDF")
    print("%s\tFigure plotted: %s.pdf" % (current_time(), args.outprefix))


def main():
    args = get_args()
    base_record = read_gbk(args)
    pseudos_record = read_gff(args)
    make_diagram(args, base_record, pseudos_record)

if __name__ == '__main__':
    main()
