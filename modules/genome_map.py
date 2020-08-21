#!/usr/bin/env python3

from . import common
import re
from reportlab.lib import colors
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature


def read_gbk(genome: str) -> SeqIO.SeqRecord:
    """Reads the input genome file and concatenates all contigs into a single SeqRecord."""
    whole_record = SeqIO.SeqRecord(seq="", id="", name="", features=None)  # Blank SeqRecord that will be added to
    for record in SeqIO.parse(handle=genome, format='genbank'):  # Merge all contigs into one large SeqRecord
        whole_record += record
    return whole_record


def read_gff(gff: str) -> SeqIO.SeqRecord:
    """Reads the input GFF file and creates a single SeqRecord with all features."""

    feature_list = []   # Blank list that will be added to the SeqRecord at the end
    contig_dict = {}    # Stores info specific for each contig in the GFF file

    with open(gff, 'r') as gff:
        lines = gff.readlines()
        previous_absolute_end = 0   # GFF is written with every contig starting at 1, so we use this to track distance
        for line in lines:
            if re.match("##sequence-region", line):  # These lines just list the lengths of all the contigs
                name = line.split(sep=" ")[1]
                length = int(line.split(sep=" ")[3])
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
def make_diagram(genome_record: SeqIO.SeqRecord, pseudo_record: SeqIO.SeqRecord, outfile: str):
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
    diagram.write(filename=outfile, output="PDF")


def main():
    args = common.get_args('genome_map')
    base_record = read_gbk(args.genome)
    pseudos_record = read_gff(args.gff)
    make_diagram(base_record, pseudos_record, args.outprefix)


# genome_map.full() allows this module to be called from another module, which is what happens in annotate.main()
def full(genome: str, gff: str, outfile: str):
    base_record = read_gbk(genome)
    pseudos_record = read_gff(gff)
    make_diagram(base_record, pseudos_record, outfile)
