#!/usr/bin/env python3
from . import annotate, genome_map

import argparse
import sys
import re


def get_args():
    parser = argparse.ArgumentParser(
        usage='\033[1m' + "[pseudofinder.py reannotate -g GENOME -p BLASTP -x BLASTX -hc HITCAP -op OUTPREFIX] or "
                          "[pseudofinder.py reannotate --help] for more options." + '\033[0m')

    always_required = parser.add_argument_group('\033[1m' + 'Required arguments' + '\033[0m')

    always_required.add_argument('-g', '--genome', help='Provide your genome file in genbank format.',
                                 required=True)
    always_required.add_argument('-p', '--blastp', help='Specify an input blastp file.',
                                 required=True)
    always_required.add_argument('-x', '--blastx', help='Specify an input blastx file.',
                                 required=True)
    always_required.add_argument('-log', '--logfile', required=True,
                                 help='Provide the log file from the run that generated the blast files.')
    always_required.add_argument('-op', '--outprefix', help='Specify an output prefix.',
                                 required=True)

    optional = parser.add_argument_group('\033[1m' + 'Adjustable parameters' + '\033[0m')

    optional.add_argument('-l', '--length_pseudo', default=None, type=float,
                          help='Please provide percentage of length for pseudo candidates, '
                               'default is 0.60 (60%%). \nExample: \"-l 0.50\" will consider genes that are '
                               'less than 50%% of the average length of similar genes.')
    optional.add_argument('-s', '--shared_hits', default=None, type=float,
                          help='Percentage of blast hits that must be shared in order to join two nearby regions,'
                               ' default is 0.30 (30%%). \nExample: \"-s 0.50\" will merge nearby regions if '
                               'they shared 50%% of their blast hits.')
    optional.add_argument('-it', '--intergenic_threshold', default=None, type=float,
                          help='Number of BlastX hits needed to annotate an intergenic region as a pseudogene.\n'
                               'Calculated as a percentage of maximum number of allowed hits (--hitcap).\n'
                               'Default is %(default)s.')
    optional.add_argument('-d', '--distance', default=None, type=int,
                          help='Maximum distance between two regions to consider joining them. Default is %(default)s.')

    # parse_known_args will create a tuple of known arguments in the first position and unknown in the second.
    # We only care about the known arguments, so we take [0].
    args = parser.parse_known_args()[0]

    return args


def parse_log(logfile: str):

    with open(logfile, 'r') as log:
        for line in log.readlines():
            if re.match("Distance:", line):
                distance = int(line.split(sep="\t")[1])
            elif re.match("hitcap", line):
                hitcap = int(line.split(sep="\t")[1])
            elif re.match("Intergenic_length", line):
                intergenic_length = int(line.split(sep="\t")[1])
            elif re.match("Intergenic_threshold", line):
                intergenic_threshold = float(line.split(sep="\t")[1])
            elif re.match("Length_pseudo", line):
                length_pseudo = float(line.split(sep="\t")[1])
            elif re.match("Shared_hits", line):
                shared_hits = float(line.split(sep="\t")[1])
            elif re.match("Database", line):
                database = line.split(sep="\t")[1]

    log_dict = {
        'distance': distance,
        'hitcap': hitcap,
        'intergenic_length': intergenic_length,
        'intergenic_threshold': intergenic_threshold,
        'length_pseudo': length_pseudo,
        'shared_hits': shared_hits,
        'database': database
    }

    return log_dict


def fix_args(command_line_args, logged_args):
    """This function is required to merge arguments given from the command line with arguments parsed
    from the provided log file. """

    args = command_line_args
    # If optional parameters are not declared, use the information from the log file to file them in.
    if args.distance is None:
        args.distance = logged_args['distance']

    if args.length_pseudo is None:
        args.length_pseudo = logged_args['length_pseudo']

    if args.shared_hits is None:
        args.shared_hits = logged_args['shared_hits']

    if args.intergenic_threshold is None:
        args.intergenic_threshold = logged_args['intergenic_threshold']

    args.hitcap = logged_args['hitcap']
    args.database = logged_args['database']
    args.intergenic_length = logged_args['intergenic_length']

    return args


def reannotate(args):

    base_outfile_name = args.outprefix + "_"
    file_dict = {
        'proteome_filename': args.blastp.replace('.blastP_output.tsv', ''),
        'intergenic_filename': args.blastx.replace('.blastX_output.tsv', ''),
        'blastp_filename': args.blastp,
        'blastx_filename': args.blastx,
        'pseudos_gff': base_outfile_name + "pseudos.gff",
        'pseudos_fasta': base_outfile_name + "pseudos.fasta",
        'functional_gff': base_outfile_name + "functional.gff",
        'functional_faa': base_outfile_name + "functional.faa",
        'chromosome_map': base_outfile_name + "map.pdf",
        'log': base_outfile_name + "log.txt"
    }

    # Collect everything from the blast files
    orfs = annotate.parse_blast(fasta_file=file_dict['proteome_filename'], blast_file=file_dict['blastp_filename'], blast_format='blastp')
    intergenic_regions = annotate.parse_blast(fasta_file=file_dict['intergenic_filename'], blast_file=file_dict['blastx_filename'], blast_format='blastx')
    all_regions = orfs + intergenic_regions

    # Sorted list of contigs containing only orfs, no intergenic regions
    orfs_by_contig = annotate.sort_contigs(loc=annotate.split_regions_into_contigs(lori=orfs))
    # Sorted list of contigs containing orfs and intergenic regions
    all_regions_by_contig = annotate.sort_contigs(loc=annotate.split_regions_into_contigs(lori=all_regions))

    pseudogenes = []
    functional_genes = []

    for contig_index, contig in enumerate(all_regions_by_contig):
        print('\033[1m'+'%s\tChecking contig %s / %s for pseudogenes.\033[0m' % (annotate.current_time(),
                                                                                 contig_index+1,
                                                                                 len(all_regions_by_contig))),
        sys.stdout.flush()

        pseudos_on_contig = annotate.annotate_pseudos(args=args, contig=contig)  # Returns 'Contig' data type
        pseudogenes.extend(pseudos_on_contig.regions)  # List of regions

        try:
            functional_genes_on_contig = annotate.get_functional_genes(contig=orfs_by_contig[contig_index],
                                                                       pseudos=pseudos_on_contig.regions)
            functional_genes.extend(functional_genes_on_contig.regions)
        except IndexError:  # If there are no orfs on a small contig, an error will be thrown when checking that contig.
            continue

        print('\t\t\tNumber of ORFs on this contig: %s\n'
              '\t\t\tNumber of pseudogenes flagged: %s' % (
                  len([region for region in contig.regions if region.region_type == annotate.RegionType.ORF]),
                  len(pseudos_on_contig.regions))),
        sys.stdout.flush()

    # Write all output files
    annotate.write_genes_to_gff(args, lopg=pseudogenes, gff=file_dict['pseudos_gff'])
    annotate.write_genes_to_gff(args, lopg=functional_genes, gff=file_dict['functional_gff'])
    annotate.write_pseudos_to_fasta(args, pseudofinder_regions=pseudogenes, outfile=file_dict['pseudos_fasta'])
    # TODO: Activate this feature once you finish writing it
    # write_functional_to_fasta(infile=file_dict['proteome_filename'], outfile=file_dict['functional_faa'],
    #                           contigs=functional_genes)
    genome_map.full(genome=args.genome, gff=file_dict['pseudos_gff'], outfile=file_dict['chromosome_map'])
    annotate.write_summary_file(args=args, file_dict=file_dict)
    annotate.reset_statistics_dict()


def main():
    # Declare variables used throughout the rest of the program
    command_line_args = get_args()
    logged_args = parse_log(command_line_args.logfile)
    args = fix_args(command_line_args, logged_args)

    # Do the reannotation
    reannotate(args)


if __name__ == '__main__':
    main()
