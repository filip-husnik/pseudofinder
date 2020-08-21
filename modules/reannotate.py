#!/usr/bin/env python3
from . import common, annotate, genome_map

import sys
import re
import ast


def parse_log(logfile: str):
    log_dict = {}
    with open(logfile, 'r') as log:
        for line in log.readlines():
            line = line.replace('\n', '')
            sep = '\t'
            for arg in common.get_args(names_only=True):
                if re.match(arg, line, re.IGNORECASE):
                    item = str(line.split(sep=sep)[1])
                    try:
                        log_dict[arg] = ast.literal_eval(item)
                    except SyntaxError:
                        log_dict[arg] = item

    log_dict['log_outprefix'] = log_dict['blastp'].replace('_proteome.faa.blastP_output.tsv', '')
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
    args.log_outprefix = logged_args['log_outprefix']
    args.reference = logged_args['reference']

    return args


def reannotate(args):

    file_dict = common.file_dict(args)
    log_file_dict = common.file_dict(args, outprefix=args.log_outprefix)

    # Collect everything from the blast files
    orfs = annotate.parse_blast(fasta_file=log_file_dict['proteome_filename'], blast_file=log_file_dict['blastp_filename'], blast_format='blastp')
    intergenic_regions = annotate.parse_blast(fasta_file=log_file_dict['intergenic_filename'], blast_file=log_file_dict['blastx_filename'], blast_format='blastx')
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
    command_line_args = common.get_args('reannotate')
    logged_args = parse_log(command_line_args.logfile)
    args = fix_args(command_line_args, logged_args)

    # Do the reannotation
    reannotate(args)


if __name__ == '__main__':
    main()
