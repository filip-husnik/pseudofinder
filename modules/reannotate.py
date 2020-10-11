#!/usr/bin/env python3
from . import common, annotate, genome_map, dnds


def reannotate(args):
    file_dict = common.file_dict(args)
    log_file_dict = common.file_dict(args, outprefix=args.log_outprefix)

    # Collect everything from the blast files
    orfs = annotate.parse_blast(fasta_file=log_file_dict['proteome_filename'],
                                blast_file=log_file_dict['blastp_filename'],
                                blast_format='blastp')

    intergenic_regions = annotate.parse_blast(fasta_file=log_file_dict['intergenic_filename'],
                                              blast_file=log_file_dict['blastx_filename'],
                                              blast_format='blastx')

    all_regions = orfs + intergenic_regions

    # Sorted list of contigs containing only orfs, no intergenic regions
    orfs_by_contig = annotate.sort_contigs(loc=annotate.split_regions_into_contigs(lori=orfs))
    # Sorted list of contigs containing orfs and intergenic regions
    all_regions_by_contig = annotate.sort_contigs(loc=annotate.split_regions_into_contigs(lori=all_regions))

    pseudogenes = []
    intact_genes = []

    for contig_index, contig in enumerate(all_regions_by_contig):
        common.print_with_time('Checking contig %s / %s for pseudogenes.' % (contig_index + 1, len(all_regions_by_contig)))
        pseudos_on_contig = annotate.annotate_pseudos(args=args, contig=contig)  # Returns 'Contig' data type
        pseudogenes.extend(pseudos_on_contig.regions)  # List of regions

        try:
            intact_genes_on_contig = annotate.get_intact_genes(contig=orfs_by_contig[contig_index],
                                                              pseudos=pseudos_on_contig.regions)
            intact_genes.extend(intact_genes_on_contig.regions)
        except IndexError:  # If there are no orfs on a small contig, an error will be thrown when checking that contig.
            continue

    if args.dnds_out:
        dnds.full(skip=True,
                  ref=args.reference,
                  nucOrfs=file_dict['cds_filename'],
                  pepORFs=file_dict['proteome_filename'],
                  referenceNucOrfs=file_dict['ref_cds_filename'],
                  referencePepOrfs=file_dict['ref_proteome_filename'],
                  c=file_dict['ctl'],
                  dndsLimit=args.max_dnds,
                  M=args.max_ds,
                  m=args.min_ds,
                  threads=1,
                  search="blast",
                  out=file_dict['dnds_out'])

    else:
        for contig_index, contig in enumerate(all_regions_by_contig):
            common.print_with_time('\tNumber of ORFs on this contig: %s\n'
                                   '\t\t\t\tNumber of pseudogenes flagged: %s'
                                   % (len([region for region in contig.regions if region.region_type == annotate.RegionType.ORF]),
                                      len(pseudos_on_contig.regions)))

            annotate.write_summary_file(args=args, file_dict=file_dict)

    # Write all output files
    annotate.write_genes_to_gff(args, lopg=pseudogenes, gff=file_dict['pseudos_gff'])
    annotate.write_genes_to_gff(args, lopg=intact_genes, gff=file_dict['intact_gff'])
    annotate.write_pseudos_to_fasta(args, pseudofinder_regions=pseudogenes, outfile=file_dict['pseudos_fasta'])
    # TODO: Activate this feature once you finish writing it
    # write_intact_to_fasta(infile=file_dict['proteome_filename'], outfile=file_dict['intact_faa'],
    #                           contigs=intact_genes)
    genome_map.full(genome=args.genome, gff=file_dict['pseudos_gff'], outfile=file_dict['chromosome_map'])

    if args.dnds_out:
        annotate.integrate_dnds(func_gff=file_dict['intact_gff'], pseudo_gff=file_dict['pseudos_gff'],
                       dnds_out=file_dict['dnds_out'] + "/dnds-summary.csv", func_faa=file_dict['intact_faa'],
                       func_ffn=file_dict['intact_faa'], cds=file_dict['cds_filename'],
                       proteome=file_dict['proteome_filename'], pseudo_fasta=file_dict['pseudos_fasta'],
                                max_ds=args.max_ds, min_ds=args.min_ds, max_dnds=args.max_dnds)

    annotate.write_summary_file(args=args, file_dict=file_dict)
    annotate.reset_statistics_dict()


def main():
    command_line_args = common.get_args('reannotate')
    logged_args = common.parse_log(command_line_args.logfile)
    args = common.reconcile_args(command_line_args, logged_args)
    reannotate(args)


if __name__ == '__main__':
    main()
