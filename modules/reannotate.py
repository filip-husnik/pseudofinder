#!/usr/bin/env python3
from . import common, annotate, sleuth
import shutil

StatisticsDict = annotate.StatisticsDict


def prepare_data_for_analysis(args, file_dict, log_file_dict):
    """
    Converts genbank file to one seqrecord entry per contig, and makes the following changes:
    - creates seqfeatures for each intergenic region
    - adds extra qualifiers to each feature (see annotate.add_qualifiers_to_features)
    - adds dnds / blast results to each feature
    - updates intergenic regions to reflect where blast information shows potential pseudogenes to be located
    """
    genome = annotate.gbk_to_seqrecord_list(args, args.genome)
    proteome = annotate.extract_features_from_genome(args, genome, 'CDS')

    StatisticsDict['NumberOfContigs'] = len(genome)
    StatisticsDict['ProteomeOrfs'] = len(proteome)

    if args.reference:
        if file_dict is not None:  # happens when this function is called from reannotate module and forces sleuth to be redone
            common.print_with_time("Starting Sleuth...")
            ref_genome = annotate.gbk_to_seqrecord_list(args, args.reference)
            ref_proteome = annotate.extract_features_from_genome(args, ref_genome, 'CDS')
            annotate.write_fasta(seqs=ref_proteome, outfile=file_dict['ref_cds_filename'], seq_type='nt')
            annotate.write_fasta(seqs=ref_proteome, outfile=file_dict['ref_proteome_filename'], seq_type='aa')
            sleuth.full(args, file_dict, log_file_dict=log_file_dict)

        else:  # If this function is called from visualize module, sleuth will not be redone
            file_dict = log_file_dict

        sleuth_dict = sleuth.relate_sleuth_data_to_locus_tags(args, file_dict, log_file_dict=log_file_dict)
        annotate.add_sleuth_data_to_genome(args, genome, sleuth_dict)
        shutil.rmtree(file_dict['temp_dir'])  # Delete temp directory now that we are done.

    annotate.add_blasthits_to_genome(args, genome, log_file_dict['blastp_filename'], 'blastp')
    annotate.add_blasthits_to_genome(args, genome, log_file_dict['blastx_filename'], 'blastx')
    annotate.update_intergenic_locations(args, genome)
    return genome


def reannotate(args, genome, file_dict, visualize=False):
    annotate.find_pseudos_on_genome(args, genome)
    common.print_with_time("Collecting run statistics and writing output files.")
    annotate.analysis_statistics(args, genome)
    annotate.write_all_outputs(args, genome, file_dict, 'visualize' if visualize else 'reannotate')


def main():
    command_line_args = common.get_args('reannotate')
    logged_args = common.parse_log_args(command_line_args.logfile)
    args = common.reconcile_args(command_line_args, logged_args)
    file_dict = common.file_dict(args)
    log_file_dict = common.file_dict(args, outprefix=args.log_outprefix)

    genome = prepare_data_for_analysis(args, file_dict, log_file_dict)
    reannotate(args, genome, file_dict)


if __name__ == '__main__':
    main()
