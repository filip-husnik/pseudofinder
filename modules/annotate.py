#!/usr/bin/env python3

import os
import argparse
import re
import sys
from typing import NamedTuple, List
from enum import Enum
from time import localtime, strftime

from Bio.Blast.Applications import NcbiblastpCommandline, NcbiblastxCommandline
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


# Data definitions
# An individual blast hit to a region.
BlastHit = NamedTuple('BlastHit', [('accession', str),
                                   ('slen', int),
                                   ('s_start', int),
                                   ('s_end', int),
                                   ('eval', float)])

# All possible types of regions
RegionType = Enum('RegionType', ['ORF',
                                 'intergenic',
                                 'shortpseudo',
                                 'fragmentedpseudo',
                                 'intergenicpseudo'])

# All information about a given region (either an ORF or intergenic region).
RegionInfo = NamedTuple('RegionInfo', [('contig', str),
                                       ('query', str),
                                       ('start', int),
                                       ('end', int),
                                       ('strand', str),
                                       ('hits', List[BlastHit]),
                                       ('note', str),
                                       ('region_type', RegionType)])


# A collection of regions (ORFs and intergenic regions) on the same contig.
Contig = NamedTuple('Contig', [('regions', List[RegionInfo]),
                               ('name', str),
                               ('number', int)])

# Global dictionary, which will be called to write the log file
StatisticsDict = {
                    'BlastpFilename': '',
                    'BlastxFilename': '',
                    'NumberOfContigs': 0,
                    'ProteomeOrfs': 0,
                    'FragmentedOrfs': 0,
                    'PseudogenesTotal': 0,
                    'PseudogenesShort': 0,
                    'PseudogenesFragmented': 0,
                    'PseudogenesIntergenic': 0
                  }


# Functions
def current_time() -> str:
    """Returns the current time. When this function was executed."""
    return str(strftime("%Y-%m-%d %H:%M:%S", localtime()))


def get_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     usage='\033[1m'+"[pseudo_finder.py annotate -g GENOME -op OUTPREFIX] or "
                                           "[pseudo_finder.py annotate --help] for more options."+'\033[0m',
                                     epilog='\033[1m'+'Explanation of output choices\n'+'\033[0m'
                                            '\'--outformat 1\': GFF file containing only ORFs flagged as pseudogenes.\n'
                                            '\'--outformat 2\': GFF file containing only ORFs not flagged as pseudogenes. *Not yet implemented*\n'
                                            '\'--outformat 3\': FASTA file containing only ORFs flagged as pseudogenes. *Not yet implemented*\n'
                                            '\'--outformat 4\': FAA file containing only ORFs not flagged as pseudogenes.\n'
                                            'Note: Outputs can be combined. If \'-of 12\' '
                                            'is specified, files \'1\' and \'2\' will be written.')

    # TODO: I think this should be removed as the version will now be stored in pseudo_finder.py
    # parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    # Always required
    always_required = parser.add_argument_group('\033[1m' + 'Required arguments' + '\033[0m')
    always_required.add_argument('-g', '--genome',
                                 help='Please provide your genome file in the genbank format.',
                                 required=True)
    always_required.add_argument('-op', '--outprefix',
                                 help='Specify an output prefix.',
                                 required=True)

    # Sometimes required (if blast files are not provided)
    required = parser.add_argument_group('\033[1m' + 'Required arguments if BLAST files are not provided' + '\033[0m')
    required.add_argument('-db', '--database',
                          help='Please provide name (if $BLASTB is set on your system)'
                               ' or absolute path of your blast database.')

    # Sometimes required (if blast files ARE provided)
    sometimes_required = parser.add_argument_group('\033[1m'+'Required arguments if providing BLAST files'+'\033[0m')

    sometimes_required.add_argument('-p', '--blastp',
                                    help='Specify an input blastp file.')

    sometimes_required.add_argument('-x', '--blastx',
                                    help='Specify an input blastx file.')
    
    # Optional arguments
    optional = parser.add_argument_group('\033[1m' + 'Adjustable parameters' + '\033[0m')

    optional.add_argument('-of', '--outformat',
                          default=1,
                          type=int,
                          help='Specifies which style of output to write. Default is %(default)s.\n'
                               'See below for explanation.')

    optional.add_argument('-t', '--threads',
                          help='Please provide total number of threads to use for blast, default is 4.',
                          default=4)

    optional.add_argument('-i', '--intergenic_length',
                          help='Please provide length of intergenic regions to check, default is 30 bp.',
                          default=30,
                          type=int)

    optional.add_argument('-l', '--length_pseudo',
                          help='Please provide percentage of length for pseudo candidates, '
                               'default is 0.60 (60%%). \nExample: \"-l 0.50\" will consider genes that are '
                               'less than 50%% of the average length of similar genes.',
                          default=0.65,
                          type=float)

    optional.add_argument('-s', '--shared_hits',
                          help='Percentage of blast hits that must be shared in order to join two nearby regions,'
                               ' default is 0.30 (30%%). \nExample: \"-s 0.50\" will merge nearby regions if '
                               'they shared 50%% of their blast hits.',
                          default=0.50,
                          type=float)

    optional.add_argument('-e', '--evalue',
                          help='Please provide e-value for blast searches. Default is 1e-4.',
                          default='1e-4')

    optional.add_argument('-d', '--distance',
                          default=1000,
                          type=int,
                          help='Maximum distance between two regions to consider joining them. Default is %(default)s.')

    optional.add_argument('-hc', '--hitcap',
                          default=15,
                          type=int,
                          help='Maximum number of allowed hits for BLAST. Default is %(default)s.\n\n')

    optional.add_argument('-ce', '--contig_ends',
                          default=False,
                          action='store_true',
                          help='Forces the program to include intergenic regions at contig ends. If not specified,\n'
                               'the program will ignore any sequence before the first ORF '
                               'and after the last ORF on a contig.')

    optional.add_argument('-it', '--intergenic_threshold',
                          default=0.30,
                          help='Number of BlastX hits needed to annotate an intergenic region as a pseudogene.\n'
                               'Calculated as a percentage of maximum number of allowed hits (--hitcap).\n'
                               'Default is %(default)s.')

    # parse_known_args will create a tuple of known arguments in the first position and unknown in the second.
    # We only care about the known arguments, so we take [0].
    args = parser.parse_known_args()[0]

    if args.blastx is None and args.blastp is None and args.database is None:
        parser.error("Pseudofinder requires a database input unless blast results are provided.")

    if (args.blastx is not None and args.blastp is None) or (args.blastp is not None and args.blastx is None):
        parser.error("Pseudofinder requires both blastP and blastX inputs.")

    return args


def get_proteome(args, faa: str) -> None:
    """Parse genbank input file for coding sequences (CDSs) and write them to the output file with coordinates."""

    with open(args.genome, "r") as input_handle:
        with open(faa, "w") as output_handle:
            for seq_record in SeqIO.parse(input_handle, "genbank"):
                for seq_feature in seq_record.features:
                    if seq_feature.type == "CDS":
                        assert len(seq_feature.qualifiers['translation']) == 1
                        output_handle.write(">%s %s %s\n%s\n" % (seq_feature.qualifiers['locus_tag'][0],
                                                                 seq_record.name,
                                                                 seq_feature.location,
                                                                 seq_feature.qualifiers['translation'][0]))

    print('%s\tProteome extracted from:\t\t%s\n'
          '\t\t\tWritten to file:\t\t\t%s.' % (current_time(), args.genome, faa,)),
    sys.stdout.flush()


def get_intergenic_regions(args, fasta: str) -> None:
    """Parse genbank input file for intergenic regions and write them to the output file with coordinates.

    Copied/modified from "get_interregions" by Iddo Friedberg & Ian MC Fleming
    Released under Biopython license. http://www.biopython.org/DIST/LICENSE
    The original code extracts all regions strand-dependently, even if there is a gene on the other strand
    Such strand information is not needed here, so I arbitrarily select plus strand sequence."""
    
    # Resets 'fasta' if it contains content already
    open(fasta, 'w').close()

    # Parse all contigs in the multicontig genbank
    for contig in SeqIO.parse(args.genome, "genbank"):  # contig = all information for an entire contig
        gene_list = []  # List of coding regions extracted from genbank file.
        intergenic_records = []  # List of intergenic regions that has been extracted from in between coding regions.
        
        for feature in contig.features:  # Loop over the contig, get the gene features on each of the strands
            if feature.type == "gene":  # Only present if prokka was run with --compliant flag
                start_position = feature.location._start.position
                end_position = feature.location._end.position
                gene_list.append((start_position, end_position))

        if args.contig_ends is True:
            # Put 'gene' at the start of the contig (position 0). This will force the next 'for loop' 
            # to consider intergenic space between position '0' and the beginning of the first gene.
            gene_list.insert(0, (0, 0))

            contig_end = len(contig.seq)  # Apparently this is the fastest way to retrieve the end of a contig
            # Append a 'gene' the end of the contig. This will force the next 'for loop' to consider 
            # intergenic space between the last gene and the end of the contig.
            gene_list.append((contig_end, contig_end))

        for i, gene in enumerate(gene_list):  # Compare current start position to previous end position
            last_end = gene_list[i - 1][1]
            this_start = gene_list[i][0]

            if this_start - last_end >= args.intergenic_length:  # Default 30bp.

                intergenic_region = SeqRecord(seq=contig.seq[last_end:this_start],        # Nucleotide sequence in range
                                              id="%s_ign_%d" % (contig.name, i),          # Individual ID
                                              description="%s %d-%d %s" % (contig.name,   # Description including name,
                                                                           last_end + 1,  # start position
                                                                           this_start,    # end position
                                                                           "+"))          # strand (default +)

                intergenic_records.append(intergenic_region)

        # Write to the intergenic records file
        SeqIO.write(intergenic_records, open(fasta, "a"), "fasta")

    print('%s\tIntergenic regions extracted from:\t%s\n'
          '\t\t\tWritten to file:\t\t\t%s.' % (current_time(), args.genome, fasta,)),
    sys.stdout.flush()


def run_blastp(args, faa: str) -> None:
    """"Run BLASTP with FAA file against DB of your choice."""

    print('%s\tBlastP executed with %s threads.' % (current_time(), args.threads)),
    sys.stdout.flush()

    blastp_cline = NcbiblastpCommandline(query=faa,
                                         num_threads=args.threads,
                                         db=args.database,
                                         max_target_seqs=args.hitcap,
                                         max_hsps=1,
                                         evalue=args.evalue,
                                         outfmt="\'7 qseqid sseqid pident length mismatch gapopen qstart qend "
                                                "sstart send slen evalue bitscore frames stitle\'",
                                         out=faa + ".blastP_output.tsv")
    blastp_cline()


def run_blastx(args, fasta: str) -> None:
    """Run BLASTX with FNA file against DB of your choice."""

    print('%s\tBlastX executed with %s threads.' % (current_time(), args.threads)),
    sys.stdout.flush()

    blastx_cline = NcbiblastxCommandline(query=fasta,
                                         num_threads=args.threads,
                                         db=args.database,
                                         max_target_seqs=args.hitcap,
                                         max_hsps=1,
                                         evalue=args.evalue,
                                         outfmt="\'7 qseqid sseqid pident length mismatch gapopen qstart qend "
                                                "sstart send slen evalue bitscore frames stitle\'",
                                         out=fasta + ".blastX_output.tsv")
    blastx_cline()


# try: parse locus and size from everything
# TODO: find the last locus tag number and continue from there
#  locus_tag
def make_gff_header(args, gff: str) -> None:
    """Prepare a GFF file header (a list of sequence identifiers)"""
    
    # open file specified by input 'gff'
    with open(gff, "w") as gff_output_handle:
        
        gff_output_handle.write("##gff-version 3\n#!annotation-date\t%s\n" % (current_time()))  # first line
        
        for i, seq_record in enumerate(SeqIO.parse(args.genome, "genbank")):  # writes one line for each contig

            gff_output_handle.write("%s %s %s %s\n" % (
                                                        "##sequence-region",
                                                        "gnl|Prokka|%s" % seq_record.id,   # contig seqid
                                                        1,                                 # contig start
                                                        (len(seq_record)))                 # contig end
                                    )                


def collect_query_ids(filename: str) -> List[str]:
    """Reads a TSV file and returns a list of Query names, to be used later."""

    loq = []  # list of Query names

    with open(filename) as tsvfile:
        lines = tsvfile.readlines()

        for line in lines:
            if re.match("^# Query:", line):
                fields_in_line = list(filter(None, re.split("\s|\[|\]|:|\(|\)", line)))
                query = fields_in_line[2]
                loq.append(query)

    return loq


def parse_blast(filename: str, blast_format: str) -> List[RegionInfo]:
    """This function needs to take a blast query and extract the relevant information (RegionInfo)."""

    print('%s\tExtracting information from %s file.' % (current_time(), blast_format)),
    sys.stdout.flush()

    query_dict = {}  # Dictionary of information relating to each query
    region_list = []  # the final list of regions

    with open(filename, 'r') as tsvfile:
        lines = tsvfile.readlines()

        # This will soon be replaced by an actual query, just have to get past the first line
        query = "Placeholder query that wont matching anything because it's way too long"

        for line_number, line in enumerate(lines):
            # matching line example: "# Query: COGCCIIJ_00001 COGCCIIJ_1 [115:223](+)"
            if re.match("^# Query:", line):
                # fields_in_line splits all fields and filters unintentional whitespace
                # example: "['#', 'Query', 'COGCCIIJ_00001', 'COGCCIIJ_1', '115', '223', '+']"
                fields_in_line = list(filter(None, re.split("\s|(?<=[0-9])-|\[|\]|:|\(|\)", line)))

                query = fields_in_line[2]  # the query that is currently being looked at

                # collect contig, start, end, strand from fields, add to dictionary
                query_dict[query] = {'contig': fields_in_line[3],
                                     'query': fields_in_line[2],
                                     'start': int(fields_in_line[4])+1,
                                     'end': int(fields_in_line[5]),
                                     'strand': fields_in_line[6],
                                     'hits': []}

                # If you're parsing a BlastP file, keep track of how many ORFs are in the file
                if blast_format == "BlastP":
                    StatisticsDict['ProteomeOrfs'] += 1
                else:
                    pass

            # Matches the current query at the front of the line
            # match example: "COGCCIIJ_00002	sp|P86052|CYC4_THIRO	47.929	169	81	5	61	225	25	190	192	1.33e-40	140"
            elif re.match("^%s" % query, line):
                # fields_in_line acts the same as above
                # example: "['COGCCIIJ_00002', 'sp|P86052|CYC4_THIRO', '47.929', '169', '81', '5', '61', '225', '25', '190', '192', '1.33e-40', '140']"
                fields_in_line = list(filter(None, re.split("\s|\[|\]|:|\(|\)", line)))

                # This try/catch is needed to prevent getting an error from trying
                # to append to a dictionary key that does not exist.
                try:
                    query_dict[query]['hits']  # Check if the list exists
                except KeyError:
                    query_dict[query]['hits'] = []  # If it does not, make it an empty list

                # Append hit info to list
                query_dict[query]['hits'].append(BlastHit(accession=fields_in_line[1],
                                                          slen=int(fields_in_line[10])*3,
                                                          s_start=int(fields_in_line[6]),
                                                          s_end=int(fields_in_line[7]),
                                                          eval=float(fields_in_line[11])))

    # Once all lines have been checked, write the results to a final list in the form of RegionInfo
    for key in query_dict:
        if blast_format == "BlastP":
            region_list.append(RegionInfo(contig=query_dict[key]['contig'],
                                          query=query_dict[key]['query'],
                                          start=(query_dict[key]['start']),
                                          end=(query_dict[key]['end']),
                                          strand=query_dict[key]['strand'],
                                          hits=query_dict[key]['hits'],
                                          note='From BlastP;colour=51 153 102',
                                          region_type=RegionType.ORF))

        # Have to modify range for intergenic regions
        if blast_format == "BlastX":
            try:  # retrieve actual intergenic range based on blast hits
                region_start, region_end = get_intergenic_query_range(query_dict[key]['hits'], query_dict[key]['start'])

            except ValueError:  # If there are no blast hits, this region will not be considered
                region_start, region_end = (0, 0)

            region_list.append(RegionInfo(contig=query_dict[key]['contig'],
                                          query=query_dict[key]['query'],
                                          start=region_start,
                                          end=region_end,
                                          strand=query_dict[key]['strand'],
                                          hits=query_dict[key]['hits'],
                                          note='From BlastX',
                                          region_type=RegionType.intergenic))

    return region_list


def get_intergenic_query_range(lobh: List[BlastHit], start_position: int) -> tuple:
    """Calculates the range of an intergenic region, based on the location of blast hits within the whole intergenic region.
    This is necessary because the start and end positions of hits are defined locally in the blast output - this function
    converts them to absolute positions on the contig."""

    # Collect all start and end positions in a list of blast hits
    all_values = [bh.s_start for bh in lobh] + [bh.s_end for bh in lobh]

    region_start = start_position + min(all_values)
    region_end = start_position + max(all_values)

    return region_start, region_end


def split_regions_into_contigs(lori: List[RegionInfo]) -> List[Contig]:
    """Takes a list of regions and splits them based on which contig it belongs to.
    Contig is defined above as 'List[RegionInfo]', so 'List[Contig]' is a list of lists."""

    # collects all contig names. Doesn't store duplicates, so a contig name will not be stored more than once.
    contig_names = list(set([ri.contig for ri in lori]))
    contig_list = []  # this will store the output

    for contig_name in contig_names:
        # Finds all numbers in the contig name (ie. '15' in EOKKIDHA_15) and returns them as a single integer
        contig_number = int("".join(re.findall('\d', str(contig_name))))
        regions_on_contig = []  # stores the List[RegionInfo] to be contained on a contig

        for ri in lori:
            if ri.contig == contig_name:  # if the region's contig name matches, it is added to that contig
                regions_on_contig.append(ri)

        # once all regions have been added, that list of regions is appended as a 'Contig' to the list of contigs.
        contig_list.append(Contig(regions=regions_on_contig, name=contig_name, number=contig_number))

    return contig_list


def annotate_pseudos(args, contig: Contig) -> List[RegionInfo]:
    """
    This function will take input blast files and return a list of all pseudogene candidates.
    """

    # 1: Look through list of regions and find individual ORFs that could be pseudogenes.
    individual_pseudos, intergenic_pseudos = check_individual_ORFs(args=args, lori=contig.regions)

    # 2: Update list of regions with any pseudogenes that were flagged from step #1.
    updated_list = replace_pseudos_in_list(pseudos=individual_pseudos+intergenic_pseudos, regions=contig.regions)

    # 3: Check adjacent regions to see if they could be pseudogene fragments.
    #   This function returns two lists: [0] = Individual pseudogenes
    #                                    [1] = Merged pseudogenes
    all_pseudos = check_adjacent_regions(args=args,
                                         lori=updated_list)

    # returns both individual and merged pseudogenes as a single list, with locus tags added.
    return add_locus_tags(lori=(all_pseudos[0] + all_pseudos[1]), contig=contig.name)


def check_individual_ORFs(args, lori: List[RegionInfo]) -> tuple:
    """This function will take an input of regions and return two lists:
    [0]: a list of individual ORFs that could be pseudogenes.
    [1]: a list of intergenic regions that could be pseudogenes."""

    initial_blastp_list = []  # This list will contain all ORFs that have enough blast hits to be considered.
    blastp_pseudos = []      # This list will contain the resulting pseudogenes
    blastx_pseudos = []

    for region in lori:
        # Only include regions that were already as genes from whichever
        # annotation software, and that have at least 2 blast hits.
        if region.region_type == RegionType.ORF and len(region.hits) > 2:
            initial_blastp_list.append(region)

        # Include a blastx hit if it meets the minimum criteria defined by args.intergenic_threshold
        # For example, if a blastx region has 5 blast hits, the blast hitcap is 15 hits, and the threshold is 0.20,
        # the region will pass. ( 5/15 >= 0.2 ) is true.
        elif region.region_type == RegionType.intergenic and len(region.hits)/args.hitcap >= args.intergenic_threshold:
            pseudo = convert_region_to_pseudo(region=region,
                                              ratio=None,  # this value is only used for BlastP-derived pseudos (below)
                                              number_of_hits=len(region.hits))
            blastx_pseudos.append(pseudo)

    for region in initial_blastp_list:

        # Retrieves lengths of genes that this region has blasted against
        list_of_database_lengths = [hit.slen for hit in region.hits]

        # Calculates the average length of genes that this region has blasted against
        average_database_length = sum(list_of_database_lengths) / len(list_of_database_lengths)

        # Calculates the length of this region
        region_length = region.end - region.start

        # ratio of the region's length to the average length of hits.
        ratio = (region_length/average_database_length)

        if ratio < args.length_pseudo:
            pseudo = convert_region_to_pseudo(region=region,
                                              ratio=ratio*100,      # Multiplied by 100 to convert to percentage
                                              number_of_hits=None)  # Not important for BlastP hits
            blastp_pseudos.append(pseudo)

    return blastp_pseudos, blastx_pseudos


def convert_region_to_pseudo(region: RegionInfo, ratio: float, number_of_hits: int) -> RegionInfo:
    """Flags a region as a pseudogene by adding a note, that will appear in the GFF file.
    Regions must be explicitly rewritten because NamedTuples are immutable."""

    if region.region_type == RegionType.ORF:
        message = 'Note=pseudogene candidate. ' \
                  'Reason: ORF is %s%% of the average length of hits to this gene.;' \
                  'colour=229 204 255' % (round(ratio, 1))  # 'colour=' makes this region appear coloured in Artemis.
        pseudo_type = RegionType.shortpseudo

    elif region.region_type == RegionType.intergenic:
        message = 'Note=pseudogene candidate. ' \
                  'Reason: Intergenic region with %s blast hits.;' \
                  'colour=229 204 255' % number_of_hits  # 'colour=' makes this region appear coloured in Artemis.
        pseudo_type = RegionType.intergenicpseudo

    pseudogene = RegionInfo(contig=region.contig,
                            query=region.query,
                            start=region.start,
                            end=region.end,
                            strand=region.strand,
                            hits=region.hits,
                            note=message,
                            region_type=pseudo_type)

    return pseudogene


def replace_pseudos_in_list(pseudos: List[RegionInfo], regions: List[RegionInfo]) -> List[RegionInfo]:
    """This function prevents duplicates of regions that would occur if a gene was
    labelled and pseudogene and the original gene was not removed from the list."""

    final_list = []

    for region in regions:
        if pseudo_present(region, pseudos)[0]:  # if a pseudogene is present at the same position, write the pseudo
            final_list.append(pseudo_present(region, pseudos)[1])
        else:
            final_list.append(region)  # if it is not, write the gene

    return final_list


def pseudo_present(region: RegionInfo, pseudos: List[RegionInfo]) -> tuple:
    """Takes a particular gene and checks if that gene has been flagged as a pseudogene.
    Returns two pieces of information.
    0. If a pseudogene has been annotated at this location (bool)
    1. The identity of the (pseudo)gene at this location (RegionInfo)"""

    for pseudo in pseudos:
        if pseudo.start == region.start:
            return True, pseudo
        else:
            pass

    return False, region


def check_adjacent_regions(args, lori: List[RegionInfo]) -> tuple:
    """This function will take input blast files and return a list of all pseudogene candidates.

    lori: List of regions you want to run through.
    contig_number: the position of the contig in a list of contigs. Used for printing information.
    cutoff: refer to arg.shared_hits. Percentage of hits shared between two regions to consider joining them."""

    sorted_lori = sorted(lori, key=lambda r: r.start)
    merged_list = []  # List of merged pseudogenes stored as RegionInfo
    individual_list = []  # List of individual pseudogenes stored as RegionInfo
    i = 0   # Iterator

    while i < len(sorted_lori)-1 and len(sorted_lori) > 1:
        new_pseudo_made = False
        try:
            # compare_regions() checks that the two regions pass certain criteria
            if compare_regions(args, r1=sorted_lori[i], r2=sorted_lori[i + 1]) is True:
                new_pseudo_made = True    # this bool will be important later on in this function
                pseudo = join_regions(sorted_lori[i], sorted_lori[i + 1])   # if they pass, create a pseudogene

                # this is to keep track of overall statistics. If the regions are plain ORFs or ORFs annotated
                # as short pseudos, the counter will increase by 1 for each of them.
                for region in [sorted_lori[i], sorted_lori[i + 1]]:
                    if region.region_type == RegionType.ORF or region.region_type == RegionType.shortpseudo:
                        StatisticsDict['FragmentedOrfs'] += 1

                del sorted_lori[i + 1]  # remove items that were joined together
                del sorted_lori[i]

            # If regions [i] and [i+1] fail to join (above), look at regions [i] and [i+2].
            elif compare_regions(args=args, r1=sorted_lori[i], r2=sorted_lori[i + 2]) is True:

                new_pseudo_made = True  # this boolean will be important later on in this function
                pseudo = join_regions(sorted_lori[i], sorted_lori[i + 2])  # if they pass, create a pseudogene

                # same as above ^
                for region in [sorted_lori[i], sorted_lori[i + 1], sorted_lori[i + 2]]:
                    if region.region_type == RegionType.ORF or region.region_type == RegionType.shortpseudo:
                        StatisticsDict['FragmentedOrfs'] += 1

                del sorted_lori[i + 2]  # remove items that were joined together, and [i+1] because it's in between them
                del sorted_lori[i + 1]
                del sorted_lori[i]

            # If the pieces were not assembled but one of them is an 'individual pseudogene',
            # it is added to the individual_list
            elif sorted_lori[i].region_type == RegionType.shortpseudo or sorted_lori[i].region_type == RegionType.intergenicpseudo:
                pseudo = sorted_lori[i]
                # Deletes an item in individual_list if it has the same start position as an individual pseudo.
                individual_list[:] = [item for item in individual_list if item.start is not pseudo.start]
                individual_list.append(pseudo)

            # If the region in question fits none of the critera, move on.
            else:
                pass

        except IndexError:  # This will be triggered when 'i' equals the length of the list of pseudos
            pass

        # 'new_pseudo_made' resets to false every loop
        # so it will only be 'True' if two regions have just been merged together
        if new_pseudo_made is True:
            # Deletes an item in merged_list if that item has the same start position as the pseudogene.
            # It works like:
            #   merged_list[:] = a new version of merged_list, that contains items from merged_list,
            #   unless that item is nested within the new pseudogene.
            merged_list[:] = [item for item in merged_list if (item.start is not pseudo.start) and (item.end is not pseudo.end)]

            # Adds the merged region to a list to keep track of all merged regions
            merged_list.append(pseudo)
            # Adds the merged region to the original list so that it will continue to be considered
            sorted_lori.append(pseudo)

            # Re-sorts the list, because two regions will have been removed and one new one added (see just above).
            sorted_lori = sorted(sorted_lori, key=lambda r: r.start)

            i = i - 1  # Resets the iterator so that new region can be tested by join_regions()

        # If new_pseudo_made is False, then the iterator moves forward in the list to keep checking new regions.
        else:
            i = i + 1

        # This will remove rare cases where a pseudogene isnt handled correctly and remains in the individual_list
        # despite being a part of a merged pseudogene in merged_list. No touchy.
        individual_list[:] = [item for item in individual_list if item.start not in [pseudo.start for pseudo in merged_list]]

    # Once the loop finishes, add all statistics to StatisticsDict for reporting in the log file.
    StatisticsDict['PseudogenesTotal'] += len(individual_list) + len(merged_list)
    StatisticsDict['PseudogenesShort'] += len([item for item in individual_list if item.region_type == RegionType.shortpseudo])
    StatisticsDict['PseudogenesIntergenic'] += len([item for item in individual_list if item.region_type == RegionType.intergenicpseudo])
    StatisticsDict['PseudogenesFragmented'] += len(merged_list)

    return individual_list, merged_list


def compare_regions(args, r1: RegionInfo, r2: RegionInfo) -> bool:
    """Takes two regions and decides if they are similar enough to join together."""

    # A list of conditions that must be met in order for two regions to be joined
    if (
        region_proximity(r1, r2) < args.distance and      # Closer than cutoff default (1000bp)
        matching_hit_critera(args, r1, r2) is True and    # Have enough matching blast hits
        r1.strand == r2.strand and                        # Same strand
        not (r1.region_type == RegionType.intergenic and r2.region_type == RegionType.intergenic)  # They are not both intergenic regions
    ):
        return True

    else:
        return False


def region_proximity(r1: RegionInfo, r2: RegionInfo) -> int:
    """Takes two regions and returns their distance from each other in # of nucleotides."""
    # sorts the two regions by starting point, so the math will always be consistent.
    sorted_by_start = sorted([r1,r2], key=lambda r: r.start)

    # substracts the end position of the first from the start position of the second
    # this value can actually be negative if a gene starts before the previous one finishes
    return sorted_by_start[1].start - sorted_by_start[0].end


def matching_hit_critera(args, r1: RegionInfo, r2: RegionInfo) -> bool:
    """This function determines if two regions meet the minimum blast hit criteria to be joined together."""

    if len(r1.hits) is not 0 and len(r2.hits) is not 0:
        # sorts the two regions based on number of blast hits.
        s = sorted([r1, r2], key=lambda r: len(r.hits))

        # math: (Number of shared hits) / (Total number of hits from the region with the least hits) >= cutoff value.
        if number_of_matching_hits(r1, r2)/len(s[0].hits) >= args.shared_hits:
            return True
        else:
            return False

    else:
        return False


def number_of_matching_hits(r1: RegionInfo, r2: RegionInfo) -> int:
    """This function returns the number of blast hits that two regions have in common."""

    r1_accessions = set([blasthit.accession for blasthit in r1.hits])
    r2_accessions = set([blasthit.accession for blasthit in r2.hits])

    return len(set(r1_accessions) & set(r2_accessions))


def join_regions(r1: RegionInfo, r2: RegionInfo) -> RegionInfo:
    """This function needs to take two regions and merge their locations."""

    # concatenates hits from both regions, discards any duplicates, and sorts them by e-value.
    merged_hits = sort_hits_by_eval(list(set(r1.hits + r2.hits)))

    merged_region = RegionInfo(contig=r1.contig,
                               query=r1.query+","+r2.query+",",
                               start=min([r1.start, r2.start]),
                               end=max([r1.end, r2.end]),
                               strand=r1.strand,
                               hits=merged_hits,
                               note='Note=pseudogene candidate. Reason: Predicted fragmentation of a single gene.;'
                                    'colour=229 204 255',  # 'colour=' makes this region appear coloured in Artemis.
                               region_type=RegionType.fragmentedpseudo)
    return merged_region


def sort_hits_by_eval(lobh: List[BlastHit]) -> List[BlastHit]:
    """Sorts a list of blasthits by e-value from low to high (returning the hit with the lowest evalue first)."""

    sorted_list = sorted(lobh, key=lambda r: r.eval)
    return sorted_list


def sort_contigs(loc: List[Contig]) -> List[Contig]:
    """Takes a list of contigs and sorts it numerically."""

    sortedlist = sorted(loc, key=lambda c: c.number)

    return sortedlist


def add_locus_tags(lori: List[RegionInfo], contig: str) -> List[RegionInfo]:
    """Adds numerically increasing locus tags to a list of regions."""

    sorted_by_start = sorted(lori, key=lambda r: r.start)

    final_list = []

    for counter, region in enumerate(sorted_by_start):
        tagged_region = RegionInfo(region.contig,
                                   region.query,
                                   region.start,
                                   region.end,
                                   region.strand,
                                   region.hits,
                                   # adds a locus tag with 4 digits.
                                   # ie, if counter = 2 and contig = 'contig1', result will be
                                   # 'locus_tag=pseudo_contig_1_0002'
                                   region.note + str(';locus_tag=%s_%04d' % (contig, counter+1)),
                                   region_type=region.region_type)

        final_list.append(tagged_region)

    return final_list


def write_genes_to_gff(lopg: List[RegionInfo], gff: str) -> None:
    """Takes an input list of genes and writes them to a GFF file in proper format."""

    with open(gff, 'a') as gff_output_handle:
        for pseudo in lopg:
            # Write a single line, tab separated string with all the information about a particular feature
            gff_output_handle.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %
                                    ("gnl|Prokka|%s" % pseudo.contig,  # 1
                                     "pseudo_finder",  # 2
                                     "gene",  # 3
                                     pseudo.start,  # 4
                                     pseudo.end,  # 5
                                     '.',  # 6
                                     pseudo.strand,  # 7
                                     '.',  # 8
                                     pseudo.note))  # 9

# Explanation of the string above:
# 1: seqname - name of the chromosome or scaffold
# 2: source - name of the program that generated this feature
# 3: feature - feature type name, e.g. Gene, Variation, Similarity
# 4: start - Start position of the feature, with sequence numbering starting at 1.
# 5: end - End position of the feature, with sequence numbering starting at 1.
# 6: score - A floating point value.
# 7: strand - defined as + (forward) or - (reverse).
# 8: frame - One of '0', '1' or '2'.
        #  '0' indicates that the first base of the feature is the first base of a codon,
        #  '1' that the second base is the first base of a codon, and so on..
# 9: attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.


def get_functional_genes(contig: Contig, pseudos: List[RegionInfo]) -> List[RegionInfo]:
    """"Inspects a contig for genes that have not been annotated as pseudogenes, and returns them."""

    # All regions on a contig, sorted by start position
    region_list = sorted(contig.regions, key=lambda r: r.start)
    # Begin with all regions on the contig
    functional_list = sorted(contig.regions, key=lambda r: r.start)

    # Iterate through regions on the contig
    for region in region_list:
        for pseudo in pseudos:
            # This will be true if the region is nested within a pseudogene
            if region.start >= pseudo.start and region.end <= pseudo.end:
                functional_list.remove(region)  # Remove that region from the final list
                break  # this will speed up the function by ending the 'for pseudo' loop if a match is successful
            else:
                pass

    return functional_list


def get_sequences_from_fasta(infile: str, outfile: str, regions: List[RegionInfo], contig: str) -> None:
    """Parses a multifasta file for regions and returns them in a list."""
    try:
        with open(infile, 'r') as input, open(outfile, 'a') as output:
            lines = input.readlines()

            region_index = 0

            for line_number, line in enumerate(lines):
                try:
                    if re.match("^>%s %s" % (regions[region_index].query, contig), line):
                        output.write("%s" % line)
                        output.write("%s" % lines[line_number+1])
                        region_index += 1
                except IndexError:
                    pass

    except FileNotFoundError:
        sys.stderr.write('\033[1m'+"Error: Pseudo_finder.py expected proteome file "
                                   "to be in the same directory as BlastP file.\n"+'\033[0m')
        exit()


def get_sequences_from_gbk(args, contig_name: str, regions: List[RegionInfo], outfile: str) -> None:
    """Parse genbank input file for the regions provided and write them to the output file with coordinates."""

    # Generates a list of [start:end] tuples that will be used to retrieve sequences
    coord_list = [(region.start, region.end) for region in regions]
    # Where the info for writing to the fasta file will be stored
    fasta_list = []
    # Parse all contigs in the multicontig genbank
    for contig in SeqIO.parse(args.genome, "genbank"):  # contig = all information for an entire contig
        if contig_name == contig.name:
            for counter, coord in enumerate(coord_list):
                fasta_list.append(
                    SeqRecord(seq=contig.seq[coord[0]:coord[1]],
                              id="%s_%04d" % (contig.name, counter+1),
                              description="%s-%s +" % (coord[0], coord[1]))
                )

    # Write to the fasta file
    SeqIO.write(fasta_list, open(outfile, "a"), "fasta")


def write_summary_file(args) -> None:
    """Writes a summary file of statistics from the pseudo_finder run."""

    name = args.outprefix + '_log.txt'

    print('%s\tWriting summary of run:\t%s' % (current_time(), name)),
    sys.stdout.flush()

    with open(name, 'w') as logfile:
        logfile.write(
            "####### Summary from annotate.py #######\n\n"
            "Date/time:\t%s\n\n"
            # TODO: add runtime
            
            "#######    Input Files   #######\n"
            "Genome:\t%s\n"
            "Database:\t%s\n"
            "BlastP:\t%s\n"
            "BlastX:\t%s\n\n"
            
            "#######   Output Files   #######\n"
            "%s\n\n" 
            
            "#######  Settings  #######\n"
            "Intergenic_length:\t%s\n"
            "Length_pseudo:\t%s\n"
            "Shared_hits:\t%s\n"
            "hitcap:\t%s\n\n"
            
            "####### Statistics #######\n"
            "#Input:\n"
            "Initial ORFs:\t%s\n"
            "Number of contigs:\t%s\n"
            "#Output:\n"
            "Inital ORFs joined:\t%s\n"
            "Pseudogenes (total):\t%s\n"
            "Pseudogenes (too short):\t%s\n"
            "Pseudogenes (fragmented):\t%s\n"
            "Pseudogenes (no predicted ORF):\t%s\n"
            "Functional genes:\t%s\n\n"

            "####### Output Key #######\n"
            "Initial ORFs joined:\t\tThe number of input open reading frames "
            "that have been merged and flagged as a fragmented pseudogene.\n"
            "Pseudogenes (too short):\tORFs smaller than the \"shared_hits\" cutoff.\n"
            "Pseudogenes (fragmented):\tPseudogenes composed of merging 2 or more input ORFs.\n"
            "Functional genes:\t\t[Initial ORFs] - [Initial ORFs joined] - [Pseudogenes (too short)]\n"


                   % (current_time(),
                      args.genome,
                      args.database,
                      StatisticsDict['BlastpFilename'],
                      StatisticsDict['BlastxFilename'],
                      "\n".join(StatisticsDict['OutputFiles']),
                      args.intergenic_length,
                      args.length_pseudo,
                      args.shared_hits,
                      args.hitcap,
                      StatisticsDict['ProteomeOrfs'],
                      StatisticsDict['NumberOfContigs'],
                      StatisticsDict['FragmentedOrfs'],
                      StatisticsDict['PseudogenesTotal'],
                      StatisticsDict['PseudogenesShort'],
                      StatisticsDict['PseudogenesFragmented'],
                      StatisticsDict['PseudogenesIntergenic'],
                      StatisticsDict['ProteomeOrfs'] - StatisticsDict['FragmentedOrfs'] - StatisticsDict['PseudogenesShort']))

# TODO: notes- Counting the number of genes going into get_functional
# adds up to the total number of ORFs found in the proteome file
# TODO: notes- Counting the number of functional genes from the get_functional
#  adds up to the number of sequences found in the FAA file


def main():
    # Collect arguments from parser
    args = get_args()

    # If blast files are not provided, must run blast.
    if args.blastp is None and args.blastx is None:

        # files generated:
        proteome_filename = args.outprefix + "_" + os.path.basename(args.genome) + "_proteome.faa"
        intergenic_filename = args.outprefix + "_" + os.path.basename(args.genome) + "_intergenic.fasta"
        blastp_filename = proteome_filename + ".blastP_output.tsv"
        blastx_filename = intergenic_filename + ".blastX_output.tsv"

        # Collect sequences
        get_proteome(args=args, faa=proteome_filename)
        get_intergenic_regions(args=args, fasta=intergenic_filename)

        # Run blast
        run_blastp(args=args, faa=proteome_filename)
        run_blastx(args=args, fasta=intergenic_filename)

    # If blast files are provided, use them.
    else:
        blastp_filename = args.blastp
        blastx_filename = args.blastx
        proteome_filename = blastp_filename.replace(".blastP_output.tsv", "")
        intergenic_filename = blastx_filename.replace(".blastX_output.tsv", "")

    # BlastP and BlastX files have just been formally declared, so now we will add their names to the StatisticsDict
    StatisticsDict['blastp_filename'] = blastp_filename
    StatisticsDict['blastx_filename'] = blastx_filename

    # Collect everything from the blast files
    orfs = parse_blast(filename=blastp_filename, blast_format='BlastP')
    intergenic_regions = parse_blast(filename=blastx_filename, blast_format='BlastX')
    # Add all the regions together
    all_regions = orfs + intergenic_regions

    # Sorted list of contigs containing only orfs, no intergenic regions
    orfs_by_contig = sort_contigs(loc=split_regions_into_contigs(lori=orfs))
    # Sorted list of contigs containing orfs and intergenic regions
    all_regions_by_contig = sort_contigs(loc=split_regions_into_contigs(lori=all_regions))

    # Add number of contigs to the StatisticsDict
    StatisticsDict['NumberOfContigs'] = len(all_regions_by_contig)

    # Prepare for writing results
    output_types = set(str(args.outformat))  # converts args.outformat into a set of integers that can be checked
    StatisticsDict['OutputFiles'] = []

    if '1' in output_types:
        pseudos_gff = args.outprefix + "_" + os.path.basename(args.genome) + "_pseudos.gff"
        make_gff_header(args=args, gff=pseudos_gff)
        StatisticsDict['OutputFiles'].append("Output [1]: %s" % pseudos_gff)

    if '2' in output_types:  # TODO: Not yet properly tested
        print('\033[1m'+"Notice! --outformat 2 specified but this output has not yet been implemented."+'\033[0m'),
        sys.stdout.flush()
        functional_gff = args.outprefix + "_" + os.path.basename(args.genome) + "_functional.gff"
        make_gff_header(args=args, gff=functional_gff)
        StatisticsDict['OutputFiles'].append("Output [2]: %s" % functional_gff)

    if '3' in output_types:  # TODO: not yet properly tested
        print('\033[1m' + "Notice! --outformat 3 specified but this output has not yet been tested." + '\033[0m'),
        sys.stdout.flush()
        fasta_filename_3 = args.outprefix + "_" + os.path.basename(args.genome) + "_pseudos.fasta"
        open(fasta_filename_3, 'w').close()

    if '4' in output_types:
        functional_faa = args.outprefix + "_" + os.path.basename(args.genome) + "_functional.faa"
        # Clear file contents if this file already exists:
        open(functional_faa, 'w').close()
        StatisticsDict['OutputFiles'].append("Output [4]: %s" % functional_faa)

    # For each contig:
    for contig_index, contig in enumerate(all_regions_by_contig):

        # Reports the contig number being examined. Indices are 0 based so I added 1 to make it more intuitive
        print('\033[1m'+'%s\tChecking contig %s / %s for pseudogenes.\033[0m' % (current_time(),
                                                                                 contig_index+1,
                                                                                 len(all_regions_by_contig))),
        sys.stdout.flush()

        # Annotate pseudogenes
        pseudo_genes = annotate_pseudos(args=args, contig=contig)
        print('\t\t\tNumber of orfs on this contig: %s\n'
              '\t\t\tNumber of pseudogenes flagged: %s' % (
                  len([region for region in contig.regions if region.region_type == RegionType.ORF]), len(pseudo_genes))),
        sys.stdout.flush()

        # Retrieve functional genes by comparing orfs to annotated pseudogenes
        try:
            functional_genes = get_functional_genes(contig=orfs_by_contig[contig_index], pseudos=pseudo_genes)
        except IndexError:  # If there are no orfs on a small contig, an error will be thrown when checking that contig.
            continue

        # Write the appropriate output types
        if '1' in output_types:  # Write GFF file with pseudogenes in it
            write_genes_to_gff(lopg=pseudo_genes, gff=pseudos_gff)

        if '2' in output_types:  # TODO: implement this if it is ever needed
            write_genes_to_gff(lopg=functional_genes, gff=functional_gff)

        if '3' in output_types:  # TODO: implement this if it is ever needed
            get_sequences_from_gbk(args,
                                   contig_name=contig.name,
                                   outfile=fasta_filename_3,  #generated just above
                                   regions=pseudo_genes)

        if '4' in output_types:  # Write FAA file with protein coding genes in it
            get_sequences_from_fasta(infile=proteome_filename,
                                     outfile=functional_faa,
                                     regions=functional_genes,
                                     contig=contig.name)

    # Finish by writing a summary file
    write_summary_file(args=args)


if __name__ == '__main__':
    main()


def notes():
    # ADDITIONAL FEATURES TO INCORPORATE?
    # TODO: sometimes ORFs are predicted by mistake on the opposite strand (e.g. in GC-rich genomes), check regions with ORFS with no blastp hits by blastx
    # TODO: include an analysis of cryptic pseudogenes based on dN/dS ratios (PAML?...)
    # TODO: include an optional analysis when RNA-Seq data are available
    # TODO: exclude mobile elements such as transposases? or treat them differently
    # TODO: visualize results by matplotlib with a simple scatter plot of all genes/pseudogenes (dN/dS, GC content, length ratio, ...)
    # TODO: put all output files into a separate folder, include folder name among arguments to parse
    # TODO: check if files are already in the folder before starting -- rewrite/append/ask

    pass
