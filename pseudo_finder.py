#!/usr/bin/env python3

"""
pseudo_finder.py: A script to find pseudogene candidates in annotated genbank files.
Tested mostly on .gbk files annotated by Prokka with the --compliant flag.

Installation requirements:
    python3
    3rd party libraries: biopython
    libraries from the Python standard library: see below

 """

from Bio.Blast.Applications import NcbiblastpCommandline, NcbiblastxCommandline
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from typing import NamedTuple, List
from time import localtime, strftime
import os
import argparse
import re
import csv
import sys

__author__ = "Mitch Syberg-Olsen & Filip Husnik"
__version__ = "0.07"
__maintainer__ = "Filip Husnik"
__email__ = "filip.husnik@gmail.com"


#Data definitions

#An individual blast hit to a region.
BlastHit = NamedTuple('BlastHit', [('accession',str),
                                   ('slen', int),
                                   ('s_start',int),
                                   ('s_end',int),
                                   ('eval',float)])

#All information about a given region (either an ORF or intergenic region).
RegionInfo = NamedTuple('RegionInfo', [('contig', str),
                                       ('query', str),
                                       ('start', int),
                                       ('end', int),
                                       ('strand', str),
                                       ('hits', List[BlastHit]),
                                       ('note', str)])

#A collection of regions (ORFs and intergenic regions) on the same contig.
Contig = NamedTuple('Contig', [('regions', List[RegionInfo]),
                               ('name', str),
                               ('number',int)])

#Global variables, which will be called to write the log file
StatisticsDict = {
                    'BlastpFilename':'',
                    'BlastxFilename':'',
                    'NumberOfContigs':0,
                    'ProteomeOrfs':0,
                    'FragmentedOrfs':0,
                    'PseudogenesTotal':0,
                    'PseudogenesShort':0,
                    'PseudogenesFragmented':0
                  }

def current_time()-> str:
    '''
    Returns the current time. When this function was executed. 
    '''
    return str(strftime("%Y-%m-%d %H:%M:%S", localtime()))


def get_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     epilog='\033[1m'+'Explanation of output choices\n'+'\033[0m'
                                            '\'--outformat 1\': GFF file containing only ORFs flagged as pseudogenes.\n'
                                            '\'--outformat 2\': GFF file containing only ORFs not flagged as pseudogenes. *Not yet implemented*\n'
                                            '\'--outformat 3\': FASTA file containing only ORFs flagged as pseudogenes. *Not yet implemented*\n'
                                            '\'--outformat 4\': FAA file containing only ORFs not flagged as pseudogenes.\n'
                                            'Note: Outputs can be combined. If \'-of 12\' is specified, files \'1\' and \'2\' will be written.')

    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    ##########################  Always required #################################
    ##############################################################################
    always_required = parser.add_argument_group('\033[1m' + 'Required arguments' + '\033[0m')

    always_required.add_argument('-g', '--genome',
                          help='Please provide your genome file in the genbank format.',
                          required=True)
    always_required.add_argument('-op', '--outprefix',
                        help='Specify an output prefix.',
                        required=True)

    ##############################################################################
    ##############################################################################
    required = parser.add_argument_group('\033[1m' + 'Required arguments if BLAST files are not provided' + '\033[0m')

    required.add_argument('-db', '--database',
                        help='Please provide name (if $BLASTB is set on your system) or absolute path of your blast database.')

    ##############################################################################
    ##############################################################################
    sometimes_required = parser.add_argument_group('\033[1m' + 'Required arguments if providing BLAST files' + '\033[0m')

    sometimes_required.add_argument('-p', '--blastp',
                        help='Specify an input blastp file.')

    sometimes_required.add_argument('-x', '--blastx',
                        help='Specify an input blastx file.')

    ##############################################################################
    ##############################################################################
    optional = parser.add_argument_group('\033[1m' + 'Adjustable parameters' + '\033[0m')



    optional.add_argument('-of', '--outformat',
                          default=1,
                          type=int,
                          help='Specifies which style of output to write. Default is %(default)s.\nSee below for explanation.')

    optional.add_argument('-t', '--threads',
                        help='Please provide total number of threads to use for blast, default is 4.',
                        default=4)

    optional.add_argument('-i', '--intergenic_length',
                        help='Please provide length of intergenic regions to check, default is 30 bp.',
                        default=30,
                        type=int)

    optional.add_argument('-l', '--length_pseudo',
                        help='Please provide percentage of length for pseudo candidates, default is 0.60 (60%%). \nExample: \"-l 0.50\" will consider genes that are less than 50%% of the average length of similar genes.',
                        default=0.60,
                        type=float)

    optional.add_argument('-s', '--shared_hits',
                        help='Percentage of blast hits that must be shared in order to join two nearby regions, default is 0.30 (30%%). \nExample: \"-s 0.50\" will merge nearby regions if they shared 50%% of their blast hits.',
                        default=0.30,
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
                               'the program will ignore any sequence after the last ORF on a contig.')

    # TODO: implement in code
    optional.add_argument('-it', '--intergenic_threshold',
                          default=0.30,
                          help='Number of BlastX hits needed to annotate an intergenic region as a pseudogene.\n'
                               'Calculated as a percentage of maximum number of allowed hits (--hitcap).\n'
                               'Default is %(default)s.')

    args = parser.parse_args()

    if args.blastx is None and args.blastp is None and args.database is None:
        parser.error("Pseudofinder requires a database input unless blast results are provided.")

    if (args.blastx is not None and args.blastp is None) or (args.blastp is not None and args.blastx is None):
        parser.error("Pseudofinder requires both blastP and blastX inputs.")

    return args


def get_proteome(args, faa: str) -> None:
    '''
    Parse genbank input file for coding sequences (CDSs) and write them to the output file with coordinates.
    '''

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
    '''
    Parse genbank input file for intergenic regions and write them to the output file with coordinates.

    Copied/modified from "get_interregions" by Iddo Friedberg & Ian MC Fleming
    Released under Biopython license. http://www.biopython.org/DIST/LICENSE
    The original code extracts all regions strand-dependently, even if there is a gene on the other strand
    Such strand information is not needed here, so I arbitrarily select plus strand sequence
    '''
    # Resets 'fasta' if it contains content already
    open(fasta, 'w').close()

    # Parse all contigs in the multicontig genbank
    for contig in SeqIO.parse(args.genome, "genbank"):  # contig = all information for an entire contig
        gene_list = []  # List of coding regions extracted from genbank file.
        intergenic_records = []  # List of intergenic regions that has been extracted from in between coding regions.
        
        for feature in contig.features:  # Loop over the contig, get the gene features on each of the strands
            if feature.type == 'CDS':  # Only present if prokka was run with --compliant flag
                start_position = feature.location._start.position
                end_position = feature.location._end.position
                gene_list.append((start_position, end_position))

        if args.contig_ends is True:
            #Put 'gene' at the start of the contig (position 0). This will force the next 'for loop' to consider intergenic space
            #between position '0' and the beginning of the first gene.
            gene_list.insert(0, (0, 0))

            contig_end = len(contig.seq) #Apparently this is the fastest way to retrieve the end of a contig
            #Append a 'gene' the end of the contig. This will force the next 'for loop' to consider intergenic space between
            #the last gene and the end of the contig.
            gene_list.append((contig_end, contig_end))

        for i, gene in enumerate(gene_list):
            # Compare current start position to previous end position
            last_end = gene_list[i - 1][1]
            this_start = gene_list[i][0]

            if this_start - last_end >= args.intergenic_length:  # Default 30bp.

                IntergenicRegion = SeqRecord(seq=contig.seq[last_end:this_start],       # Nucleotide sequence in range
                                             id="%s_ign_%d" % (contig.name, i),         # Individual ID
                                             description="%s %d-%d %s" % (contig.name,  # Description including name,
                                                                          last_end + 1, # start position
                                                                          this_start,   # end position
                                                                          "+"))         # strand (default +)

                intergenic_records.append(IntergenicRegion)

        # Write to the intergenic records file
        SeqIO.write(intergenic_records, open(fasta, "a"), "fasta")

    print('%s\tIntergenic regions extracted from:\t%s\n'
          '\t\t\tWritten to file:\t\t\t%s.' % (current_time(), args.genome, fasta,)),
    sys.stdout.flush()


def run_blastp(args, faa: str) -> None:
    '''
    run BLASTP with FAA file against DB of your choice.
    '''

    print('%s\tBlastP executed with %s threads.' % (current_time(), args.threads)),
    sys.stdout.flush()

    blastp_cline = NcbiblastpCommandline(query=faa,
                                         num_threads=args.threads,
                                         db=args.database,
                                         num_alignments=args.hitcap,
                                         evalue=args.evalue,
                                         outfmt= "\'7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send slen evalue bitscore sscinames\'",
                                         out=faa + ".blastP_output.tsv")
    blastp_cline()


def run_blastx(args, fasta: str) -> None:
    '''
    run BLASTX with FNA file against DB of your choice.
    '''

    print('%s\tBlastX executed with %s threads.' % (current_time(), args.threads)),
    sys.stdout.flush()

    blastx_cline = NcbiblastxCommandline(query=fasta,
                                         num_threads=args.threads,
                                         db=args.database,
                                         num_alignments=args.hitcap,
                                         evalue=args.evalue,
                                         outfmt="\'7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send slen evalue bitscore sscinames\'",
                                         out=fasta + ".blastX_output.tsv")
    blastx_cline()


#try: parse locus and size from everything
#TODO: find the last locus tag number and continue from there
    #locus_tag
def make_gff_header(args, gff: str) -> None:
    '''
    Prepare a GFF file header (a list of sequence identifiers)
    '''

    #open file specified by input 'gff'
    with open(gff, "w") as gff_output_handle:

        #first line
        gff_output_handle.write("##gff-version 3\n"
                                "#!annotation-date\t%s\n" % (current_time()))

        #writes one line for each contig
        for i, seq_record in enumerate(SeqIO.parse(args.genome, "genbank")):

            gff_output_handle.write(
                "%s %s %s %s\n" % ("##sequence-region",
                                      "gnl|Prokka|%s" % seq_record.id,   #seqid
                                      1,                      #start
                                      (len(seq_record))))     #end


def collect_query_ids(filename: str) -> List[str]:
    '''
    Reads a TSV file and returns a list of Query names, to be used later.
    '''
    loq = []  # list of Query names

    with open(filename) as tsvfile:
        lines = tsvfile.readlines()

        for line in lines:
            if re.match("^# Query:", line):
                FieldsInLine = list(filter(None, re.split("\s|\[|\]|:|\(|\)", line)))
                query = FieldsInLine[2]
                loq.append(query)

    return loq


def parse_blast(filename: str, blast_format: str) -> List[RegionInfo]:
    '''
    This function needs to take a blast query and extract the relevant information (RegionInfo).
    '''
    print('%s\tExtracting information from %s file.' % (current_time(),blast_format)),
    sys.stdout.flush()

    #Dictionary of information relating to each query
    QueryDict = {}
    #the final list of regions
    regionList = []

    with open(filename, 'r') as tsvfile:
        lines = tsvfile.readlines()

        #This will soon be replaced by an actual query, just have to get past the first line
        query = "Placeholder query that wont matching anything because it's way too long"

        for line_number, line in enumerate(lines):
            #matching line example: "# Query: COGCCIIJ_00001 COGCCIIJ_1 [115:223](+)"
            if re.match("^# Query:", line):
                #FieldsInLine splits all fields and filters unintentional whitespace
                #example: "['#', 'Query', 'COGCCIIJ_00001', 'COGCCIIJ_1', '115', '223', '+']"
                FieldsInLine = list(filter(None, re.split("\s|(?<=[0-9])-|\[|\]|:|\(|\)", line)))

                #This will be the name of the query that is currently being looked at
                query = FieldsInLine[2]

                # collect contig, start, end, strand from fields, add to dictionary
                QueryDict[query] = {'contig':FieldsInLine[3],
                                    'query':FieldsInLine[2],
                                    'start':int(FieldsInLine[4]),
                                    'end':int(FieldsInLine[5]),
                                    'strand':FieldsInLine[6],
                                    'hits':[]}

                # If you're parsing a BlastP file, keep track of how many ORFs are in the file
                if blast_format == "BlastP":
                    StatisticsDict['ProteomeOrfs'] += 1
                else:
                    pass

            # Matches the current query at the front of the line
            # matching line example: "COGCCIIJ_00002	sp|P86052|CYC4_THIRO	47.929	169	81	5	61	225	25	190	192	1.33e-40	140"
            elif re.match("^%s" % query, line):

                #FieldsInLine acts the same as above
                #example: "['COGCCIIJ_00002', 'sp|P86052|CYC4_THIRO', '47.929', '169', '81', '5', '61', '225', '25', '190', '192', '1.33e-40', '140']"
                FieldsInLine = list(filter(None, re.split("\s|\[|\]|:|\(|\)", line)))

                #This chunk of code is needed to prevent getting an error from trying to append to a dictionary key that does not exist.
                #Check if the list exists
                try:
                    QueryDict[query]['hits']
                #If it does not, make it an empty list
                except KeyError:
                    QueryDict[query]['hits'] = []

                #Append hit info to list
                QueryDict[query]['hits'].append(BlastHit(accession=FieldsInLine[1],
                                                         slen=int(FieldsInLine[10])*3,
                                                         s_start=int(FieldsInLine[6]),
                                                         s_end=int(FieldsInLine[7]),
                                                         eval=float(FieldsInLine[11])))

    #Once all lines have been checked, write the results to a final list in the form of RegionInfo
    for key in QueryDict:
        if blast_format == "BlastP":
            regionList.append(RegionInfo(contig=QueryDict[key]['contig'],
                                         query=QueryDict[key]['query'],
                                         start=(QueryDict[key]['start']),
                                         end=(QueryDict[key]['end']),
                                         strand=QueryDict[key]['strand'],
                                         hits=QueryDict[key]['hits'],
                                         note='From BlastP;colour=51 153 102'))

        #Have to modify range for intergenic regions
        if blast_format == "BlastX":
            #retrieve actual intergenic range based on blast hits
            try:
                regionStart,regionEnd = get_intergenic_query_range(QueryDict[key]['hits'],QueryDict[key]['start'])
            #If there are no blast hits, this region will not be considered
            except ValueError:
                regionStart,regionEnd = (0,0)

            regionList.append(RegionInfo(contig=QueryDict[key]['contig'],
                                         query=QueryDict[key]['query'],
                                         start=regionStart,
                                         end=regionEnd,
                                         strand=QueryDict[key]['strand'],
                                         hits=QueryDict[key]['hits'],
                                         note='From BlastX'))

    return regionList


def get_intergenic_query_range(lobh: List[BlastHit], start_position: int) -> tuple:
    '''
    Calculates the range of an intergenic region, based on the location of blast hits within the whole intergenic region.
    This is necessary because the start and end positions of hits are defined locally in the blast output - this function
    converts them to absolute positions on the contig.
    '''
    # Collect all start and end positions in a list of blast hits
    AllValues = [bh.s_start for bh in lobh] + [bh.s_end for bh in lobh]

    regionStart = start_position + min(AllValues)
    regionEnd = start_position + max(AllValues)

    return (regionStart, regionEnd)


def split_regions_into_contigs(lori: List[RegionInfo]) -> List[Contig]:
    '''
    Takes a list of regions and splits them based on which contig it belongs to.
    Contig is defined above as 'List[RegionInfo]', so 'List[Contig]' is a list of lists.
    '''
    #collects all contig names. Doesn't store duplicates, so a contig name will not be stored more than once.
    contig_names = list(set([ri.contig for ri in lori]))

    #this will store the output
    contig_list = []

    for contig_name in contig_names:
        #splits a contig name (ie. EOKKIDHA_1) by the underscore, and takes the number.
        contig_number = int(contig_name.split(sep='_')[1])
        #this will store the List[RegionInfo] to be contained on a contig
        regions_on_contig = []

        for ri in lori:
            #if the region's contig name matches, it is added to that contig
            if ri.contig == contig_name:
                regions_on_contig.append(ri)

        #once all regions have been added, that list of regions is appended as a 'Contig' to the list of contigs.
        contig_list.append(Contig(regions=regions_on_contig, name=contig_name, number=contig_number))

    return contig_list


def annotate_pseudos(args, contig: Contig) -> List[RegionInfo]:
    '''
    This function will take input blast files and return a list of all pseudogene candidates.
    '''

    #1: Look through list of regions and find individual ORFs that could be pseudogenes.
    IndividualPseudos = check_individual_ORFs(args=args, lori=contig.regions)

    #2: Update list of regions with any pseudogenes that were flagged from step #1.
    UpdatedList = replace_pseudos_in_list(pseudos=IndividualPseudos, genes=contig.regions)

    #3: Check adjacent regions to see if they could be pseudogene fragments.
    #   This function returns two lists: [0] = Individual pseudogenes
    #                                    [1] = Merged pseudogenes
    AllPseudos = check_adjacent_regions(args=args,
                                        lori=UpdatedList,
                                        distance_cutoff=args.distance,
                                        hit_cutoff=args.shared_hits)

    #returns both individual and merged pseudogenes as a single list, with locus tags added.
    return add_locus_tags(lori=(AllPseudos[0] + AllPseudos[1]), contig=contig.name)


def check_individual_ORFs(args, lori: List[RegionInfo]) -> List[RegionInfo]:
    '''
    This function will take an input of regions and return a list of individual ORFs that could be pseudogenes.
    '''

    InitialBlastpList = [] #TODO: probably delete the list below
    InitialList = [] #This list will contain all ORFs that have enough blast hits to be considered.
                     #If less than 3 blast hits, its hard to calculate a reliable "AverageDatabaseLength"

    lopg = [] #This list will contain the resulting pseudogenes

    for region in lori:
        #Only include regions that were already as genes from whichever annotation software, and that have at least 2 blast hits.
        if 'BlastP' in region.note and len(region.hits) > 2:
            InitialBlastpList.append(region)
        #TODO: fix this
        # elif 'BlastX' in region.note and len(region.hits)/hitcap > intergenic_cutoff:
        #     pseudo = convert_region_to_pseudo(region=region, pseudo_type='intergenic', ratio=0)
        #     lopg.append(pseudo)


    for region in InitialBlastpList:

        #Retrieves lengths of genes that this region has blasted against
        ListOfDatabaseLengths = [hit.slen for hit in region.hits]

        #Calculates the average length of genes that this region has blasted against
        AverageDatabaseLength = sum(ListOfDatabaseLengths) / len(ListOfDatabaseLengths)

        #Calculates the length of this region
        RegionLength = region.end - region.start

        #Ratio of the region's length to the average length of hits.
        Ratio = (RegionLength/AverageDatabaseLength)

        if Ratio < args.length_pseudo:
            pseudo = convert_region_to_pseudo(region=region,
                                              pseudo_type='ORF',
                                              ratio=Ratio*100)  #Multiplied by 100 to convert to percentage
            lopg.append(pseudo)

    return lopg


def convert_region_to_pseudo(region: RegionInfo, pseudo_type: str, ratio: float) -> RegionInfo:
    '''
    Flags a region as a pseudogene by adding a note, that will appear in the GFF file.
    '''

    #TODO: finish implementing this
    if pseudo_type == 'ORF':
        message = 'Note=pseudogene candidate. ' \
                  'Reason: ORF is %s%% of the average length of hits to this gene.;' \
                  'colour=229 204 255' % (round(ratio, 1))
    elif pseudo_type == 'intergenic':
        message = 'Note=pseudogene candidate. ' \
                  'Reason: Intergenic region with'

    pseudogene = RegionInfo(region.contig,
                            region.query,
                            region.start,
                            region.end,
                            region.strand,
                            region.hits,
                            'Note=pseudogene candidate. Reason: ORF is %s%% of the average length of hits to this gene.;colour=229 204 255' % (round(ratio,1)))
                            #'colour=' makes this region appear coloured in Artemis.
    return pseudogene


def replace_pseudos_in_list(pseudos: List[RegionInfo], genes: List[RegionInfo]) -> List[RegionInfo]:
    '''
    This function prevents duplicates of regions that would occur if a gene was labelled and pseudogene and the original
    gene was not removed from the list.
    '''

    FinalList = []

    #if a pseudogene is present at the same position, write the pseudo
    #if it is not, write the gene
    for gene in genes:

        if pseudo_present(gene, pseudos)[0]:
            FinalList.append(pseudo_present(gene, pseudos)[1])
        else:
            FinalList.append(gene)

    return FinalList


def pseudo_present(gene: RegionInfo, pseudos: List[RegionInfo]) -> tuple:
    '''
    Takes a particular gene and checks if that gene has been flagged as a pseudogene.
    Returns two pieces of information.
    0. If a pseudogene has been annotated at this location
    1. The identity of the (pseudo)gene at this location
    '''

    for pseudo in pseudos:
        if pseudo.start == gene.start:
            return (True, pseudo)
        else:
            pass

    return False, gene


def check_adjacent_regions(args, lori: List[RegionInfo], distance_cutoff: int, hit_cutoff: float) -> tuple:
    '''
    This function will take input blast files and return a list of all pseudogene candidates.

    lori: List of regions you want to run through.
    contig_number: the position of the contig in a list of contigs. Used for printing information.
    cutoff: refer to arg.shared_hits. Percentage of hits shared between two regions to consider joining them.
    '''

    sorted_lori = sorted(lori, key=lambda r: r.start)

    MergedList = []  # List of merged pseudogenes stored as RegionInfo
    IndividualList = [] #List of individual pseudogenes stored as RegionInfo

    i = 0   # Iterator

    while i < len(sorted_lori)-1:

        NewPseudoMade = False

        #compare_regions() checks that the two regions pass certain criteria
        if compare_regions(args, r1=sorted_lori[i], r2=sorted_lori[i + 1]) is True:

            #this boolean will be important later on in this function
            NewPseudoMade = True

            #if they pass, create a pseudogene
            pseudo = join_regions(sorted_lori[i], sorted_lori[i + 1])

            #this is to keep track of overall statistics. If the regions have not yet been annotated by this program,
            #the counter will increase by 1 for each of them.
            for region in [sorted_lori[i], sorted_lori[i + 1]]:
                if 'Predicted fragmentation of a single gene' not in region.note and 'From BlastX' not in region.note:
                    StatisticsDict['FragmentedOrfs'] += 1

            #remove items that were joined together
            del sorted_lori[i + 1]
            del sorted_lori[i]

        #If regions [i] and [i+1] fail to join (above), look at regions [i] and [i+2].
        elif i < len(sorted_lori) - 2:
            if compare_regions(args=args, r1=sorted_lori[i], r2=sorted_lori[i + 2]) is True:

                # this boolean will be important later on in this function
                NewPseudoMade = True

                #if they pass, create a pseudogene
                pseudo = join_regions(sorted_lori[i], sorted_lori[i + 2])

                # this is to keep track of overall statistics. If the regions are not intergenic regions and
                #  have not yet been annotated by this program, the counter will increase by 1 for each of them.
                for region in [sorted_lori[i], sorted_lori[i + 1], sorted_lori[i + 2]]:
                    if 'Predicted fragmentation of a single gene' not in region.note and 'From BlastX' not in region.note:
                        StatisticsDict['FragmentedOrfs'] += 1

                #remove items that were joined together (and [i+1] because it's in between them)
                del sorted_lori[i + 2]
                del sorted_lori[i + 1]
                del sorted_lori[i]

            # If the pieces were not assembled, but one of them is an 'individual pseudogene', it is added to the IndividualList
            elif 'Note=pseudogene candidate. Reason: ORF' in sorted_lori[i].note:
                pseudo = sorted_lori[i]
                IndividualList[:] = [item for item in IndividualList if item.start is not pseudo.start]
                IndividualList.append(pseudo)

        #TODO: check if this cane be changed to an 'if' statement and then remove the statement above?
        #This piece of code will only be accessed if 'i' is almost at the end of the list, otherwise it will be captured immediately above &
        elif 'Note=pseudogene candidate. Reason: ORF' in sorted_lori[i].note:
            pseudo = sorted_lori[i]
            IndividualList[:] = [item for item in IndividualList if item.start is not pseudo.start]
            IndividualList.append(pseudo)

        # If the region in question fits none of the critera, move on.
        else:
            pass

        #this boolean resets to false every loop, so it will only be 'True' if two regions have just been merged together
        if NewPseudoMade is True:
            #This code deletes an item in MergedList if that item has the same start position as the pseudogene.
            #It works like:
            #   MergedList[:] = a new version of MergedList, that contains items from MergedList,
            #   unless that item's start point is the same as the pseudo's start point.
            MergedList[:] = [item for item in MergedList if item.start is not pseudo.start]

            #Adds the merged region to a list to keep track of all merged regions
            MergedList.append(pseudo)

            # Adds the merged region to the original list so that it will continue to be considered
            sorted_lori.append(pseudo)

            #Re-sorts the list, because two regions will have been removed and one new one added (see just above).
            sorted_lori = sorted(sorted_lori, key=lambda r: r.start)

            #Resets the iterator so that new region can be tested by join_regions()
            if i > 0:
                i = i - 1

        #If NewPseudoMade is False, then the iterator moves forward in the list to keep checking new regions.
        else:
            i = i + 1

    #Once the loop finishes, add all statistics to StatisticsDict for reporting in the log file.
    StatisticsDict['PseudogenesTotal'] += len(IndividualList) + len(MergedList)
    StatisticsDict['PseudogenesShort'] += len(IndividualList)
    StatisticsDict['PseudogenesFragmented'] += len(MergedList)

    return (IndividualList, MergedList)


def compare_regions(args, r1: RegionInfo, r2: RegionInfo) -> bool:
    '''
    Takes two regions and decides if they are similar enough to join together.
    '''
    #This if statement is a list of conditions that must be met in order for two regions to be joined
    if (
        region_proximity(r1, r2) < args.distance and      #Closer than cutoff default (1000bp)
        matching_hit_critera(args, r1, r2) is True and    #Have enough matching blast hits
        r1.strand == r2.strand and                        #Same strand
        not ("ign" in r1.query and "ign" in r2.query)     #They are not both intergenic regions
    ):
        return True

    else:
        return False


def region_proximity(r1: RegionInfo, r2: RegionInfo) -> int:
    '''
    Takes two regions and returns their distance from each other in # of nucleotides.
    '''
    #sorts the two regions by starting point, so the math will always be consistent.
    sorted_by_start = sorted([r1,r2], key=lambda r: r.start)

    # substracts the end position of the first from the start position of the second
    # this value can actually be negative if a gene starts before the previous one finishes
    return sorted_by_start[1].start - sorted_by_start[0].end


def matching_hit_critera(args, r1: RegionInfo, r2: RegionInfo) -> bool:
    '''
    This function determines if two regions meet the minimum blast hit criteria to be joined together.
    '''

    if len(r1.hits) is not 0 and len(r2.hits) is not 0:
        #sorts the two regions based on number of blast hits.
        s = sorted([r1, r2], key=lambda r: len(r.hits))

        # math is : (Number of shared hits) / (Total number of hits from the region with the least hits) >= cutoff value.
        if number_of_matching_hits(r1,r2)/len(s[0].hits) >= args.shared_hits:
            return True
        else:
            return False

    else:
        return False


def number_of_matching_hits(r1: RegionInfo, r2: RegionInfo) -> int:
    '''
    This function returns the number of blast hits that two regions have in common.
    '''

    r1_accessions = set([blasthit.accession for blasthit in r1.hits])
    r2_accessions = set([blasthit.accession for blasthit in r2.hits])

    return len(set(r1_accessions) & set(r2_accessions))


def join_regions(r1: RegionInfo, r2: RegionInfo) -> RegionInfo:
    '''
    This function needs to take two regions and merge their locations.
    '''

    # concatenates hits from both regions, discards any duplicates, and sorts them by e-value.
    merged_hits = sort_hits_by_eval(list(set(r1.hits + r2.hits)))

    merged_region = RegionInfo(contig=r1.contig,
                               query=r1.query+","+r2.query+",",
                               start=min([r1.start, r2.start]),
                               end=max([r1.end, r2.end]),
                               strand=r1.strand,
                               hits=merged_hits,
                               note='Note=pseudogene candidate. Reason: Predicted fragmentation of a single gene.;colour=229 204 255') #'colour=' makes this region appear coloured in Artemis.
    return merged_region


def sort_hits_by_eval(lobh: List[BlastHit]) -> List[BlastHit]:
    '''
    Sorts a list of blasthits by e-value from low to high (returning the hit with the lowest evalue first).
    '''

    sorted_list = sorted(lobh, key=lambda r: r.eval)

    return sorted_list


def add_locus_tags(lori: List[RegionInfo], contig: str) -> List[RegionInfo]:
    '''
    Adds numerically increasing locus tags to a list of regions.
    '''

    sorted_by_start = sorted(lori, key=lambda r: r.start)

    FinalList = []

    for counter, region in enumerate(sorted_by_start):
        TaggedRegion = RegionInfo(region.contig,
                                  region.query,
                                  region.start,
                                  region.end,
                                  region.strand,
                                  region.hits,
                                  #adds a locus tag with 4 digits.
                                  # ie, if counter = 2 and contig = 'contig1', result will be 'locus_tag=pseudo_contig_1_0002'
                                  region.note + str(';locus_tag=%s_%04d' % (contig, counter+1)))

        FinalList.append(TaggedRegion)

    return FinalList


def sort_contigs(loc: List[Contig]) -> List[Contig]:
    '''
    Takes a list of contigs and sorts it numerically.
    '''

    sortedlist = sorted(loc, key=lambda c: c.number)

    return sortedlist


def write_genes_to_gff(lopg: List[RegionInfo], gff: str) -> None:
    '''
    Takes an input list of genes and writes them to a GFF file in proper format.
    '''

    with open(gff, 'a') as gff_output_handle:
        for pseudo in lopg:
            gff_output_handle.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %
                                    ("gnl|Prokka|%s" % pseudo.contig,                 #1: seqname - name of the chromosome or scaffold
                                     "pseudo_finder",               #2: source - name of the program that generated this feature
                                     "gene",                        #3: feature - feature type name, e.g. Gene, Variation, Similarity
                                     pseudo.start,                  #4: start - Start position of the feature, with sequence numbering starting at 1.
                                     pseudo.end,                    #5: end - End position of the feature, with sequence numbering starting at 1.
                                     '.',                           #6: score - A floating point value.
                                     pseudo.strand,                 #7: strand - defined as + (forward) or - (reverse).
                                     '.',                           #8: frame - One of '0', '1' or '2'.
                                                                        # '0' indicates that the first base of the feature is the first base of a codon,
                                                                        # '1' that the second base is the first base of a codon, and so on..
                                     pseudo.note))                  #9: attribute - A semicolon-separated list of tag-value pairs,
                                                                        #  providing additional information about each feature.

def get_functional_genes(contig: Contig, pseudos: List[RegionInfo]) -> List[RegionInfo]:
    """"Inspects a contig for genes that have not been annotated as pseudogenes, and returns them."""

    # All regions on a contig, sorted by start position
    RegionList = sorted(contig.regions, key=lambda r: r.start)
    # Begin with all regions on the contig
    FunctionalList = sorted(contig.regions, key=lambda r: r.start)

    #Iterate through regions on the contig
    for region in RegionList:
        for pseudo in pseudos:
            #This will be true if the region is nested within a pseudogene
            if region.start >= pseudo.start and region.end <= pseudo.end:

                #Remove that region from the final list
                FunctionalList.remove(region)
                break  #this will speed up the function by ending the 'for pseudo' loop if a match is successful
            else:
                pass

    return FunctionalList

def get_sequences_from_fasta(infile: str, outfile: str, regions: List[RegionInfo], contig:str) -> None:
    """Parses a multifasta file for regions and returns them in a list."""

    with open(infile, 'r') as input, open(outfile, 'a') as output:
        lines = input.readlines()

        regionIndex = 0

        for line_number, line in enumerate(lines):
            try:
                if re.match("^>%s %s" % (regions[regionIndex].query, contig), line):
                    output.write("%s" % line)
                    output.write("%s" % lines[line_number+1])
                    regionIndex += 1
            except IndexError:
                pass


def write_summary_file(args) -> None:
    '''
    Writes a summary file of statistics from the pseudo_finder run.
    '''

    print('%s\tWriting summary of run.' % (current_time())),
    sys.stdout.flush()

    name = args.outprefix + '_log.txt'

    with open(name, 'w') as logfile:
        logfile.write(
            "####### Summary from pseudo_finder.py #######\n\n"
            "Date/time:\t%s\n\n"
            #TODO: add runtime
            
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
            "Shared_hits:\t%s\n\n"
            
            "####### Statistics #######\n"
            "#Input:\n"
            "Initial ORFs:\t%s\n"
            "Number of contigs:\t%s\n"
            "#Output:\n"
            "Inital ORFs joined:\t%s\n"
            "Pseudogenes (total):\t%s\n"
            "Pseudogenes (too short):\t%s\n"
            "Pseudogenes (fragmented):\t%s\n"
            "Pseudogenes (no predicted ORF):\t\n" #TODO: add dict entry to track this
            "Functional genes:\t%s\n\n"

            "####### Output Key #######\n"
            "Initial ORFs joined:\t\tThe number of input open reading frames that have been merged and flagged as a fragmented pseudogene.\n"
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
                      StatisticsDict['ProteomeOrfs'],
                      StatisticsDict['NumberOfContigs'],
                      StatisticsDict['FragmentedOrfs'],
                      StatisticsDict['PseudogenesTotal'],
                      StatisticsDict['PseudogenesShort'],
                      StatisticsDict['PseudogenesFragmented'],
                      StatisticsDict['ProteomeOrfs'] - StatisticsDict['FragmentedOrfs'] - StatisticsDict['PseudogenesShort']))

#TODO: notes- Counting the number of genes going into get_functional adds up to the total number of ORFs found in the proteome file
#TODO: notes- Counting the number of functional genes from the get_functional adds up to the number of sequences found in the FAA file


def main():
    #Collect arguments from parser
    args = get_args()

    #If blast files are not provided, must run blast.
    if args.blastp is None and args.blastx is None:

        #files generated:
        ProteomeFilename = args.outprefix + "_" + os.path.basename(args.genome) + "_proteome.faa"
        IntergenicFilename = args.outprefix + "_" + os.path.basename(args.genome) + "_intergenic.fasta"
        BlastpFilename = ProteomeFilename + ".blastP_output.tsv"
        BlastxFilename = IntergenicFilename + ".blastX_output.tsv"

        #Collect sequences
        get_proteome(args=args, faa=ProteomeFilename)
        get_intergenic_regions(args=args, fasta=IntergenicFilename)

        #Run blast
        run_blastp(args=args, faa=ProteomeFilename)
        run_blastx(args=args, fasta=IntergenicFilename)

    #If blast files are provided, use them.
    else:
        BlastpFilename = args.blastp
        BlastxFilename = args.blastx
        ProteomeFilename = BlastpFilename.replace(".blastP_output.tsv","")
        IntergenicFilename = BlastxFilename.replace(".blastX_output.tsv","")

    #BlastP and BlastX files have just been formally declared, so now we will add their names to the StatisticsDict
    StatisticsDict['BlastpFilename'] = BlastpFilename
    StatisticsDict['BlastxFilename'] = BlastxFilename

    # Collect everything from the blast files
    ORFs = parse_blast(filename=BlastpFilename, blast_format='BlastP')
    IntergenicRegions = parse_blast(filename=BlastxFilename, blast_format='BlastX')
    # Add all the regions together
    AllRegions = ORFs + IntergenicRegions

    #Split into contigs
    ORFsByContig = sort_contigs(loc=split_regions_into_contigs(lori=ORFs))              # Sorted list of contigs containing only ORFs, no intergenic regions
    AllRegionsByContig = sort_contigs(loc=split_regions_into_contigs(lori=AllRegions))  # Sorted list of contigs containing ORFs and intergenic regions

    #Add number of contigs to the StatisticsDict
    StatisticsDict['NumberOfContigs'] = len(AllRegionsByContig)

    ######## PREPARE FOR WRITING RESULTS ########
    OutputTypes = set(str(args.outformat))  #converts args.outformat into a set of integers that can be checked
    StatisticsDict['OutputFiles'] = []

    if '1' in OutputTypes:
        PseudosGff = args.outprefix + "_" + os.path.basename(args.genome) + "_pseudos.gff"
        make_gff_header(args=args, gff=PseudosGff)
        StatisticsDict['OutputFiles'].append("Output [1]: %s" % PseudosGff)

    if '2' in OutputTypes:  #TODO: Not yet properly implemented
        print('\033[1m'+"Notice! --outformat 2 specified but this output has not yet been implemented."+'\033[0m'),
        sys.stdout.flush()
        FunctionalGff = args.outprefix + "_" + os.path.basename(args.genome) + "_functional.gff"
        make_gff_header(args=args, gff=FunctionalGff)
        StatisticsDict['OutputFiles'].append("Output [2]: %s" % FunctionalGff)

    if '3' in OutputTypes: #TODO: not yet properly implemented
        print('\033[1m' + "Notice! --outformat 3 specified but this output has not yet been implemented." +'\033[0m'),
        sys.stdout.flush()
        # faa_filename_3 = args.outprefix + "_" + os.path.basename(args.genome) + "_pseudos.faa"

    if '4' in OutputTypes:
        FunctionalFaa = args.outprefix + "_" + os.path.basename(args.genome) + "_functional.faa"
        #Clear file contents if this file already exists:
        open(FunctionalFaa, 'w').close()
        StatisticsDict['OutputFiles'].append("Output [4]: %s" % FunctionalFaa)

    #For each contig:
    for contig_index, contig in enumerate(AllRegionsByContig):

        # Reports the contig number being examined.
        print('%s\tChecking contig %s / %s for pseudogenes.' % (current_time(),
                                                                contig_index+1,  #indices are 0 based so I added 1 to make it more intuitive
                                                                len(AllRegionsByContig))),
        sys.stdout.flush()

        #Annotate pseudogenes
        PseudoGenes = annotate_pseudos(args=args, contig=contig)

        #Retrieve functional genes by comparing ORFs to annotated pseudogenes
        try:
            FunctionalGenes = get_functional_genes(contig=ORFsByContig[contig_index], pseudos=PseudoGenes)
        except IndexError: #If there are no ORFs on a small contig, an error will be thrown when trying to check that contig.
            break

        #Write the appropriate output types
        if '1' in OutputTypes: # Write GFF file with pseudogenes in it
            write_genes_to_gff(lopg=PseudoGenes, gff=PseudosGff)

        if '2' in OutputTypes: # TODO: implement this if it is ever needed
            write_genes_to_gff(lopg=FunctionalGenes, gff=FunctionalGff)

        # if '3' in OutputTypes: #TODO: implement this if it is ever needed
        #     get_sequences_from_fasta(infile=faa_filename, #from the beginning of main()
        #                              outfile=faa_filename_3, #generated just above
        #                              regions=pseudos)

        if '4' in OutputTypes: #Write FAA file with protein coding genes in it
            get_sequences_from_fasta(infile=ProteomeFilename,
                                     outfile=FunctionalFaa,
                                     regions=FunctionalGenes,
                                     contig=contig.name)

    #Finish by writing a summary file
    write_summary_file(args=args)


if __name__ == '__main__':
    main()

def notes():
#################################################################################################################
# ADDITIONAL FEATURES TO INCORPORATE?
#################################################################################################################
#TODO: in the future :)
    #TODO: sometimes ORFs are predicted by mistake on the opposite strand (e.g. in GC-rich genomes), check regions with ORFS with no blastp hits by blastx
    #TODO: include an analysis of cryptic pseudogenes based on dN/dS ratios (PAML?...)
    #TODO: include an optional analysis when RNA-Seq data are available
    #TODO: exclude mobile elements such as transposases? or treat them differently
    #TODO: visualize results by matplotlib with a simple scatter plot of all genes/pseudogenes (dN/dS, GC content, length ratio, ...)
    #TODO: put all output files into a separate folder, include folder name among arguments to parse
    #TODO: check if files are already in the folder before starting -- rewrite/append/ask

    pass
