#!/usr/bin/env python3

"""
pseudo_finder.py: A script to find pseudogene candidates in annotated genbank files.
Tested mostly on .gbk files annotated by Prokka with the --compliant flag.

Installation requirements:
    python3
    Libraries: argparse, biopython

 """

from Bio.Blast.Applications import NcbiblastpCommandline, NcbiblastxCommandline
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from typing import NamedTuple, List
from time import localtime, strftime
from os import path
import argparse
import re
import csv
import sys

__author__ = "Mitch Syberg-Olsen & Filip Husnik"
__version__ = "0.06"
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
                               ('name', str)])

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

#TODO: figure out --help.
    # current error: [python pseudo_finder.py -h] "ValueError: unsupported format character ')' (0x29) at index 79"
def get_args():
    parser = argparse.ArgumentParser()
        #TODO: tab-complete doesnt work with this argument for some reason.
    parser.add_argument('-g', '--genome', help='Please provide your genome file in the genbank format.', required=True)

    parser.add_argument('-o', '--output', help='Specify an output prefix.', required=True)

    parser.add_argument('-d', '--database', help='Please provide name (if $BLASTB is set on your system) or absolute path of your blast database.')

    #TODO: discontinue this argument. XML is not supported in code
    parser.add_argument('-b', '--blast_format', help='Specify a blast output format, either \'xml\' or \'tsv\'.', default='tsv')

    parser.add_argument('-p', '--blastp', help='Specify an input blastp file.')

    parser.add_argument('-x', '--blastx', help='Specify an input blastx file.')

    parser.add_argument('-t', '--threads', help='Please provide total number of threads to use for blast, default is 4.', default=4)

    parser.add_argument('-i', '--intergenic_length', help='Please provide length of intergenic regions to check, default is 30 bp.', type=int, default=30)

    parser.add_argument('-l', '--length_pseudo',
                        help='Please provide percentage of length for pseudo candidates, default is 0.60 (60%).\nexample: \"-l 50\" will consider genes that are less than 50% of the average length of similar genes.',
                        type=float,
                        default=0.60)

    parser.add_argument('-s', '--shared_hits',
                        help='Percentage of blast hits that must be shared in order to join two nearby regions, default is 0.30 (30%).\nexample: \"-s 0.50\" will merge nearby regions if they shared 50% of their blast hits.',
                        type=float,
                        default=0.30)

    parser.add_argument('-e', '--evalue', help='Please provide e-value for blast searches, default is 1e-4.', default='1e-4')

    args = parser.parse_args()

    if args.blastx is None and args.blastp is None and args.database is None:
        parser.error("Pseudofinder requires a database input unless blast results are provided.")

    if (args.blastx is not None and args.blastp is None) or (args.blastp is not None and args.blastx is None):
        parser.error("Pseudofinder requires both blastP and blastX inputs.")

    return args


def get_proteome(gbk: str, faa: str) -> None:
    '''
    Parse genbank input file for coding sequences (CDSs) and write them to the output file with coordinates.
    '''

    with open(gbk, "r") as input_handle:
        with open(faa, "w") as output_handle:
            for seq_record in SeqIO.parse(input_handle, "genbank"):
                for seq_feature in seq_record.features:
                    if seq_feature.type == "CDS":
                        assert len(seq_feature.qualifiers['translation']) == 1
                        output_handle.write(">%s %s %s\n%s\n" % (seq_feature.qualifiers['locus_tag'][0],
                                                                 seq_record.name,
                                                                 seq_feature.location,
                                                                 seq_feature.qualifiers['translation'][0]))

    print('%s\tProteome extracted from:\t%s\n'
          '\t\t\tWritten to file:\t\t%s.' % (current_time(), gbk, faa,)),
    sys.stdout.flush()

#TODO: check how this works exactly, then clean
def get_intergenic_regions(gbk: str, fasta: str, igl: int) -> None:
    '''
    Parse genbank input file for intergenic regions and write them to the output file with coordinates.

    Copied/modified from "get_interregions" by Iddo Friedberg & Ian MC Fleming
    Released under Biopython license. http://www.biopython.org/DIST/LICENSE
    The original code extracts all regions strand-dependently, even if there is a gene on the other strand
    Such strand information is not needed here, so I arbitrarily select plus strand sequence
    '''
    #Parse all files in the multiple-file genbank
    for seq_record in SeqIO.parse(gbk, "genbank"):
        gene_list = []
        intergenic_records = []

        # Loop over the genome file, get the gene features on each of the strands
        for feature in seq_record.features:
            if feature.type == 'gene': #doesnt equal ribsomal and equals 'gene'
                mystart = feature.location._start.position
                myend = feature.location._end.position
                gene_list.append((mystart, myend, 1))

        for i, pospair in enumerate(gene_list[1:]):
            # Compare current start position to previous end position
            last_end = gene_list[i][1]
            this_start = pospair[0]
            strand = pospair[2]
            if this_start - last_end >= igl:
                intergene_seq = seq_record.seq[last_end:this_start]
                strand_string = "+"
                intergenic_records.append(SeqRecord(intergene_seq, id="%s-ign-%d" % (seq_record.name, i),
                                                    description="%s %d-%d %s" % (
                                                    seq_record.name, last_end + 1, this_start, strand_string)))

        SeqIO.write(intergenic_records, open(fasta, "a"), "fasta")

    print('%s\tIntergenic regions extracted from:\t%s\n'
          '\t\t\tWritten to file:\t\t%s.' % (current_time(), gbk, fasta,)),
    sys.stdout.flush()


def run_blastp(bf: str, faa: str, t: str, db: str, eval: str) -> None:
    '''
    run BLASTP with FAA file against DB of your choice.
    '''

    print('%s\tBlastP executed with %s threads.' % (current_time(),t)),
    sys.stdout.flush()

    if bf == 'xml':
        blastp_cline = NcbiblastpCommandline(query=faa,
                                             num_threads=t,
                                             db=db,
                                             num_alignments=15,
                                             evalue=eval,
                                             outfmt=5,
                                             out=faa + ".blastP_output.xml")
        blastp_cline()

    elif bf == 'tsv':
        blastp_cline = NcbiblastpCommandline(query=faa,
                                             num_threads=t,
                                             db=db,
                                             num_alignments=15,
                                             evalue=eval,
                                             outfmt= "\'7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send slen evalue bitscore\'",
                                             out=faa + ".blastP_output.tsv")
        blastp_cline()

    else:
        print('Function run_blastp() expected either \'xml\' or \'tsv\' as --blast_format, but received \'%s\'.\n' % bf)


def run_blastx(bf: str, fasta: str, t: str, db: str, eval: str) -> None:
    '''
    run BLASTX with FNA file against DB of your choice.
    '''

    print('%s\tBlastX executed with %s threads.' % (current_time(), t)),
    sys.stdout.flush()

    if bf == 'xml':
        blastx_cline = NcbiblastpCommandline(query=fasta,
                                             num_threads=t,
                                             db=db,
                                             num_alignments=15,
                                             evalue=eval,
                                             outfmt=5,
                                             out=fasta + ".blastX_output.xml")
        blastx_cline()

    elif bf == 'tsv':
        blastx_cline = NcbiblastxCommandline(query=fasta,
                                             num_threads=t,
                                             db=db,
                                             num_alignments=15,
                                             evalue=eval,
                                             outfmt="\'7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send slen evalue bitscore\'",
                                             out=fasta + ".blastX_output.tsv")
        blastx_cline()

    else:
        print('Function run_blastx() expected either \'xml\' or \'tsv\' as --blast_format.\n')


#try: parse locus and size from everything
#TODO: find the last locus tag number and continue from there
    #locus_tag

def make_gff_header(gbk: str, gff: str, blastp: str) -> None:
    '''
    Prepare a GFF file header (a list of sequence identifiers)
    '''

    #open file specified by input 'gff'
    with open(gff, "w") as gff_output_handle:

        #first line
        gff_output_handle.write("##gff-version 3\n"
                                "!annotation-date\t%s" % (current_time()))

        #writes one line for each contig
        for i, seq_record in enumerate(SeqIO.parse(gbk, "genbank")):

            gff_output_handle.write(
                "%s\t%s\t%s\t%s\n" % ("##sequence-region",
                                     seq_record.id,
                                      (i+1),
                                      (len(seq_record)+1)))


def collect_query_ids(filename: str) -> List[str]:
    '''
    Reads a TSV file and returns a list of Query names, to be used later.
    '''
    loq = []  # list of Query names

    with open(filename) as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            # Regex to match only lines that contain the Query ID
            if re.match("# Query:", row[0]):
                # Access Query ID by splitting a whole row based on whitespace, then selecting the 3rd field.
                query = row[0].split()[2]

                loq.append(query)
    return loq


def collect_hits(query: str, filename: str) -> List[BlastHit]:
    '''
    Parses a TSV file for hits that match the Query ID and stores them as a list.
    '''
    loh = []  # list of hits.

    with open(filename) as csvfile:

        reader = csv.reader(csvfile, delimiter='\t')

        for row in reader:
            # Regex to match only lines that contain the Query ID
            if re.match(query, row[0]):
                #In tabular blast files, rows that contain a given Query ID at the start will be followed by blast hits to that region.

                hit = BlastHit(row[1],          #BlastHit.accession
                               int(row[10])*3,  #BlastHit.slen - This is multiplied by 3 because the BlastP/BlastX files display slen in aa length.
                               int(row[8]),     #BlastHit.s_start
                               int(row[9]),     #BlastHit.s_end
                               float(row[11]))  #BlastHit.eval

                loh.append(hit)

    return loh


def get_contig(query: str, filename: str) -> str:
    '''
    Parses a TSV file for the range of a Query. The range is stored as [start:end](strand), ie [10:367](+).
    '''
    with open(filename) as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            # Regex to match a Query info line with a specific Query ID
            if re.match("# Query: %s" % (query), row[0]):
                # Access Query contig by splitting a whole row based on whitespace, then selecting the 4th field.
                query_contig = row[0].split()[3]
                return query_contig


def get_database_name(filename: str) -> str:
    '''
    Parses a TSV file for the name of the database used.
    '''
    with open(filename) as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            # Regex to match the line that says "Database"
            if re.match("# Database:", row[0]):
                # Access the database name by splitting a whole row based on whitespace, then selecting the 3rd field.
                database_name = row[0].split()[2]
                return database_name


def get_range(query: str, filename: str) -> str:
    '''
    Parses a TSV file for the range of a Query. The range is stored as [start:end](strand), ie [10:367](+).
    '''
    with open(filename) as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            # Regex to match a Query info line with a specific Query ID
            if 'blastP' in filename:
                if re.match("# Query: %s" % (query), row[0]):
                    # Access Query range by splitting a whole row based on whitespace, then selecting the 5th field.
                    query_range = row[0].split()[4]
                    return query_range

            elif 'blastX' in filename:
                if re.match("# Query: %s" % (query), row[0]):
                    # blastX output files look different than blastP, so this is a messy line to reformat it the same way.
                    query_range = str("[" + row[0].split()[4] + "]" + "(" + row[0].split()[5] + ")")
                    query_range = query_range.replace("-",":")
                    return query_range


def get_intergenic_query_range(lobh: List[BlastHit], coordinate: str, start_position: int) -> int:
    '''
    Calculates the range of an intergenic region, based on the location of blast hits within the whole intergenic region.
    This is necessary because the start and end positions of hits are defined locally in the blast output - this function
    converts them to absolute positions on the contig.
    '''

    #If there are no blast hits, the calculations will return an error.
    #This ensures that calculations are only done on lists with at least 1 item,
    if len(lobh) is not 0:

        #Collect all start positions in a list of blast hits
        start_list = [bh.s_start for bh in lobh]
        #Collect all end positions in a list of blast hits
        end_list = [bh.s_end for bh in lobh]

        if coordinate is 'start':
            # return the start position of the region itself + smallest start position of any blast hit to that region
            # this gives the absolute start position of this particular list of hits on the contig
            # Add 1 to convert being 0 based and 1 based numbering systems
            return start_position + min(start_list) + 1

        elif coordinate is 'end':
            # return the start position of the region itself + largest end position of any blast hit to that region
            # this gives the absolute end position of this particular list of hits on the contig
            # Add 1 to convert being 0 based and 1 based numbering systems
            return start_position + max(end_list) + 1

    else:
        return 0


def breakdown_range(range: str, i: str) -> int or str:
    '''
    Takes a str that looks like '[start:end](strand)' and extracts either the start coord, end coord, or strand.
    '''
    # Regex matches only the [start:end]
    coords = re.match("\[([0-9]*):([0-9]*)\]", range)
    # Regex matches only the '+' or '-'
    strand = re.search("\+|\-", range)

    if i is 'start_coord':
        # Add 1 to convert being 0 based and 1 based numbering systems
        return int(coords.group(1)) + 1
    elif i is 'end_coord':
        # Add 1 to convert being 0 based and 1 based numbering systems
        return int(coords.group(2)) + 1
    elif i is 'strand':
        return strand.group()


def get_regioninfo(filename: str, blast_format: str) -> List[RegionInfo]:
    '''
    This function needs to take a blast query and extract the relevant information (RegionInfo).
    '''

    print('%s\tExtracting information from %s file.' % (current_time(),blast_format)),
    sys.stdout.flush()

    query_ids = collect_query_ids(filename)
    lori = []  # Stores a list of RegionInfo

    if blast_format == 'BlastP':

        StatisticsDict['ProteomeOrfs'] += len(query_ids)

        for query in query_ids:
            #The following collects all necessary information for a single region, following RegionInfo data definition.

            region = RegionInfo(get_contig(query, filename),                                    #contig name
                                query,                                                          #query name
                                breakdown_range(get_range(query, filename), 'start_coord'),     #start position of region
                                breakdown_range(get_range(query, filename), 'end_coord'),       #end position of region
                                breakdown_range(get_range(query, filename), 'strand'),          #(+) or (-) strand
                                collect_hits(query, filename),                                  #List of blast hits for this region
                                'note: ORF not yet annotated as a pseudogene.')                     #note for storing info

            if len(region.hits) is not 0:
                lori.append(region)

    elif blast_format == 'BlastX':
        for query in query_ids:
            query_start = breakdown_range(get_range(query, filename), 'start_coord')
            lobh = collect_hits(query, filename)

            # The following collects all necessary information for a single region, following RegionInfo data definition.

            region = RegionInfo(get_contig(query, filename),                                #contig name
                                query,                                                      #query name
                                get_intergenic_query_range(lobh, 'start', query_start),     #start position of region
                                get_intergenic_query_range(lobh, 'end', query_start),       #end position of region
                                breakdown_range(get_range(query, filename), 'strand'),      #(+) or (-) strand
                                lobh,                                                       #List of blast hits for this region
                                'note: Intergenic region not yet annotated as a pseudogene.') #note for storing info

            if len(region.hits) is not 0:
                lori.append(region)

    return lori


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
        #this will store the List[RegionInfo] to be contained on a contig
        regions_on_contig = []

        for ri in lori:
            #if the region's contig name matches, it is added to that contig
            if ri.contig == contig_name:
                regions_on_contig.append(ri)

        #once all regions have been added, that list of regions is appended as a 'Contig' to the list of contigs.
        contig_list.append(Contig(regions_on_contig, contig_name))

    return contig_list


def number_of_matching_hits(r1: RegionInfo, r2: RegionInfo) -> int:
    '''
    This function returns the number of blast hits that two regions have in common.
    '''

    r1_accessions = set([blasthit.accession for blasthit in r1.hits])
    r2_accessions = set([blasthit.accession for blasthit in r2.hits])

    return len(set(r1_accessions) & set(r2_accessions))


def matching_hit_critera(r1: RegionInfo, r2: RegionInfo, cutoff: float) -> bool:
    '''
    This function determines if two regions meet the minimum blast hit criteria to be joined together.
    '''

    if len(r1.hits) is not 0 and len(r2.hits) is not 0:
        #sorts the two regions based on number of blast hits.
        s = sorted([r1, r2], key=lambda r: len(r.hits))

        # math is : (Number of shared hits) / (Total number of hits from the region with the least hits) >= cutoff value.
        if number_of_matching_hits(r1,r2)/len(s[0].hits) >= cutoff:
            return True
        else:
            return False

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


def compare_regions(r1: RegionInfo, r2: RegionInfo, cutoff: float) -> bool:
    '''
    Takes two regions and decides if they are similar enough to join together.
    '''

    if region_proximity(r1, r2) < 1000 and matching_hit_critera(r1, r2, cutoff) is True and r1.strand == r2.strand:
        return True

    else:
        return False


def join_regions(r1: RegionInfo, r2: RegionInfo) -> RegionInfo:
    '''
    This function needs to take two regions and merge their locations.
    '''

    # concatenates hits from both regions, discards any duplicates, and sorts them by e-value.
    merged_hits = sort_hits_by_eval(list(set(r1.hits + r2.hits)))

    merged_region = RegionInfo(r1.contig,                       #RegionInfo.contig
                               'locus_tag=pseudo_uniqueID',     #RegionInfo.query
                               r1.start,                        #RegionInfo.start
                               r2.end,                          #RegionInfo.end
                               r1.strand,                       #RegionInfo.strand
                               merged_hits,                     #RegionInfo.hits
                               'Note=pseudogene candidate. Reason: Predicted fragmentation of a single gene.;colour=229 204 255') #'colour=' makes this region appear coloured in Artemis.
    return merged_region


def sort_hits_by_eval(lobh: List[BlastHit]) -> List[BlastHit]:
    '''
    Sorts a list of blasthits by e-value from low to high (returning the hit with the lowest evalue first).
    '''

    sorted_list = sorted(lobh, key=lambda r: r.eval)

    return sorted_list


def annotate_pseudos(contig: Contig, contig_number: int, hits_cutoff: float, length_cutoff: float) -> List[RegionInfo]:
    '''
    This function will take input blast files and return a list of all pseudogene candidates.
    '''

    #1: Look through list of regions and find individual ORFs that could be pseudogenes.
    IndividualPseudos = check_individual_ORFs(contig.regions, contig_number, length_cutoff)

    #2: Update list of regions with any pseudogenes that were flagged from step #1.
    UpdatedList = replace_pseudos_in_list(IndividualPseudos, contig.regions)

    #3: Check adjacent regions to see if they could be pseudogene fragments.
    #   This function returns two lists: [0] = Individual pseudogenes
    #                                    [1] = Merged pseudogenes
    AllPseudos = check_adjacent_regions(UpdatedList, contig_number, hits_cutoff)

    #returns both individual and merged pseudogenes as a single list, with locus tags added.
    return add_locus_tags(AllPseudos[0] + AllPseudos[1], contig.name)


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

    return (False, gene)


def check_individual_ORFs(lori: List[RegionInfo], contig_number: int, length_cutoff: float) -> List[RegionInfo]:
    '''
    This function will take an input of regions and return a list of individual ORFs that could be pseudogenes.
    '''


    InitialList = [] #This list will contain all ORFs that have enough blast hits to be considered.
                     #If less than 3 blast hits, its hard to calculate a reliable "AverageDatabaseLength"

    lopg = [] #This list will contain the resulting pseudogenes

    for region in lori:
        if len(region.hits) > 2:
            InitialList.append(region)

    for region in InitialList:

        #Retrieves lengths of genes that this region has blasted against
        ListOfDatabaseLengths = [hit.slen for hit in region.hits]

        #Calculates the average length of genes that this region has blasted against
        AverageDatabaseLength = sum(ListOfDatabaseLengths) / len(ListOfDatabaseLengths)

        #Calculates the length of this region
        RegionLength = region.end - region.start

        #Ratio of the region's length to the average length of hits.
        Ratio = (RegionLength/AverageDatabaseLength)

        if Ratio < length_cutoff:
            print('%s\tIndividual gene flagged on contig %s, location %s-%s.' % (current_time(),
                                                                                 contig_number,
                                                                                 region.start,
                                                                                 region.end)),
            sys.stdout.flush()

            pseudo = convert_region_to_pseudo(region, Ratio*100)  #Multiplied by 100 to convert to percentage

            lopg.append(pseudo)

    return lopg


def convert_region_to_pseudo(region: RegionInfo, ratio: float) -> RegionInfo:
    '''
    Flags a region as a pseudogene by adding a note, that will appear in the GFF file.
    '''

    pseudogene = RegionInfo(region.contig,
                            'locus_tag=pseudo_uniqueID',
                            region.start,
                            region.end,
                            region.strand,
                            region.hits,
                            'Note=pseudogene candidate. Reason: ORF is %s%% of the average length of hits to this gene.;colour=229 204 255' % (round(ratio,1)))
                            #'colour=' makes this region appear coloured in Artemis.
    return pseudogene


def check_adjacent_regions(lori: List[RegionInfo], contig_number: int, cutoff: float) -> tuple:
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
        if compare_regions(sorted_lori[i], sorted_lori[i + 1], cutoff) is True:

            #this boolean will be important later on in this function
            NewPseudoMade = True

            #if they pass, create a pseudogene
            pseudo = join_regions(sorted_lori[i], sorted_lori[i + 1])

            #this is to keep track of overall statistics. If the regions have not yet been annotated by this program,
            #the counter will increase by 1 for each of them.
            for region in [sorted_lori[i], sorted_lori[i + 1]]:
                if 'Predicted fragmentation of a single gene' not in region.note:
                    StatisticsDict['FragmentedOrfs'] += 1

            #remove items that were joined together
            del sorted_lori[i + 1]
            del sorted_lori[i]

        #If regions [i] and [i+1] fail to join (above), look at regions [i] and [i+2].
        elif i < len(sorted_lori) - 2:
            if compare_regions(sorted_lori[i], sorted_lori[i + 2], cutoff) is True:

                # this boolean will be important later on in this function
                NewPseudoMade = True

                #if they pass, create a pseudogene
                pseudo = join_regions(sorted_lori[i], sorted_lori[i + 2])

                # this is to keep track of overall statistics. If the regions are have not yet been annotated by this program,
                # the counter will increase by 1 for each of them.
                for region in [sorted_lori[i], sorted_lori[i + 1], sorted_lori[i + 2]]:
                    if 'Predicted fragmentation of a single gene' not in region.note:
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
            print('%s\tRegions merged and flagged on contig %s, location %s-%s.' % (current_time(),
                                                                                    contig_number,
                                                                                    pseudo.start,
                                                                                    pseudo.end)),
            sys.stdout.flush()

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


def sort_contigs(loc: List[Contig]) -> List[Contig]:
    '''
    Takes a list of contigs and sorts it numerically.
    '''

    sortedlist = sorted(loc, key=lambda c: c.name)

    return sortedlist


def write_pseudos_to_gff(lopg: List[RegionInfo], gff: str) -> None:
    '''
    Takes an input list of Pseudogenes and writes them to a GFF file in proper format.
    '''

    print('%s\tWriting pseudogene candidates to GFF file.' % (current_time())),
    sys.stdout.flush()

    with open(gff, 'a') as gff_output_handle:
        for pseudo in lopg:
            gff_output_handle.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %
                                    (pseudo.contig,                 #1: seqname - name of the chromosome or scaffold
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


def write_summary_file(output_prefix, args) -> None:
    '''
    Writes a summary file of statistics from the pseudo_finder run.
    '''

    print('%s\tWriting summary of run.' % (current_time())),
    sys.stdout.flush()

    name = output_prefix + '_log.txt'

    with open(name, 'w') as logfile:
        logfile.write(
            "####### Summary from pseudo_finder.py #######\n\n"
            "Date/time:\t\t%s\n\n"
            #TODO: add runtime
            
            "#######    Files   #######\n"
            "Genome:\t\t\t%s\n"
            "Database:\t\t%s\n"
            "BlastP:\t\t\t%s\n"
            "BlastX:\t\t\t%s\n\n"
            
            "#######  Settings  #######\n"
            "Intergenic_length:\t%s\n"
            "Length_pseudo:\t\t%s\n"
            "Shared_hits:\t\t%s\n\n"
            
            "####### Statistics #######\n"
            "#Input:\n"
            "Initial ORFs:\t\t\t%s\n"
            "Number of contigs:\t\t%s\n"        #TODO: currently reporting smaller value than expected due to small contigs with no ORFs.
                                                #TODO: when intergenic parsing is fixed, check to make sure that this is resolved
            "#Output:\n"
            "Inital ORFs joined:\t\t%s\n"
            "Pseudogenes (total):\t\t%s\n"
            "Pseudogenes (too short):\t%s\n"
            "Pseudogenes (fragmented):\t%s\n"
            #'Functional genes' calculated as:  Input ORFs - [all ORFs marked as fragments] - [all ORFs marked as too short]
            "Functional genes:\t\t%s\n\n"

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


def main():
    args = get_args()

    #If blast files are not provided, must run blast.
    if args.blastp is None and args.blastx is None:

        #files generated:
        faa_filename = args.output + "_" + path.basename(args.genome) + "_proteome.faa"
        intergenic_filename = args.output + "_" + path.basename(args.genome) + "_intergenic.fasta"
        blastp_filename = faa_filename + ".blastP_output.tsv"
        blastx_filename = intergenic_filename + ".blastX_output.tsv"

        #Running blast
        get_proteome(args.genome, faa_filename)
        get_intergenic_regions(args.genome, intergenic_filename, args.intergenic_length)
        run_blastp(args.blast_format, faa_filename, args.threads, args.database, args.evalue)
        run_blastx(args.blast_format, intergenic_filename, args.threads, args.database, args.evalue)

    #If blast files are provided, use them.
    else:
        blastp_filename = args.blastp
        blastx_filename = args.blastx

    #BlastP and BlastX files have just been formally declared, so now we will add their names to the StatisticDict
    StatisticsDict['BlastpFilename'] = blastp_filename
    StatisticsDict['BlastxFilename'] = blastx_filename

    #Collect everything from the blast files
    #TODO: this was not actually parsing blastx files before, and when it does, things actually get buggy.
    #Investigate the bugs. Seen in one case to cause one fragment to be successfully joined but have multiple overlapping chunks also appear in gff.
    all_regions = get_regioninfo(blastp_filename, 'BlastP') #+ get_regioninfo(blastx_filename, 'BlastX')

    #Split into contigs
    all_contigs = sort_contigs(split_regions_into_contigs(all_regions))
    StatisticsDict['NumberOfContigs'] = len(all_contigs)

    # Write header
    gff_filename = args.output + "_" + path.basename(args.genome) + "_pseudo_finder.gff"
    make_gff_header(args.genome, gff_filename, blastp_filename)

    #For each contig:
    for contig in all_contigs:

        # Reports the contig number being examined.
        print('%s\tChecking contig %s / %s for pseudogenes.' % (current_time(),
                                                                all_contigs.index(contig)+1,  #indices are 0 based so I added 1 to make it more intuitive
                                                                len(all_contigs))),
        sys.stdout.flush()

        #Annotate pseudos then write them to the GFF file
        write_pseudos_to_gff(annotate_pseudos(contig, all_contigs.index(contig)+1, args.shared_hits, args.length_pseudo), gff_filename)

    write_summary_file(args.output, args)



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
