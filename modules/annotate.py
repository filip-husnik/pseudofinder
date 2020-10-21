#!/usr/bin/env python3

from . import common
import re
import sys
import os
from copy import deepcopy, copy
from enum import Enum
from typing import NamedTuple, List, OrderedDict
from time import localtime, strftime
import statistics

from Bio.Blast.Applications import NcbiblastpCommandline, NcbiblastxCommandline
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
from collections import defaultdict

# This try block was added to stop a circular import error that occurs when this module is called from reannotate.py
try:
    from . import genome_map
except ImportError:
    pass

# importing dnds module
try:
    from . import dnds
except ImportError:
    pass

# Data definitions
# An individual blast hit to a region.
BlastHit = NamedTuple('BlastHit', [('blast_type', str),
                                   ('query', str),
                                   ('subject_accession', str),
                                   ('percent_ident', float),
                                   ('length', int), #reported in aa for blastp and nt for blastx. see blasthit_length()
                                   ('mismatch', int),
                                   ('gapopen', int),
                                   ('q_start', int),
                                   ('q_end', int),
                                   ('s_start', int),
                                   ('s_end', int),
                                   ('evalue', float),
                                   ('bitscore', float),
                                   ('stitle', str)])

dnds_data = NamedTuple('dnds_data', [('locus_tag', str),
                           ('reference_tag', str),
                           ('dn', float),
                           ('ds', float),
                           ('dnds', float)])

# All possible types of regions
RegionType = Enum('RegionType', ['ORF',
                                 'intergenic',
                                 'shortpseudo',
                                 'fragmentedpseudo',
                                 'intergenicpseudo'])

PseudoType = Enum('PseudoType', ['short',
                                 'fragmented',
                                 'intergenic',
                                 'consumed',
                                 'dnds'])

RegionInfo = NamedTuple('RegionInfo', [('contig', str),
                                       ('query', str),
                                       ('genbank_locus_tags', list),
                                       ('pseudo_locus_tag', str),
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
                    'PseudogenesIntergenic': 0,
                    'OutputFiles': [],
                    'dnds': 0,
                    'IntactORFs': 0
                  }


# Functions
def current_time() -> str:
    """Returns the current time when this function was executed."""
    return str(strftime("%Y-%m-%d %H:%M:%S", localtime()))


def add_intergenic_to_seqrecord(args, seqrecord):
    """
    Parses a seqrecord for intergenic space and adds each intergenic region as a feature on the record.
    """
    gene_list = []  # List of coding regions extracted from genbank file.

    for feature in seqrecord.features:
        if feature.type == 'gene':
            start_position = feature.location.start
            end_position = feature.location.end
            gene_list.append((start_position, end_position))

    if args.contig_ends is True:    # Forces the program to consider space between the start and the first gene and same for the end of the contig
        gene_list.insert(index=0, object=(0,0))
        contig_end = len(seqrecord.seq)
        gene_list.append((contig_end, contig_end))

    for i, gene in enumerate(gene_list):
        last_end = gene_list[i - 1][1]
        this_start = gene_list[i][0]

        if this_start - last_end >= args.intergenic_length:
            intergenic_region = SeqFeature(location=FeatureLocation(last_end, this_start),
                                           type='intergenic',
                                           strand=0)

            seqrecord.features.append(intergenic_region)


def add_qualifiers_to_features(args, seqrecord):
    """
    Adds additional qualifiers to each feature contained in the seqrecord, necessary for further analysis in the pipeline.
    """
    intergenic_counter = 1
    for seq_feature in seqrecord.features:
        if seq_feature.type != 'source':
            seq_feature.qualifiers['nucleotide_seq'] = seq_feature.extract(seqrecord.seq)
            seq_feature.qualifiers['contig_id'] = seqrecord.name
            seq_feature.qualifiers['hits'] = []
            seq_feature.qualifiers['pseudo_type'] = None
            seq_feature.qualifiers['note'] = ''

        if seq_feature.type == 'intergenic':
            seq_feature.qualifiers['locus_tag'] = ['%s_ign_%05d' % (seqrecord.name, intergenic_counter)]
            intergenic_counter += 1

        if seq_feature.type == 'CDS':
            seq_feature.qualifiers['dnds'] = []



def gbk_to_seqrecord_list(args, gbk: str):
    """
    Uses biopython SeqIO to convert a genbank file into a list of SeqRecord, one for each contig then modified.
    - intergenic space is added as features
    - additional qualifiers added to features
    """
    seqrecord_list = []
    with open(gbk, "r") as input_handle:
        for record in SeqIO.parse(input_handle, "genbank"):
            seqrecord_list.append(record)

    for contig in seqrecord_list:
        add_intergenic_to_seqrecord(args, contig)
        add_qualifiers_to_features(args, contig)
        contig.features = sorted(contig.features, key=lambda x: x.location.start)

    return seqrecord_list


def parse_features_from_record(record: SeqRecord, feature_type: str) -> list:
    """
    Returns all features of a particular type from a record as a list of seq_features.
    """
    feature_list = []
    for seq_feature in record.features:
        if seq_feature.type == feature_type:
            feature_list.append(seq_feature)
    return feature_list


def extract_features_from_genome(args, genome, feature_type: str) -> list:
    """
    Returns all features from genome as a list.
    """
    feature_list = []
    for record in genome:
        feature_list.extend(parse_features_from_record(record, feature_type))
    return feature_list


def write_fasta(seqs: list, outfile: str, seq_type='nt') -> None:
    """Takes a list of sequences in seq_feature format and writes them to fasta file"""

    with open(outfile, "w") as output_handle:
        for seq in seqs:
            if seq_type == 'nt':
                seq_string = seq.qualifiers['nucleotide_seq']
            elif seq_type == 'aa':
                seq_string = seq.qualifiers['translation'][0]
            else:
                print("Invalid seq_type. Please check your write_fasta() function call.")
                exit()

            output_handle.write(">%s %s %s\n%s\n" % (seq.qualifiers['locus_tag'][0],
                                                     seq.qualifiers['contig_id'],
                                                     seq.location,
                                                     seq_string))


def gff_entry(seq_id, source, feature_type, feature_start, feature_end, score, strand, phase, attributes):
    return "\t".join([seq_id, source, feature_type, feature_start, feature_end, score, strand, phase, attributes])


def format_strand(strand, format='gff'):
    """Formats strand information namely for GFF3 compliance."""
    if format == 'gff':
        if strand == 1:
            return "+"
        elif strand == -1:
            return "-"
        elif strand == 0 or strand is None:
            return "."

    elif format == 'biopython':
        if strand == '+':
            return 1
        elif strand == "-":
            return -1
        elif strand == ".":
            return 0

    else:
        raise RuntimeError("incorrect format for format_strand().")


def attribute_string(feature):
    """Generates a string to fill attribute section of a GFF file."""
    list = []
    if feature.qualifiers['note']:
        list.append("note=" + feature.qualifiers['note'])
    list.append("locus_tag=" + feature.qualifiers['locus_tag'][0])
    return ";".join(list)


def write_gff(args, genome, outfile: str, seq_type):
    """Writes a GFF3 format compliant GFF file with the desired sequence type."""

    with open(outfile, 'w') as outfile:
        # Header
        outfile.write("##gff-version 3\n#!annotation-date\t%s\n" % (common.current_time()))
        for seqrecord in genome:
            entry = ["##sequence-region",
                     "gnl|Prokka|%s" % seqrecord.name,
                     str(1),
                     str(len(seqrecord))]
            outfile.write(" ".join(entry) + '\n')

        seqs = extract_features_from_genome(args, genome, seq_type)
        for seq in seqs:
            seq_info = {'seq_id': "gnl|Prokka|%s" % seq.qualifiers['contig_id'],
                        'source': "pseudofinder",
                        'feature_type': seq_type,
                        'feature_start': str(seq.location.start),
                        'feature_end': str(seq.location.end),
                        'score': '.',
                        'strand': format_strand(seq.strand),
                        'phase': '.',   # TODO: Understand what phase is and how to properly declare this
                        'attributes': attribute_string(seq)}

            outfile.write(gff_entry(**seq_info) + "\n")



def write_gbk(args, genome, outfile):
    """Writes a genbank file, excluding any analysis qualifiers present in the data structure."""

    cleared_genome = clear_analysis_qualifiers(args, genome)

    with open(outfile, "w") as output:
        SeqIO.write(cleared_genome, output, "genbank")


def get_CDSs(gbk: str, out_fasta: str) -> None:
    """Parse genbank input file for coding sequences (CDSs) and write the
    nucleotide sequences to the output file with coordinates."""

    with open(gbk, "r") as input_handle:
        with open(out_fasta, "w") as output_handle:
            for seq_record in SeqIO.parse(input_handle, "genbank"):
                for seq_feature in seq_record.features:
                    if seq_feature.type == "CDS":
                        assert len(seq_feature.qualifiers['translation']) == 1
                        output_handle.write(">%s %s %s\n%s\n" % (seq_feature.qualifiers['locus_tag'][0],
                                                                 seq_record.name,
                                                                 seq_feature.location,
                                                                 seq_feature.extract(seq_record.seq)))

    print('%s\tCDS extracted from:\t\t\t%s\n'
          '\t\t\tWritten to file:\t\t\t%s.' % (current_time(), gbk, out_fasta,)),
    sys.stdout.flush()


def get_proteome(gbk: str, out_faa: str) -> None:
    """Parse genbank input file for coding sequences (CDSs) and write them to the output file with coordinates."""

    with open(gbk, "r") as input_handle:
        with open(out_faa, "w") as output_handle:
            for seq_record in SeqIO.parse(input_handle, "genbank"):
                for seq_feature in seq_record.features:
                    if seq_feature.type == "CDS":
                        assert len(seq_feature.qualifiers['translation']) == 1
                        output_handle.write(">%s %s %s\n%s\n" % (seq_feature.qualifiers['locus_tag'][0],
                                                                 seq_record.name,
                                                                 seq_feature.location,
                                                                 seq_feature.qualifiers['translation'][0]))

    print('%s\tProteome extracted from:\t\t%s\n'
          '\t\t\tWritten to file:\t\t\t%s.' % (current_time(), gbk, out_faa,)),
    sys.stdout.flush()


def get_intergenic_regions(args, out_fasta: str) -> None:
    """Parse genbank input file for intergenic regions and write them to the output file with coordinates.

    Copied/modified from "get_interregions" by Iddo Friedberg & Ian MC Fleming
    Released under Biopython license. http://www.biopython.org/DIST/LICENSE
    The original code extracts all regions strand-dependently, even if there is a gene on the other strand
    Such strand information is not needed here, so I arbitrarily select plus strand sequence."""

    # Resets 'fasta' if it contains content already
    open(out_fasta, 'w').close()

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
        SeqIO.write(intergenic_records, open(out_fasta, "a"), "fasta")

    print('%s\tIntergenic regions extracted from:\t%s\n'
          '\t\t\tWritten to file:\t\t\t%s.' % (current_time(), args.genome, out_fasta,)),
    sys.stdout.flush()


def run_blast(args, search_type: str, in_fasta: str, out_tsv: str) -> None:
    """"Run BLASTP or BLASTX with fasta file against DB of your choice."""

    common.print_with_time("%s executed with %s threads." % (search_type, args.threads))
    # print('%s\t%s executed with %s threads.' % (current_time(), search_type, args.threads)),
    # sys.stdout.flush()

    blast_dict = {'query': in_fasta,
                  'num_threads': args.threads,
                  'db': args.database,
                  'max_target_seqs': args.hitcap,
                  'max_hsps': 1,
                  'evalue': args.evalue,
                  'outfmt': "6 qseqid sseqid pident slen mismatch gapopen qstart qend sstart send evalue bitscore stitle",
                  'out': out_tsv}

    if search_type == 'blastp':
        NcbiblastpCommandline(**blast_dict)()
    if search_type == 'blastx':
        NcbiblastxCommandline(**blast_dict)()


def manage_diamond_db(args):
    # try:
    #     open(args.database + ".dmnd")
    #     print("Found DMND database! " + args.database + ".dmnd")
    # except FileNotFoundError:
    os.system('diamond makedb --in %s --db %s.dmnd > /dev/null 2>&1' % (args.database, args.database))


def manage_blast_db(args):
    # try:
    #     open(args.database + ".psq")
    #     open(args.database + ".phr")
    #     open(args.database + ".pin")
    #     print("Found BLAST database! " + args.database)
    # except FileNotFoundError:
    os.system('makeblastdb -dbtype prot -in %s -out %s > /dev/null 2>&1' % (args.database, args.database))


def run_diamond(args, search_type: str, in_fasta: str, out_tsv: str) -> None:
    if search_type == 'blastp' or search_type == 'blastx':

        common.print_with_time("Diamond %s executed with %s threads." % (search_type, args.threads))

        diamond_cline = ('diamond %s --quiet --query %s --out %s --threads %s --max-target-seqs %s --evalue %s --db %s '
                         '--outfmt 6 qseqid sseqid pident slen mismatch gapopen qstart qend sstart send evalue bitscore stitle --max-hsps 1'
                         % (search_type, in_fasta, out_tsv, args.threads, args.hitcap, args.evalue, args.database))
        os.system(diamond_cline)

    else:
        print("function run_diamond() can only accept \'blastp\' or \'blastx\' as search_type. Please use one of these options.")
        exit()


def convert_tsv_to_blasthits(blast_tsv, blast_type):
    """
    Parses a TSV formatted blast file and returns each hit as a BlastHit object in a list.
    """
    hit_list = []
    with open(blast_tsv, 'r') as tsv:
        for line in tsv.readlines():
            fields = [blast_type] + [common.literal_eval(x) for x in re.split("\t", line.rstrip("\n"))]
            hit_list.append(BlastHit(*fields))

    return hit_list


def blasthit_length(hit, seq_type='nt', query=False):
    """Returns the length of the given blasthit in the preferred format."""
    if query:
        if hit.blast_type == 'blastp':
            aa_length = (hit.q_end - hit.q_start)
            nt_length = aa_length*3

        elif hit.blast_type == 'blastx':
            nt_length = (hit.q_end - hit.q_start)
            aa_length = nt_length/3

        else:
            raise RuntimeError("unknown blasthit type.")

    else:
        if hit.blast_type == 'blastp':
            aa_length = hit.length
            nt_length = aa_length*3

        if hit.blast_type == 'blastx':
            nt_length = hit.length
            aa_length = nt_length/3

    if seq_type == 'nt':
        return nt_length
    elif seq_type == 'aa':
        return aa_length
    else:
        raise RuntimeError("unknown seq_type.")


def features_with_locus_tags(genome):
    """
    Returns only features that have been given a locus tag.
    """
    relevant_features = []
    for seqrecord in genome:
        relevant_features.extend([feature for feature in seqrecord.features if feature.qualifiers.get('locus_tag')])
    return relevant_features


def add_blasthits_to_genome(args, genome, blast_file, blast_type):
    """
    Reads a TSV formatted blast output and adds all hits to record.features.qualifiers['hits'], where:
    record = contig
    feature = CDS or intergenic region
    """

    hit_list = convert_tsv_to_blasthits(blast_file, blast_type)
    relevant_features = features_with_locus_tags(genome)
    feature_dict = {feature.qualifiers['locus_tag'][0]: feature for feature in relevant_features}

    for hit in hit_list:
        feature_dict[hit.query].qualifiers['hits'].append(hit)


def convert_csv_to_dnds(dnds_file):

    dnds_list = []
    with open(dnds_file, 'r') as csv:
        next(csv)
        for line in csv.readlines():
            fields = [common.literal_eval(x) for x in re.split(",", line.rstrip("\n"))]
            dnds_list.append(dnds_data(*fields[:-1]))
    return dnds_list


def add_dnds_info_to_genome(args, genome, dnds_folder):
    """Adds dN/dS info from analysis to genome"""

    dnds_list = convert_csv_to_dnds(dnds_folder+"/dnds-summary.csv")
    relevant_features = extract_features_from_genome(args, genome, 'CDS')
    feature_dict = {feature.qualifiers['locus_tag'][0]: feature for feature in relevant_features}
    for item in dnds_list:
        feature_dict[item.locus_tag].qualifiers['dnds'].append(item)


def update_intergenic_locations(args, genome):
    """Uses blasthit information to update where relevant intergenic space is on the genome."""

    intergenic = extract_features_from_genome(args, genome, 'intergenic')

    for feature in intergenic:
        hits = feature.qualifiers['hits']
        if len(hits) > 0:
            all_starts = [hit.q_start-1 for hit in hits]
            all_ends = [hit.q_end-1 for hit in hits]
            all_vals = all_starts + all_ends
            absolute_start = feature.location.start
            new_start = min(all_vals)
            new_end = max(all_vals)
            feature.location = FeatureLocation(absolute_start + new_start, absolute_start + new_end)
            feature.qualifiers['nucleotide_seq'] = feature.qualifiers['nucleotide_seq'][new_start:new_end]


def feature_length_relative_to_hits(feature) -> float:
    """
    Returns a float that is equal to the ratio of (feature length / average length of hits)
    """
    db_lengths = [blasthit_length(hit, 'nt') for hit in feature.qualifiers['hits']]
    average_db_len = sum(db_lengths) / len(db_lengths)
    return len(feature) / average_db_len


def find_individual_pseudos(args, seqrecord):

    if args.reference:
        find_dnds_pseudos(args, seqrecord)

    for feature in seqrecord.features:
        if feature.type == 'CDS' and len(feature.qualifiers['dnds']) > 0:
            dnds_val = feature.qualifiers['dnds'][0].dnds

            if dnds_val >= args.max_dnds:
                feature.type = 'pseudogene'
                feature.qualifiers['pseudo_type'] = PseudoType.dnds
                feature.qualifiers['note'] = 'Pseudogene candidate. Reason: Elevated dN/dS (%s).' % round(dnds_val, 3)

        if feature.type == 'CDS' and len(feature.qualifiers['hits']) > 0:
            if feature_length_relative_to_hits(feature) <= args.length_pseudo:
                feature.type = 'pseudogene'
                feature.qualifiers['pseudo_type'] = PseudoType.short
                feature.qualifiers['note'] = 'Pseudogene candidate. Reason: ORF is %s%% of the average length of ' \
                                             'hits to this gene.' % (round(feature_length_relative_to_hits(feature)*100, 1))

        if feature.type == 'intergenic':
            if len(feature.qualifiers['hits']) / args.hitcap >= args.intergenic_threshold:
                feature.type = 'pseudogene'
                feature.qualifiers['pseudo_type'] = PseudoType.intergenic
                feature.qualifiers['note'] = 'Pseudogene candidate. Reason: Intergenic region with %s blast hits.' % len(feature.qualifiers['hits'])
        else:
            pass


def find_dnds_pseudos(args, seqrecord):
    return False


def matching_hits(f1, f2):
    """Returns percentage of hits shared by two regions."""
    f1_hits = f1.qualifiers['hits']
    f2_hits = f2.qualifiers['hits']

    if len(f1_hits) == 0 or len(f2_hits) == 0:
        return False

    else:
        f1_set = set([hit.subject_accession for hit in f1_hits])
        f2_set = set([hit.subject_accession for hit in f2_hits])

        num_matching_hits = len(f1_set & f2_set)
        least_hits_in_comparison = min(len(x) for x in (f1_hits, f2_hits))

    return num_matching_hits / least_hits_in_comparison


def matching_strands(f1, f2):
    """
    Returns whether two features are on the same strand. Always return true if one feature is intergenic
    """
    if f1.type == 'intergenic' or f2.type == 'intergenic':
        return True
    else:
        return f1.strand == f2.strand


def adjacent_fragments_match(args, f1, f2):
    """Checks if features are fragments of a single pseudogene"""

    if (
        f1.location.end - f2.location.start < args.distance and
        matching_hits(f1, f2) >= args.shared_hits and
        matching_strands(f1, f2)
    ):
        return True

    else:
        return False


def create_fragmented_pseudo(args, fragments, seqrecord):
    """Takes a list of features are concatenates them into a single pseudogene feature"""

    start = min([feature.location.start for feature in fragments])
    end = max([feature.location.end for feature in fragments])

    strands = [feature.strand for feature in fragments if (feature.strand is not None and feature.strand != 0)]
    if len(strands) == 0:   # Occurs if two intergenic regions are being merged
        strand = 0
    elif all([strand == strands[0] for strand in strands]):
        strand = strands[0]
    else:
        raise RuntimeError("Trying to combine genes on opposite strands.")

    hits = []
    for feature in fragments:
        hits.extend(feature.qualifiers['hits'])

    parents = []
    for fragment in fragments:
        frag_parents = fragment.qualifiers.get('parents')
        if frag_parents is not None:
            parents.extend(frag_parents)
        parents.extend([feature.qualifiers['locus_tag'][0] for feature in fragments])
        parents = list(set(parents))


    pseudo = SeqFeature(location=FeatureLocation(start, end),
                        type='pseudogene',
                        strand=strand)

    pseudo.qualifiers['nucleotide_seq'] = pseudo.extract(seqrecord.seq)
    pseudo.qualifiers['contig_id'] = seqrecord.name
    pseudo.qualifiers['hits'] = hits
    pseudo.qualifiers['locus_tag'] = ''
    pseudo.qualifiers['parents'] = parents
    pseudo.qualifiers['pseudo_type'] = PseudoType.fragmented
    pseudo.qualifiers['note'] = "Pseudogene candidate. Reason: Predicted fragmentation of a single gene."

    seqrecord.features.append(pseudo)


def manage_parent_fragments(fragments):
    """
    Adds an identifier to parent fragments in the list so that they will not be used again in the analysis.
    """
    for feature in fragments:
        feature.type = 'consumed'
        feature.qualifiers['pseudo_type'] = PseudoType.consumed


def find_fragmented_pseudos(args, seqrecord):
    """Finds pseudogenes composed of nearby fragments"""

    features = [feature for feature in seqrecord.features if feature.type == 'CDS' or feature.type == 'intergenic' or feature.type == 'pseudogene']  # initial count

    if len(features) > 2:   # contig must have at least 3 features in order to run this algorithm
        i = 0
        while i < len(features)-2:
            f1 = features[i]
            f2 = features[i+1]
            f3 = features[i+2]

            if adjacent_fragments_match(args, f1, f2):
                create_fragmented_pseudo(args, (f1, f2), seqrecord)
                manage_parent_fragments((f1, f2))
                i -= 1

            elif adjacent_fragments_match(args, f1, f3):
                create_fragmented_pseudo(args, (f1, f2, f3), seqrecord)
                manage_parent_fragments((f1, f2, f3))
                i -= 1

            else:
                i += 1
                continue

            features = [feature for feature in features if feature.qualifiers.get('pseudo_type') is not PseudoType.consumed]
            features = sorted(features, key=lambda x: x.location.start)

        # final sort when finished
        seqrecord.features = sorted(seqrecord.features, key=lambda x: x.location.start)


def update_locus_tags(args, genome):
    """Gives all pseudogenes a unique locus tag."""

    pseudos = extract_features_from_genome(args, genome, 'pseudogene')
    i = 1
    for feature in pseudos:
        contig_id = feature.qualifiers['contig_id']
        feature.qualifiers['locus_tag'] = ['%s_pseudo_%05d' % (contig_id, i)]
        i += 1


def find_pseudos_on_genome(args, genome):
    """
    The main pseudogene function which performs individual pseudogene analyses.
    """

    for i, seqrecord in enumerate(genome):
        common.print_with_time("Checking contig %s / %s for pseudogenes." % (i+1, len(genome)))
        find_individual_pseudos(args, seqrecord)
        find_fragmented_pseudos(args, seqrecord)

        num_ORFs = len(parse_features_from_record(seqrecord, 'CDS'))
        pseudos = [x for x in parse_features_from_record(seqrecord, 'pseudogene') if x.qualifiers.get("pseudo_type") != PseudoType.consumed]
        num_pseudos = len(pseudos)

        common.print_with_time("Number of ORFs on this contig: %s\n"
                               "\t\t\tNumber of pseudogenes flagged: %s" % (num_ORFs, num_pseudos))

    update_locus_tags(args, genome)




def write_test_genome_output(file_dict, genome): #TODO: maybe move this to common or sandbox when done
    with open(file_dict['base_filename']+"test_genome.gbk", "w") as output_handle:
        SeqIO.write(genome, output_handle, "genbank")


def parse_blast(fasta_file: str, blast_file: str, blast_format: str) -> List[RegionInfo]:
    """This function needs to take a blast query and extract the relevant information (RegionInfo)."""
    print('%s\tExtracting information from %s file.' % (current_time(), blast_format)),
    sys.stdout.flush()

    query_dict = {}  # Dictionary of information relating to each query

    # Read fasta file, build query dict without knowledge of blasthits
    with open(fasta_file, 'r') as fasta:
        lines = [line for line in fasta.readlines() if re.match("^>", line)]
        if blast_format == 'blastp':
            StatisticsDict['ProteomeOrfs'] = len(lines)
        for line in lines:
            fields_in_line = list(filter(None, re.split("\s|(?<=[0-9])-|\[|\]|:|\(|\)", line)))
            query = fields_in_line[0][1:]
            query_dict[query] = {'contig': fields_in_line[1],
                                 'query': query,
                                 'start': int(fields_in_line[2]),
                                 'end': int(fields_in_line[3]),
                                 'strand': fields_in_line[4],
                                 'hits': []}

    # Iterate through the tsv file, adding information from each line to the appropriate query
    with open(blast_file, 'r') as tsv:
        lines = tsv.readlines()
        for line in lines:
            fields_in_line = list(filter(None, re.split("\s|(?<=[0-9])-|\[|\]|:|\(|\)", line)))
            line_query = fields_in_line[0]
            blast_hit = BlastHit(query=line_query,
                                 subject_accession=fields_in_line[1],
                                 percent_ident=float(fields_in_line[2]),
                                 aa_length=int(fields_in_line[3])*3, # convert aa to nucleotide
                                 mismatch=int(fields_in_line[4]),
                                 gapopen=int(fields_in_line[5]),
                                 q_start=int(fields_in_line[6]),
                                 q_end=int(fields_in_line[7]),
                                 s_start=int(fields_in_line[8]),
                                 s_end=int(fields_in_line[9]),
                                 evalue=float(fields_in_line[10]),
                                 bitscore=float(fields_in_line[11]),
                                 stitle=" ".join(fields_in_line[12:]))

            query_dict[line_query]['hits'].append(blast_hit)

    # Write all of the information to a RegionInfo, which is then added to a list of RegionInfo
    if blast_format == 'blastp':
        note = "From BlastP;colour=51 153 102"
        region_type = RegionType.ORF
    elif blast_format == 'blastx':
        note = "From BlastX"
        region_type = RegionType.intergenic

    region_list = []
    for query in query_dict:
        region_list.append(RegionInfo(contig=query_dict[query]['contig'],
                                      query=query,
                                      genbank_locus_tags=[query],
                                      pseudo_locus_tag="",
                                      start=query_dict[query]['start'],
                                      end=query_dict[query]['end'],
                                      strand=query_dict[query]['strand'],
                                      hits=query_dict[query]['hits'],
                                      note=note,
                                      region_type=region_type))
    return region_list


def split_regions_into_contigs(lori: List[RegionInfo]) -> List[Contig]:
    """Takes a list of regions and splits them based on which contig it belongs to.
    Contig is defined above as 'List[RegionInfo]', so 'List[Contig]' is a list of lists."""

    # collects all contig names. Doesn't store duplicates, so a contig name will not be stored more than once.
    contig_names = list(set([ri.contig for ri in lori]))
    StatisticsDict['NumberOfContigs'] = len(contig_names)
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


def annotate_pseudos(args, contig: Contig) -> Contig:
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
    all_pseudos = check_adjacent_regions(args=args, lori=updated_list)

    final_regions = add_locus_tags(lori=(all_pseudos[0] + all_pseudos[1]), contig=contig.name)

    # returns both individual and merged pseudogenes as a single list, with locus tags added.
    return Contig(regions=final_regions, name=contig.name, number=contig.number)


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
        elif region.region_type == RegionType.intergenic and len(region.hits)/int(args.hitcap) >= args.intergenic_threshold:
            pseudo = convert_region_to_pseudo(region=region,
                                              ratio=None,  # this value is only used for BlastP-derived pseudos (below)
                                              number_of_hits=len(region.hits))
            blastx_pseudos.append(pseudo)

    for region in initial_blastp_list:

        # Retrieves lengths of genes that this region has blasted against
        list_of_database_lengths = [hit.nucleotide_length for hit in region.hits]

        # Calculates the average length of genes that this region has blasted against
        average_database_length = sum(list_of_database_lengths) / len(list_of_database_lengths)

        # Calculates the length of this region
        region_length = (region.end - region.start)

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
                            genbank_locus_tags=region.genbank_locus_tags,
                            pseudo_locus_tag="",
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
    sorted_by_start = sorted([r1, r2], key=lambda r: r.start)

    # substracts the end position of the first from the start position of the second
    # this value can actually be negative if a gene starts before the previous one finishes
    return sorted_by_start[1].start - sorted_by_start[0].end


def matching_hit_critera(args, r1: RegionInfo, r2: RegionInfo) -> bool:
    """This function determines if two regions meet the minimum blast hit criteria to be joined together."""

    if len(r1.hits) != 0 and len(r2.hits) != 0:
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

    r1_accessions = set([blasthit.subject_accession for blasthit in r1.hits])
    r2_accessions = set([blasthit.subject_accession for blasthit in r2.hits])

    return len(set(r1_accessions) & set(r2_accessions))


def join_regions(r1: RegionInfo, r2: RegionInfo) -> RegionInfo:
    """This function needs to take two regions and merge their locations."""

    # concatenates hits from both regions, discards any duplicates, and sorts them by e-value.
    merged_hits = sort_hits_by_eval(list(set(r1.hits + r2.hits)))

    merged_region = RegionInfo(contig=r1.contig,
                               query=r1.query+","+r2.query+",",
                               genbank_locus_tags=r1.genbank_locus_tags + r2.genbank_locus_tags,
                               pseudo_locus_tag="",
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

    sorted_list = sorted(lobh, key=lambda r: r.evalue)
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
        tagged_region = RegionInfo(contig=region.contig,
                                   query=region.query,
                                   genbank_locus_tags=region.genbank_locus_tags,
                                   pseudo_locus_tag=str("%s_%04d" % (contig, counter+1)),
                                   start=region.start,
                                   end=region.end,
                                   strand=region.strand,
                                   hits=region.hits,
                                   note=region.note,
                                   region_type=region.region_type)

        final_list.append(tagged_region)

    return final_list


def write_genes_to_gff(args, lopg: List[RegionInfo], gff: str) -> None:
    """Takes an input list of genes and writes them to a GFF file in proper format."""

    with open(gff, 'w') as gff_output_handle:
        # write header
        gff_output_handle.write("##gff-version 3\n#!annotation-date\t%s\n" % (current_time()))  # first line
        for i, seq_record in enumerate(SeqIO.parse(args.genome, "genbank")):  # writes one line for each contig
            entry_elements = ["##sequence-region",                # Necessary to comply with GFF3 formatting
                              "gnl|Prokka|%s" % seq_record.name,    # contig seqid
                              1,                                  # contig start
                              len(seq_record)]

            gff_output_handle.write(' '.join(map(str, entry_elements))+'\n')

        # write genes
        for region in lopg:
            locus_tag = "locus_tag=%s" % region.pseudo_locus_tag
            genbank_locus_tags = "gbk_locus_tags=%s" % ",".join(region.genbank_locus_tags)
            attributes = ";".join([region.note, locus_tag, genbank_locus_tags])
            entry_elements = ["gnl|Prokka|%s" % region.contig,
                              "pseudofinder",
                              "gene",
                              region.start,
                              region.end,
                              '.',
                              region.strand,
                              '.',
                              attributes]

            gff_output_handle.write('\t'.join(map(str, entry_elements))+'\n')


def get_intact_genes(contig: Contig, pseudos: List[RegionInfo]) -> Contig:
    """"Inspects a contig for genes that have not been annotated as pseudogenes, and returns them."""

    # All regions on a contig, sorted by start position
    region_list = sorted(contig.regions, key=lambda r: r.start)
    # Begin with all regions on the contig
    intact_list = sorted(contig.regions, key=lambda r: r.start)

    # Iterate through regions on the contig
    for region in region_list:
        for pseudo in pseudos:
            # This will be true if the region is nested within a pseudogene
            if region.start >= pseudo.start and region.end <= pseudo.end:
                intact_list.remove(region)  # Remove that region from the final list
                break  # this will speed up the function by ending the 'for pseudo' loop if a match is successful
            else:
                pass

    intact_genes = Contig(regions=intact_list, name=contig.name, number=contig.number)

    return intact_genes


#TODO: FINISH THIS BOY
def write_intact_to_fasta(infile: str, outfile: str, contigs: List[Contig]) -> None:
    """Parses a multifasta file for regions and returns them in a list."""
    #
    # parsed = SeqIO.parse(infile, 'fasta')
    # print([fasta for fasta in parsed][0].description)
    # print([fasta for fasta in parsed][0].description)
    # exit()
    with open(infile, 'r') as infile, open(outfile, 'a') as output:
        lines = infile.readlines()
        for contig in contigs:
            region_index = 0
            for line_number, line in enumerate(lines):
                try:
                    if re.match("^>%s %s" % (regions[region_index].query, contig), line):
                        output.write("%s\n%s" % (line, lines[line_number+1]))
                        region_index += 1
                except IndexError:
                    pass


def write_pseudos_to_fasta(args, pseudofinder_regions: List[RegionInfo], outfile: str) -> None:
    """Parse genbank input file for the regions provided and write them to the output file in fasta format."""

    fasta_list = []

    for contig in SeqIO.parse(args.genome, "genbank"):
        try:
            pseudofinder_regions_on_contig = [region for region in pseudofinder_regions if region.contig == contig.name]
        except IndexError:
            continue

        coord_list = [(region.start, region.end) for region in pseudofinder_regions_on_contig]
        for counter, coordinate in enumerate(coord_list):
            fasta_list.append(SeqRecord(seq=contig.seq[coordinate[0]:coordinate[1]], id="%s_%04d" % (contig.name, counter + 1),
                                        description="%s-%s +" % (coordinate[0], coordinate[1])))

    SeqIO.write(fasta_list, open(outfile, "w"), "fasta")


def write_summary_file(args, outfile, file_dict) -> None:
    """Writes a summary file of statistics from the pseudo_finder run."""

    printable_args = common.convert_args_to_str(args)
    printable_stats = {k: str(v) for k, v in StatisticsDict.items()}

    with open(outfile, 'w') as logfile:

        logfile.write(
            "####### Summary from annotate/reannotate #######\n\n"
            "Date/time:\t" + current_time() + "\n\n"
            "#######    Files   #######\n"
            "Genome:\t" + printable_args.genome + "\n"
            "Database:\t" + printable_args.database + "\n"
            "Reference Genome:\t" + printable_args.reference + "\n"
            "BlastP:\t" + file_dict['blastp_filename'] + "\n"
            "BlastX:\t" + file_dict['blastx_filename'] + "\n"
            "Pseudogenes (GFF):\t" + file_dict['pseudos_gff'] + "\n"
            "Pseudogenes (Fasta):\t" + file_dict['pseudos_fasta'] + "\n"
            "Intact genes (GFF):\t" + file_dict['intact_gff'] + "\n"
            "Intact genes (protein seq):\t" + file_dict['intact_faa'] + "\n"
            "Intact genes (nucleotide seq):\t" + file_dict['intact_ffn'] + "\n"
            "Chromosome map:\t" + file_dict['chromosome_map'] + "\n\n"

            "#######  Settings  #######\n"
            "Distance:\t" + printable_args.distance + "\n"
            "hitcap:\t" + printable_args.hitcap + "\n"
            "Intergenic_length:\t" + printable_args.intergenic_length + "\n"
            "Intergenic_threshold:\t" + printable_args.intergenic_threshold + "\n"
            "Length_pseudo:\t" + printable_args.length_pseudo + "\n"
            "Shared_hits:\t" + printable_args.shared_hits + "\n"
            "contig_ends:\t" + printable_args.contig_ends + "\n"
            "max_dnds:\t" + printable_args.max_dnds + "\n"
            "max_ds:\t" + printable_args.max_ds + "\n"
            "min_ds:\t" + printable_args.min_ds + "\n\n"

            "####### Statistics #######\n"
            "#Input:\n"
            "Initial ORFs:\t" + printable_stats['ProteomeOrfs'] + "\n"
            "Number of contigs:\t" + printable_stats['NumberOfContigs'] + "\n"
            "#Output:\n"
            "Inital ORFs joined:\t" + printable_stats['FragmentedOrfs'] + "\n"
            "Pseudogenes (total):\t" + printable_stats['PseudogenesTotal'] + "\n"
            "Pseudogenes (too short):\t" + printable_stats['PseudogenesShort'] + "\n"
            "Pseudogenes (fragmented):\t" + printable_stats['PseudogenesFragmented'] + "\n"
            "Pseudogenes (no predicted ORF):\t" + printable_stats['PseudogenesIntergenic'] + "\n"
            "Pseudogenes (high dN/dS):\t" + printable_stats["dnds"] + "\n"
            "Intact genes:\t" + printable_stats['IntactORFs'] + "\n\n"

            "####### Output Key #######\n"
            "Initial ORFs joined:\t\tThe number of input open reading frames "
            "that have been merged and flagged as a fragmented pseudogene.\n"
            "Pseudogenes (too short):\tORFs smaller than the \"shared_hits\" cutoff.\n"
            "Pseudogenes (fragmented):\tPseudogenes composed of merging 2 or more input ORFs.\n"
            "Pseudogenes (high dN/dS):\tIncipient pseudogenes that look intact, but have an elevated dN/dS value compared to a reference gene.\n"
            "Intact genes:\t\t[Initial ORFs] - [Initial ORFs joined] - [Pseudogenes (too short)] - [Pseudogenes (high dN/dS)]\n"
        )


def reset_statistics_dict():
    """This function is needed because visualize.py was continuously accumulating values in StatisticsDict.
    Calling this function at the end of reannotate.py runs will ensure that the values are reset properly."""
    for key, value in StatisticsDict.items():
        if common.is_int(value):
            StatisticsDict[key] = 0

    # StatisticsDict['ProteomeOrfs'] = 0
    # StatisticsDict['NumberOfContigs'] = 0
    #
    # StatisticsDict['FragmentedOrfs'] = 0
    # StatisticsDict['PseudogenesTotal'] = 0
    # StatisticsDict['PseudogenesShort'] = 0
    # StatisticsDict['PseudogenesIntergenic'] = 0
    # StatisticsDict['PseudogenesFragmented'] = 0


def lastItem(ls):
    x = ''
    for i in ls:
        if i != "":
            x = i
    return x


def allButTheLast(iterable, delim):
    x = ''
    length = len(iterable.split(delim))
    for i in range(0, length - 1):
        x += iterable.split(delim)[i]
        x += delim
    return x[0:len(x) - 1]


def fastaReader(fasta_file):
    count = 0
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            if len(seq) > 0:
                Dict[header] = seq
                header = i[1:]
                # header = header.split(" ")[0]
                seq = ''
            else:
                header = i[1:]
                # header = header.split(" ")[0]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    return Dict


def fastaReader2(fasta_file):
    count = 0
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            if len(seq) > 0:
                Dict[header] = seq
                header = i[1:]
                header = header.split(" ")[0]
                seq = ''
            else:
                header = i[1:]
                header = header.split(" ")[0]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    return Dict


def integrate_dnds(func_gff: str, pseudo_gff: str, dnds_out: str, func_faa: str, func_ffn: str, proteome: str,
                   cds: str, pseudo_fasta: str, max_ds: float, min_ds: float, max_dnds: float):
    # args = common.get_args('dnds')
    funcDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    func = open(func_gff)  # reading original Intact gff outfile to dictionary
    for i in func:
        if not re.match(r'#', i):
            ls = i.rstrip().split("\t")
            locus = (lastItem(ls[8].split(";")))
            locus = locus.split("=")[1]
            funcDict[locus]["list"] = ls
            funcDict[locus]["locus"] = locus
        else:
            pass

    pseudoDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    pseudo = open(pseudo_gff)  # reading original pseudos gff outfile to dictionary
    for i in pseudo:
        if not re.match(r'#', i):
            ls = i.rstrip().split("\t")
            locus = (lastItem(ls[8].split(";")))
            locus = locus.split("=")[1]
            for j in locus.split(","):
                pseudoDict[j]["list"] = ls
                pseudoDict[j]["locus"] = locus
        else:
            pass

    # processing dnds output
    newPseudosDict = defaultdict(list)
    dndsfile = open(dnds_out)
    for i in dndsfile:
        ls = i.rstrip().split(",")
        locus = ls[0]
        if ls[2] != "dN":
            dn = float(ls[2])
            ds = float(ls[3])
            DNDS = float(ls[4])
            ref = ls[1]
            if ds < float(max_ds) and ds > float(min_ds):
                if DNDS > float(max_dnds):
                    if locus in funcDict.keys():
                        gffList = (funcDict[locus]["list"])
                        originalLocus = (funcDict[locus]["locus"])
                        commentsList = (gffList[8].split(";"))
                        comment = "Note=pseudogene candidate. Reason: dN/dS value between this gene and reference (" + ref + ") is " + str(
                            round(DNDS, 3)) + "."
                        newCommentsList = ([comment, commentsList[1], commentsList[2], commentsList[3]])
                        newCommentsLine = ";".join(newCommentsList)
                        gffLine = (("\t".join(gffList)))
                        newgffLine = (allButTheLast(gffLine, "\t") + "\t" + newCommentsLine + "\t" + str(DNDS))

                        newPseudosDict[originalLocus].append(newgffLine)
                    elif locus in pseudoDict.keys():

                        gffList = (pseudoDict[locus]["list"])
                        originalLocus = (pseudoDict[locus]["locus"])
                        commentsList = (gffList[8].split(";"))
                        newComment = (commentsList[
                                          0] + " Another Reason: dN/dS value between this gene and reference (" + ref + ") is " + str(
                            round(DNDS, 3)) + ".")
                        newCommentsList = [newComment, commentsList[1], commentsList[2], commentsList[3]]
                        newCommentsLine = ";".join(newCommentsList)
                        gffLine = (("\t".join(gffList)))
                        newgffLine = (allButTheLast(gffLine, "\t") + "\t" + newCommentsLine + "\t" + str(DNDS))

                        newPseudosDict[originalLocus].append(newgffLine)

    newPseudosDict2 = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in newPseudosDict.keys():
        dndsList = []
        for j in newPseudosDict[i]:
            dndsList.append(float(j.split("\t")[9]))
        ls = j.split("\t")
        DNDS = (statistics.mean(dndsList))
        commentsList = (ls[8].split(";"))
        comment = (commentsList[0])
        newComment = (comment.split(") is ")[0] + ") is " + str(round(DNDS, 3)))
        newCommentsList = [newComment, commentsList[1], commentsList[2], commentsList[3]]
        newCommentsLine = ";".join(newCommentsList)
        newgffLine = ls[0] + "\t" + ls[1] + "\t" + ls[2] + "\t" + ls[3] + "\t" + ls[4] + "\t" + ls[5] + "\t" + ls[
            6] + "\t" + ls[7] + "\t" + newCommentsLine
        newPseudosDict2[i] = newgffLine

    # re-writing Intact GFF output file
    func = open(func_gff)
    newFunctionalOut = func_gff.split(".")[0] + "-new.gff"
    out = open(newFunctionalOut, "w")
    for i in func:
        if not re.match(r'#', i):
            ls = i.rstrip().split("\t")
            locus = (lastItem(ls[8].split(";")))
            locus = locus.split("=")[1]
            if locus not in newPseudosDict2.keys():
                out.write(i.rstrip() + "\n")
        else:
            out.write(i.rstrip() + "\n")
    out.close()

    # re-writing pseudos GFF output file
    pseudo = open(pseudo_gff)
    newPseudoOut = pseudo_gff.split(".")[0] + "-new.gff"
    out = open(newPseudoOut, "w")
    for i in pseudo:
        if not re.match(r'#', i):
            ls = i.rstrip().split("\t")
            locus = (lastItem(ls[8].split(";")))
            locus = locus.split("=")[1]
            if locus in newPseudosDict2.keys():
                out.write(newPseudosDict2[locus] + "\n")
                newPseudosDict2.pop(locus, None)
            else:
                out.write(i.rstrip() + "\n")
        else:
            out.write(i.rstrip() + "\n")

    StatisticsDict["dnds"] = len(newPseudosDict2.keys())
    for i in newPseudosDict2.keys():
        out.write(newPseudosDict2[i] + "\n")

    out.close()

    os.system("mv %s %s" % (newPseudoOut, pseudo_gff))
    os.system("mv %s %s" % (newFunctionalOut, func_gff))

    # resorting the new GFF output file
    sortDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    pseudo = open(pseudo_gff)
    out = open(newPseudoOut, "w")
    for i in pseudo:
        if not re.match(r'#', i):
            ls = i.rstrip().split("\t")
            start = int(ls[3])
            sortDict[start] = i.rstrip()
        else:
            out.write(i.rstrip() + "\n")

    for i in sorted(sortDict.keys()):
        out.write(sortDict[i] + "\n")

    out.close()

    os.system("mv %s %s" % (newPseudoOut, pseudo_gff))

    # write intact genes to FASTA files
    faa = open(proteome)
    faa = fastaReader(faa)
    ffn = open(cds)
    ffn = fastaReader(ffn)
    ffn2 = open(cds)
    ffn2 = fastaReader2(ffn2)

    out = open(func_faa, "w")
    for i in faa.keys():
        locus = i.split(" ")[0]
        if locus in funcDict.keys():
            out.write(">" + i + "\n")
            out.write(faa[i] + "\n")
    out.close()

    out = open(func_ffn, "w")
    for i in ffn.keys():
        locus = i.split(" ")[0]
        if locus in funcDict.keys() and locus not in newPseudosDict2.keys():
            out.write(">" + i + "\n")
            out.write(ffn[i] + "\n")
    out.close()

    pseudos = open(pseudo_fasta)
    newPseudoSeqs = pseudo_fasta + "-new.fasta"
    out = open(newPseudoSeqs, "w")
    for i in pseudos:
        out.write(i.rstrip() + "\n")

    for i in newPseudosDict2.keys():
        ls = (newPseudosDict2[i].split("\t"))
        contig = ls[0].split("|")[2]
        header = contig + " " + ls[3] + "-" + ls[4] + " " + ls[6]
        seq = (ffn2[i])
        out.write(">" + header + "\n")
        out.write(seq + "\n")

    os.system("mv %s %s" % (newPseudoSeqs, pseudo_fasta))
    return len(newPseudosDict2.keys())


def write_all_outputs(args, genome, file_dict, visualize=False):
    """After analysis of pseudogenes is done, writes all necessary files."""

    if visualize:
        write_summary_file(args=args, outfile=file_dict['log'], file_dict=file_dict)

    else:
        intact = extract_features_from_genome(args, genome, 'CDS')
        pseudogenes = extract_features_from_genome(args, genome, 'pseudogene')
        write_gff(args, genome, file_dict['pseudos_gff'], 'pseudogene')
        write_gff(args, genome, file_dict['intact_gff'], 'CDS')
        write_fasta(intact, file_dict['intact_ffn'], 'nt')
        write_fasta(intact, file_dict['intact_faa'], 'aa')
        write_fasta(pseudogenes, file_dict['pseudos_fasta'], 'nt')
        write_gbk(args, genome, file_dict['gbk_out'])
        genome_map.full(genome=args.genome, gff=file_dict['pseudos_gff'], outfile=file_dict['chromosome_map'])
        write_summary_file(args=args, outfile=file_dict['log'], file_dict=file_dict)


def analysis_statistics(args, genome):
    """Adds statistics to the main dictionary about the run."""

    # cds
    intact = len(extract_features_from_genome(args, genome, 'CDS'))

    # pseudos
    pseudos = extract_features_from_genome(args, genome, 'pseudogene')
    total = len([x for x in pseudos if x.qualifiers.get('pseudo_type') != PseudoType.consumed])
    total = total + StatisticsDict["dnds"]
    short = len([x for x in pseudos if x.qualifiers.get('pseudo_type') == PseudoType.short])
    fragmented = len([x for x in pseudos if x.qualifiers.get('pseudo_type') == PseudoType.fragmented])
    intergenic = len([x for x in pseudos if x.qualifiers.get('pseudo_type') == PseudoType.intergenic])
    fragments = len(extract_features_from_genome(args, genome, 'consumed'))
    dnds = len([x for x in pseudos if x.qualifiers.get('pseudo_type') == PseudoType.dnds])

    StatisticsDict['PseudogenesTotal'] = total
    StatisticsDict['PseudogenesShort'] = short
    StatisticsDict['PseudogenesFragmented'] = fragmented
    StatisticsDict['PseudogenesIntergenic'] = intergenic
    StatisticsDict['dnds'] = dnds
    StatisticsDict['IntactORFs'] = intact
    StatisticsDict['FragmentedOrfs'] = fragments


def clear_analysis_qualifiers(args, genome):

    copy_to_clear = copy(genome)
    for seqrecord in copy_to_clear:
        for feature in seqrecord.features:
            if feature.type == 'intergenic':
                seqrecord.features.remove(feature)
            else:
                quals = feature.qualifiers
                quals.pop('nucleotide_seq', None)
                quals.pop('contig_id', None)
                quals.pop('hits', None)
                quals.pop('pseudo_type', None)
                quals.pop('note', None)
                quals.pop('dnds', None)
                quals.pop('parents', None)

    return copy_to_clear


def main():
    args = common.get_args('annotate')
    file_dict = common.file_dict(args)  # Declare filenames used throughout the rest of the program

    genome = gbk_to_seqrecord_list(args, args.genome)
    proteome = extract_features_from_genome(args, genome, 'CDS')
    write_fasta(seqs=proteome, outfile=file_dict['cds_filename'], seq_type='nt')
    write_fasta(seqs=proteome, outfile=file_dict['proteome_filename'], seq_type='aa')

    common.print_with_time("CDS extracted from:\t\t\t%s\n"
                           "\t\t\tWritten to file:\t\t\t%s." % (args.genome, file_dict['cds_filename']))

    intergenic = extract_features_from_genome(args, genome, 'intergenic')
    write_fasta(seqs=intergenic, outfile=file_dict['intergenic_filename'], seq_type='nt')

    common.print_with_time("Intergenic regions extracted from:\t%s\n"
                           "\t\t\tWritten to file:\t\t\t%s." % (args.genome, file_dict['intergenic_filename']))

    StatisticsDict['NumberOfContigs'] = len(genome)
    StatisticsDict['ProteomeOrfs'] = len(proteome)



    # get_CDSs(gbk=args.genome, out_fasta=file_dict['cds_filename'])
    # get_proteome(gbk=args.genome, out_faa=file_dict['proteome_filename'])
    # get_intergenic_regions(args=args, out_fasta=file_dict['intergenic_filename'])

    if args.reference:  # #########################################################################################
        common.print_with_time("Starting dN/dS analysis pipeline...")

        ref_genome = gbk_to_seqrecord_list(args, args.reference)
        ref_proteome = extract_features_from_genome(args, ref_genome, 'CDS')
        write_fasta(seqs=ref_proteome, outfile=file_dict['ref_cds_filename'], seq_type='nt')
        write_fasta(seqs=ref_proteome, outfile=file_dict['ref_proteome_filename'], seq_type='aa')
        #
        # get_CDSs(gbk=args.reference, out_fasta=file_dict['ref_cds_filename'])
        # get_proteome(gbk=args.reference, out_faa=file_dict['ref_proteome_filename'])
        dnds.full(args, file_dict)

        add_dnds_info_to_genome(args, genome, file_dict['dnds_out'])

    if args.diamond:  # run diamond
        if not args.skip_makedb:
            manage_diamond_db(args)

        run_diamond(args=args, search_type='blastp', in_fasta=file_dict['proteome_filename'],
                    out_tsv=file_dict['blastp_filename'])
        run_diamond(args=args, search_type='blastx', in_fasta=file_dict['intergenic_filename'],
                    out_tsv=file_dict['blastx_filename'])

    else:  # run vanilla blast
        if not args.skip_makedb:
            manage_blast_db(args)
        run_blast(args=args, search_type='blastp', in_fasta=file_dict['proteome_filename'],
                  out_tsv=file_dict['blastp_filename'])
        run_blast(args=args, search_type='blastx', in_fasta=file_dict['intergenic_filename'],
                  out_tsv=file_dict['blastx_filename'])

    add_blasthits_to_genome(args, genome, file_dict['blastp_filename'], 'blastp')
    add_blasthits_to_genome(args, genome, file_dict['blastx_filename'], 'blastx')
    update_intergenic_locations(args, genome)
    find_pseudos_on_genome(args, genome)
    # write_test_genome_output(file_dict, genome)


    # orfs = parse_blast(fasta_file=file_dict['proteome_filename'], blast_file=file_dict['blastp_filename'], blast_format='blastp')
    #
    # intergenic_regions = parse_blast(fasta_file=file_dict['intergenic_filename'], blast_file=file_dict['blastx_filename'], blast_format='blastx')
    # all_regions = orfs + intergenic_regions
    #
    # # Sorted list of contigs containing only orfs, no intergenic regions
    # orfs_by_contig = sort_contigs(loc=split_regions_into_contigs(lori=orfs))
    # # Sorted list of contigs containing orfs and intergenic regions
    # all_regions_by_contig = sort_contigs(loc=split_regions_into_contigs(lori=all_regions))
    #
    # pseudogenes = []
    # intact_genes = []
    #
    # for contig_index, contig in enumerate(all_regions_by_contig):
    #     print('\033[1m' + '%s\tChecking contig %s / %s for pseudogenes.\033[0m' % (current_time(),
    #                                                                                contig_index + 1,
    #                                                                                len(all_regions_by_contig))),
    #     sys.stdout.flush()
    #
    #     pseudos_on_contig = annotate_pseudos(args=args, contig=contig)  # Returns 'Contig' data type
    #     pseudogenes.extend(pseudos_on_contig.regions)  # List of regions
    #
    #     try:
    #         intact_genes_on_contig = get_intact_genes(contig=orfs_by_contig[contig_index],
    #                                                           pseudos=pseudos_on_contig.regions)
    #         intact_genes.extend(intact_genes_on_contig.regions)
    #     except IndexError:  # If there are no orfs on a small contig, an error will be thrown when checking that contig.
    #         continue
    #
    #     sys.stdout.flush()
    #
    #     # Write all output files
    #     write_genes_to_gff(args, lopg=pseudogenes, gff=file_dict['pseudos_gff'])
    #     write_genes_to_gff(args, lopg=intact_genes, gff=file_dict['intact_gff'])
    #     write_pseudos_to_fasta(args, pseudofinder_regions=pseudogenes, outfile=file_dict['pseudos_fasta'])
    #     # TODO: Activate this feature once you finish writing it
    #     # write_intact_to_fasta(infile=file_dict['proteome_filename'], outfile=file_dict['intact_faa'],
    #     #                           contigs=intact_genes)
    #     genome_map.full(genome=args.genome, gff=file_dict['pseudos_gff'], outfile=file_dict['chromosome_map'])
    #
    #
    # # INTEGRATING RESULTS FROM DNDS MODULE WITH THE REST OF THE PSEUDO-FINDER OUTPUT
    # if args.reference: #TODO: make sure integration still works
    #     newPseudos = integrate_dnds(func_gff=file_dict['intact_gff'], pseudo_gff=file_dict['pseudos_gff'],
    #                    dnds_out=file_dict['dnds_out'] + "/dnds-summary.csv", func_faa=file_dict['intact_faa'],
    #                    func_ffn=file_dict['intact_faa'], cds=file_dict['cds_filename'],
    #                    proteome=file_dict['proteome_filename'], pseudo_fasta=file_dict['pseudos_fasta'],
    #                    max_ds=args.max_ds, min_ds=args.min_ds, max_dnds=args.max_dnds)

        # for contig_index, contig in enumerate(all_regions_by_contig):
        #     print('\t\t\tNumber of ORFs on this contig: %s\n'
        #           '\t\t\tNumber of pseudogenes flagged: %s' % (
        #               len([region for region in contig.regions if region.region_type == RegionType.ORF]),
        #               len(pseudos_on_contig.regions) + newPseudos)),
        #     sys.stdout.flush()

        #write_summary_file(args=args, file_dict=file_dict)
    # else:
    #     for contig_index, contig in enumerate(all_regions_by_contig):
    #
    #         print('\t\t\tNumber of ORFs on this contig: %s\n'
    #               '\t\t\tNumber of pseudogenes flagged: %s' % (
    #                   len([region for region in contig.regions if region.region_type == RegionType.ORF]),
    #                   len(pseudos_on_contig.regions))),
    #         sys.stdout.flush()

        # write the log file
        #write_summary_file(args=args, file_dict=file_dict)

    common.print_with_time("Collecting run statistics and writing output files.")
    analysis_statistics(args, genome)
    write_all_outputs(args, genome, file_dict)



if __name__ == '__main__':
    main()









