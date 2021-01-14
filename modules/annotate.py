#!/usr/bin/env python3

from . import common
import re
import os
from copy import deepcopy
from enum import Enum
from typing import NamedTuple
from pandas.core.common import flatten
from Bio.Blast.Applications import NcbiblastpCommandline, NcbiblastxCommandline
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError

# This try block was added to stop a circular import error that occurs when this module is called from reannotate.py
try:
    from . import genome_map
except ImportError:
    pass

# importing selection module
try:
    from . import selection
except ImportError:
    pass

# importing interactive module
try:
    from . import interactive
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


PseudoType = Enum('PseudoType', ['short',
                                 'fragmented',
                                 'intergenic',
                                 'consumed',
                                 'dnds',
                                 'input'])

# Global dictionary, which will be called to write the log file
StatisticsDict = {
                    'BlastpFilename': '',
                    'BlastxFilename': '',
                    'NumberOfContigs': 0,
                    'ProteomeOrfs': 0,
                    'FragmentedOrfs': 0,
                    'PseudogenesInput': 0,
                    'PseudogenesTotal': 0,
                    'PseudogenesShort': 0,
                    'PseudogenesFragmented': 0,
                    'PseudogenesIntergenic': 0,
                    'OutputFiles': [],
                    'dnds': 0,
                    'IntactORFs': 0
                  }


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


def translate_cds(args, feature, seqrecord):
    """Translates a cds, with error handling and warnings."""
    if not feature.type == 'CDS':
        common.print_with_time('translate_cds() called on a feature that is not a CDS.')
        exit()

    try:  # Checks if the CDS feature contains a translation
        feature.qualifiers['translation'][0]
        return  # If translation already exists, no need to do anything else

    except KeyError:
        common.print_with_time(
            "WARNING: %s with locus tag \"%s\" has no translation qualifier in the input genbank file."
            " Pseudofinder will generate a translation using bacterial genetic code." % (feature.type, feature.qualifiers['locus_tag'][0]))

        try:
            feature.qualifiers['translation'] = [str(feature.translate(parent_sequence=seqrecord, table='Bacterial').seq)]
        except TranslationError as e:
            common.print_with_time("WARNING: Exception encountered when translating %s. Pseudofinder cannot analyze invalid CDS features.\n"
                                   "Error: %s\n"
                                   "Please fix input file. Pseudofinder will exit now." % (feature.qualifiers['locus_tag'][0], e))
            exit()


def add_qualifiers_to_features(args, seqrecord):
    """
    Adds additional qualifiers to each feature contained in the seqrecord, necessary for further analysis in the pipeline.
    """
    intergenic_counter = 1
    for feature in seqrecord.features:
        if feature.type == 'CDS' and feature.qualifiers.get('pseudo'):    # Finds CDS that are already annotated as pseudogenes
            feature.type = 'pseudogene'
            feature.qualifiers['pseudo_type'] = PseudoType.input
            feature.qualifiers['note'] = 'Annotated as pseudogene in input genbank file.'
            feature.qualifiers['parents'] = [feature.qualifiers['locus_tag'][0]]

        if feature.type == 'CDS':
            feature.qualifiers['dnds'] = []
            translate_cds(args, feature, seqrecord)

        if feature.type == 'intergenic':    # Changes for only intergenic regions
            feature.qualifiers['locus_tag'] = ['%s_ign_%05d' % (seqrecord.name, intergenic_counter)]
            intergenic_counter += 1

        if feature.type == 'CDS' or feature.type == 'intergenic' or feature.type == 'pseudogene':
            feature.qualifiers['nucleotide_seq'] = feature.extract(seqrecord.seq)
            feature.qualifiers['contig_id'] = seqrecord.name
            feature.qualifiers['hits'] = []

        if feature.type == 'CDS' or feature.type == 'intergenic':
            feature.qualifiers['pseudo_type'] = None
            feature.qualifiers['note'] = ''
            feature.qualifiers['parents'] = []


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


def extract_features_from_genome(args, genome, feature_type: list or str) -> list:
    """
    Returns all features from genome as a list.
    feature_type can be a string or infinitely deeply nested lists of strings.
    ie, 'CDS' or ['CDS', 'intergenic'] work but something crazy like
     ['CDS', ['intergenic', 'ncRNA', ['pseudogene']]] would also be valid.
    Will throw an error if feature_type is anything other than a list or a string.
    """
    feature_list = []
    if isinstance(feature_type, list):
        for item in feature_type:
            feature_list.extend(extract_features_from_genome(args, genome, item))

    elif isinstance(feature_type, str):
        for record in genome:
            feature_list.extend(parse_features_from_record(record, feature_type))

    else:
        raise TypeError("extract_features_from_genome() received invalid type: %s" % type(feature_type))

    # remove duplicates
    feature_list = list(set(feature_list))

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

    if feature.type == 'pseudogene':
        list.append("old_locus_tag=" + ",".join(feature.qualifiers['parents']))

    return ";".join(list)


def write_gff(args, genome, outfile: str, seq_type):
    """Writes a GFF3 format compliant GFF file with the desired sequence type."""

    with open(outfile, 'w') as outfile:
        # Header
        outfile.write("##gff-version 3\n#!annotation-date\t%s\n" % (common.current_time()))
        for seqrecord in genome:
            entry = ["##sequence-region",
                     seqrecord.name,
                     #"gnl|Prokka|%s" % seqrecord.name,
                     str(1),
                     str(len(seqrecord))]
            outfile.write(" ".join(entry) + '\n')

        seqs = extract_features_from_genome(args, genome, seq_type)
        for seq in seqs:
            seq_info = {'seq_id': seq.qualifiers['contig_id'],
                        #'seq_id': "gnl|Prokka|%s" % seq.qualifiers['contig_id'],
                        'source': "pseudofinder",
                        'feature_type': "gene", # TODO: change this to 'CDS' or 'pseudogene' if given the go-ahead to do so
                        'feature_start': str(seq.location.start),
                        'feature_end': str(seq.location.end),
                        'score': '.',
                        'strand': format_strand(seq.strand),
                        'phase': '.',   # TODO: Understand what phase is and how to properly declare this
                        'attributes': attribute_string(seq)}

            outfile.write(gff_entry(**seq_info) + "\n")


# TODO: Works, but not completely. Some features are missing info such as notes for reasons unclear to me,
# TODO: and some features have extra info that should have been cleared but isn't.
def write_gbk(args, genome, outfile):
    """Writes a genbank file, excluding any analysis qualifiers present in the data structure."""

    cleared_genome = clear_analysis_qualifiers(args, genome)

    with open(outfile, "w") as output:
        SeqIO.write(cleared_genome, output, "genbank")


def run_blast(args, search_type: str, in_fasta: str, out_tsv: str) -> None:
    """"Run BLASTP or BLASTX with fasta file against DB of your choice."""

    common.print_with_time("%s executed with %s threads on %s." % (search_type, args.threads, in_fasta))

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

        common.print_with_time("Diamond %s executed with %s threads on %s." % (search_type, args.threads, in_fasta))

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
    relevant_features = [f for f in features_with_locus_tags(genome) if f.type == 'CDS' or f.type == 'intergenic' or 'pseudo' in f.type.lower()]
    feature_dict = {feature.qualifiers['locus_tag'][0]: feature for feature in relevant_features}

    for hit in hit_list:
        try:
            feature_dict[hit.query].qualifiers['hits'].append(hit)
        except KeyError:
            common.print_with_time("Potential duplicate locus tags detected in genome, unable to process. "
                                   "Offending tag: %s "
                                   "Exiting now." % hit.query)
            exit()


def convert_csv_to_dnds(dnds_file):
    """
    Reads values from the dNdS csv file and converts each row into a dnds_data entry.
    Returns a list of dnds_data
    """
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

    for feature in seqrecord.features:
        # Option 0: This will catch pseudogenes that were annotated as such in the input genbank file.
        if feature.type == 'pseudogene':
            pass

        # Option 1: DNDS pseudogene
        if feature.type == 'CDS' and len(feature.qualifiers['dnds']) > 0:
            dnds_val = feature.qualifiers['dnds'][0].dnds

            if dnds_val >= args.max_dnds:
                feature.type = 'pseudogene'
                feature.qualifiers['pseudo_type'] = PseudoType.dnds
                feature.qualifiers['note'] = 'Pseudogene candidate. Reason: Elevated dN/dS (%s).' % round(dnds_val, 3)
                feature.qualifiers['parents'].append(feature.qualifiers['locus_tag'][0])

        # Option 2: Short pseudogene
        if feature.type == 'CDS' and len(feature.qualifiers['hits']) > 0:
            if feature_length_relative_to_hits(feature) <= args.length_pseudo:
                feature.type = 'pseudogene'
                feature.qualifiers['pseudo_type'] = PseudoType.short
                feature.qualifiers['note'] = 'Pseudogene candidate. Reason: ORF is %s%% of the average length of ' \
                                             'hits to this gene.' % (round(feature_length_relative_to_hits(feature)*100, 1))
                feature.qualifiers['parents'].append(feature.qualifiers['locus_tag'][0])

        # Option 3: Intergenic pseudogene
        if feature.type == 'intergenic':
            if len(feature.qualifiers['hits']) / args.hitcap >= args.intergenic_threshold:
                feature.type = 'pseudogene'
                feature.qualifiers['pseudo_type'] = PseudoType.intergenic
                feature.qualifiers['note'] = 'Pseudogene candidate. Reason: Intergenic region with %s blast hits.' % len(feature.qualifiers['hits'])
                feature.qualifiers['parents'].append(feature.qualifiers['locus_tag'][0])
        else:
            pass


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
        # TODO: This occurs if features from + and - strand are going to be combined. Very rarely does this happen
        # TODO: but there could be biologically relevant reasons for this to occur ie two fragments on + strand
        # TODO: separated by an ORF on the - strand which is an insertion sequence.
        # TODO: We should explore options to handle these cases in the future.
        tags = [feature.qualifiers['locus_tag'][0] for feature in fragments]
        parent_tags = [feature.qualifiers.get('parents') for feature in fragments]

        # TODO: understand why 'tags + parent_tags' results in a nested list and fix the source of the issue rather than forcibly flattening.
        # Order of operations: flatten -> remove duplicates -> remove 'None' -> convert to list
        all_tags = list(filter(None, set(flatten(tags + parent_tags))))

        common.print_with_time("WARNING: Pseudogene detected which traverses features on (+) and (-) strands.\n"
                               "We recommend manual inspection of this region.\n"
                               "Features involved: %s" % all_tags)

        strand = 0
        # raise RuntimeError("Trying to combine genes on opposite strands.")

    hits = []
    for feature in fragments:
        hits.extend(feature.qualifiers['hits'])

    parents = []
    for fragment in fragments:
        frag_parents = fragment.qualifiers.get('parents')
        if frag_parents is not None:
            parents.extend(frag_parents)
        parents.extend([feature.qualifiers['locus_tag'][0] for feature in fragments])
        parents = list(set(flatten(parents)))

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
        while i < len(features) - 2:
            f1 = features[i]
            f2 = features[i+1]
            f3 = features[i+2]

            if adjacent_fragments_match(args, f1, f2):
                create_fragmented_pseudo(args, (f1, f2), seqrecord)
                manage_parent_fragments((f1, f2))
                if i > 0:
                    i -= 1

            elif adjacent_fragments_match(args, f1, f3):
                create_fragmented_pseudo(args, (f1, f2, f3), seqrecord)
                manage_parent_fragments((f1, f2, f3))
                if i > 0:
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
                               "Number of pseudogenes flagged: %s" % (num_ORFs, num_pseudos))

    update_locus_tags(args, genome)


def write_summary_file(args, outfile, file_dict) -> None:
    """Writes a summary file of statistics from the pseudo_finder run."""

    printable_args = common.convert_args_to_str(args)
    printable_stats = {k: str(v) for k, v in StatisticsDict.items()}

    with open(outfile, 'w') as logfile:

        logfile.write(
            "####### Summary from annotate/reannotate #######\n\n"
            "Date/time:\t" + common.current_time() + "\n\n"
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
            "Initial pseudogenes:\t" + printable_stats['PseudogenesInput'] + "\n"
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


def write_all_outputs(args, genome, file_dict, visualize=False):
    """After analysis of pseudogenes is done, writes all necessary files."""

    if visualize:   # Only write the necessary summary file for visualize module to speed up analysis
        write_summary_file(args=args, outfile=file_dict['log'], file_dict=file_dict)

    else:
        intact = extract_features_from_genome(args, genome, 'CDS')
        pseudogenes = extract_features_from_genome(args, genome, 'pseudogene')
        write_gff(args, genome, file_dict['pseudos_gff'], 'pseudogene')
        write_gff(args, genome, file_dict['intact_gff'], 'CDS')
        write_fasta(intact, file_dict['intact_ffn'], 'nt')
        write_fasta(intact, file_dict['intact_faa'], 'aa')
        write_fasta(pseudogenes, file_dict['pseudos_fasta'], 'nt')
        # write_gbk(args, genome, file_dict['gbk_out']) #TODO: understand why this doesn't work as it should. See notes inside function
        # common.write_test_genome_output(file_dict, genome)
        # interactive.genome_to_graphs(args=args, file_dict=file_dict, genome=genome)
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


 # TODO: This is deleting almost all, but not all things that we don't want. I am not sure why some features are slipping past this
 # TODO: since they show up when you call seqrecord.features, but then do not get acted on.
 # TODO: for example you can find some features in the output genbank that have 'nucleotide_seq' qualifiers. No idea why.
def clear_analysis_qualifiers(args, genome):
    """Clears all analysis qualifiers that we don't want to write to a properly formatted genbank file."""
    copy_to_clear = deepcopy(genome)
    for seqrecord in copy_to_clear:
        for feature in seqrecord.features:
            quals = feature.qualifiers
            quals.pop('nucleotide_seq', None)
            quals.pop('contig_id', None)
            quals.pop('hits', None)
            quals.pop('pseudo_type', None)
            quals.pop('note', None)
            quals.pop('dnds', None)
            quals.pop('parents', None)

            if feature.type == 'intergenic' or feature.type == 'consumed':
                seqrecord.features.remove(feature)

    return copy_to_clear


def main():
    args = common.get_args('annotate')
    file_dict = common.file_dict(args)  # Declare filenames used throughout the rest of the program

    genome = gbk_to_seqrecord_list(args, args.genome)

    proteome = extract_features_from_genome(args, genome, 'CDS')
    write_fasta(seqs=proteome, outfile=file_dict['cds_filename'], seq_type='nt')
    write_fasta(seqs=proteome, outfile=file_dict['proteome_filename'], seq_type='aa')
    common.print_with_time("CDS extracted from:\t\t\t%s\n"
                           "Written to file:\t\t\t%s." % (args.genome, file_dict['cds_filename']))

    intergenic = extract_features_from_genome(args, genome, 'intergenic')
    write_fasta(seqs=intergenic, outfile=file_dict['intergenic_filename'], seq_type='nt')
    common.print_with_time("Intergenic regions extracted from:\t%s\n"
                           "Written to file:\t\t\t%s." % (args.genome, file_dict['intergenic_filename']))

    input_pseudos = extract_features_from_genome(args, genome, 'pseudogene')
    if len(input_pseudos) > 0:
        write_fasta(seqs=input_pseudos, outfile=file_dict['input_pseudos_filename'], seq_type='nt')
        common.print_with_time("%s pseudogenes found in genbank file:\t%s\n"
                               "Written to file:\t\t\t%s." % (len(input_pseudos), args.genome, file_dict['input_pseudos_filename']))

    StatisticsDict['NumberOfContigs'] = len(genome)
    StatisticsDict['ProteomeOrfs'] = len(proteome)
    StatisticsDict['PseudogenesInput'] = len(input_pseudos)

    if args.reference:  # #########################################################################################
        common.print_with_time("Starting dN/dS analysis pipeline...")
        ref_genome = gbk_to_seqrecord_list(args, args.reference)
        ref_proteome = extract_features_from_genome(args, ref_genome, 'CDS')
        write_fasta(seqs=ref_proteome, outfile=file_dict['ref_cds_filename'], seq_type='nt')
        write_fasta(seqs=ref_proteome, outfile=file_dict['ref_proteome_filename'], seq_type='aa')
        selection.full(args, file_dict)
        add_dnds_info_to_genome(args, genome, file_dict['dnds_out'])

    if args.diamond:  # run diamond
        if not args.skip_makedb:
            manage_diamond_db(args)

        run_diamond(args=args, search_type='blastp', in_fasta=file_dict['proteome_filename'],
                    out_tsv=file_dict['blastp_filename'])
        run_diamond(args=args, search_type='blastx', in_fasta=file_dict['intergenic_filename'],
                    out_tsv=file_dict['blastx_filename'])
        if len(input_pseudos) > 0:
            run_diamond(args=args, search_type='blastx', in_fasta=file_dict['input_pseudos_filename'],
                        out_tsv=file_dict['blastx_pseudos_filename'])

    else:  # run vanilla blast
        if not args.skip_makedb:
            manage_blast_db(args)
        run_blast(args=args, search_type='blastp', in_fasta=file_dict['proteome_filename'],
                  out_tsv=file_dict['blastp_filename'])
        run_blast(args=args, search_type='blastx', in_fasta=file_dict['intergenic_filename'],
                  out_tsv=file_dict['blastx_filename'])
        if len(input_pseudos) > 0:
            run_blast(args=args, search_type='blastx', in_fasta=file_dict['input_pseudos_filename'],
                        out_tsv=file_dict['blastx_pseudos_filename'])

    add_blasthits_to_genome(args, genome, file_dict['blastp_filename'], 'blastp')
    add_blasthits_to_genome(args, genome, file_dict['blastx_filename'], 'blastx')
    if len(input_pseudos) > 0:
        add_blasthits_to_genome(args, genome, file_dict['blastx_pseudos_filename'], 'blastx')

    update_intergenic_locations(args, genome)
    find_pseudos_on_genome(args, genome)
    common.print_with_time("Collecting run statistics and writing output files.")
    analysis_statistics(args, genome)
    write_all_outputs(args, genome, file_dict)


if __name__ == '__main__':
    main()









