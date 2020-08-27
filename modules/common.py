#!/usr/bin/env python3

import argparse
import os
import sys
from time import localtime, strftime


def bold(x):
    start_bold = '\033[1m'
    end_bold = '\033[0m'
    return start_bold + x + end_bold


def current_time() -> str:
    """Returns the current time when this function was executed."""
    return str(strftime("%Y-%m-%d %H:%M:%S", localtime()))


def print_with_time(x):
    print(bold(current_time()) + '\t' + x)
    sys.stdout.flush()


def unpack_arg(arg_group, arg: dict):
    name = (arg['short'], arg['long'])
    del arg['short']
    del arg['long']
    arg_group.add_argument(*name, **arg)


def get_args(module='None', **kwargs):
    names_only = kwargs.get('names_only', False)
    genome = {
        'short': '-g',
        'long': '--genome',
        'help': 'Please provide your genome file in the genbank format.',
        'required': True
    }
    database = {
        'short': '-db',
        'long': '--database',
        'help': 'Please provide name (if $BLASTB is set on your system) or absolute path of your blast database.',
        'required': True
    }
    outprefix = {
        'short': '-op',
        'long': '--outprefix',
        'help': 'Specify an output prefix.',
        'required': True
    }
    threads = {
        'short': '-t',
        'long': '--threads',
        'help': 'Please provide total number of threads to use for blast, default is %(default)s.',
        'required': False,
        'default': 4
    }
    length_pseudo = {
        'short': '-l',
        'long': '--length_pseudo',
        'help': 'Please provide percentage of length for pseudo candidates, default is %(default)s. '
                '\nExample: \"-l 0.50\" will consider genes that are less than 50%% of the average '
                'length of similar genes.',
        'required': False,
        'default': 0.65,
        'type': float
    }
    intergenic_length = {
        'short': '-i',
        'long': '--intergenic_length',
        'help': 'Please provide length of intergenic regions to check, default is %(default)s bp.',
        'required': False,
        'default': 30,
        'type': int
    }
    shared_hits = {
        'short': '-s',
        'long': '--shared_hits',
        'help': 'Percentage of blast hits that must be shared in order to join two nearby regions,'
                ' default is %(default)s. \nExample: \"-s 0.50\" will merge nearby regions if '
                'they shared 50%% of their blast hits.',
        'required': False,
        'default': 0.50,
        'type': float
    }
    evalue = {
        'short': '-e',
        'long': '--evalue',
        'help': 'Please provide e-value for blast searches. Default is %(default)s.',
        'required': False,
        'default': '1e-4',
        'type': str
    }
    distance = {
        'short': '-d',
        'long': '--distance',
        'help': 'Maximum distance between two regions to consider joining them. Default is %(default)s.',
        'required': False,
        'default': 1000,
        'type': int
    }
    hitcap = {
        'short': '-hc',
        'long': '--hitcap',
        'help': 'Maximum number of allowed hits for BLAST. Default is %(default)s.',
        'required': False,
        'default': 15,
        'type': int
    }
    contig_ends = {
        'short': '-ce',
        'long': '--contig_ends',
        'help': 'Forces the program to include intergenic regions at contig ends. If not specified,\n the '
                'program will ignore any sequence before the first ORF and after the last ORF on a contig.',
        'required': False,
        'default': False,
        'action': 'store_true'
    }
    intergenic_threshold = {
        'short': '-it',
        'long': '--intergenic_threshold',
        'help': 'Number of BlastX hits needed to annotate an intergenic region as a pseudogene.\n'
                'Calculated as a percentage of maximum number of allowed hits (--hitcap).\n'
                'Default is %(default)s.',
        'required': False,
        'default': 0.30,
        'type': float
    }
    reference = {
        'short': '-ref',
        'long': '--reference',
        'help': 'Please provide a reference genome if you would like for the program to carry out\n'
                'maximum-likelihood phylogenentic analysis using PAML, and calculate dN/dS values for each \n'
                'identified ORF in your query genome.',
        'required': False,
        'type': str
    }
    max_dnds = {
        'short': '-dnds',
        'long': '--max_dnds',
        'help': 'maximum dN/dS value for gene too be considered \'intact\'. Default is %(default)s.',
        'required': False,
        'default': 0.30,
        'type': float
    }
    max_ds = {
        'short': '-M',
        'long': '--max_ds',
        'help': 'maximum dS value for dN/dS calculation. Default is %(default)s.',
        'required': False,
        'default': 3.00,
        'type': float
    }
    min_ds = {
        'short': '-m',
        'long': '--min_ds',
        'help': 'minimum dS value for dN/dS calculation. Default is %(default)s.',
        'required': False,
        'default': 0.001,
        'type': float
    }
    diamond = {
        'short': '-di',
        'long': '--diamond',
        'help': 'Use DIAMOND BLAST as the search engine. If not specified,\n standard BLAST will be used.',
        'required': False,
        'default': False,
        'action': 'store_true'
    }
    skip = {
        'short': '-sk',
        'long': '--skip',
        'help': 'If you provided a reference genome for dN/dS analysis, already ran pseudo-finder, and \n'
                'would just like to re-run with one one of the flags for --max_ds, --max_dnds, or --min_ds altered. \n'
                'This flag allows you to skip the time-consuming steps (which you don\'t need if you already \n'
                'have the codeml output) and use the previously-created output.',
        'required': False,
        'default': False,
        'action': 'store_true'
    }
    blastp = {
        'short': '-p',
        'long': '--blastp',
        'help': 'Specify an input blastp file.',
        'required': True,
        'type': str
    }
    blastx = {
        'short': '-x',
        'long': '--blastx',
        'help': 'Specify an input blastx file.',
        'required': True,
        'type': str
    }
    logfile = {
        'short': '-log',
        'long': '--logfile',
        'help': 'Provide the log file from the run that generated the blast files.',
        'required': True,
        'type': str
    }
    gff = {
        'short': '-gff',
        'long': '--gff',
        'help': 'Provide your pseudogene calls in GFF format.',
        'required': True,
        'type': str
    }
    resolution = {
        'short': '-r',
        'long': '--resolution',
        'help': 'Specifies the resolution of your 3D plot. '
                'Lowering the resolution will increase the speed of this process.'
                '\'-r 100\' means the program\n will generate a 100x100 data matrix. '
                'Default is %(default)s.',
        'required': False,
        'type': int,
        'default': 20
    }
    keep_files = {
        'short': '-k',
        'long': '--keep_files',
        'help': 'Specifies whether to keep all output files. Default is %(default)s, which will remove'
                ' all files except for the graph and data matrix after running.',
        'required': False,
        'default': False,
        'action': 'store_true'
    }
    title = {
        'short': '-ti',
        'long': '--title',
        'help': 'Specifies a title for your plot. Default is None.',
        'required': False,
        'default': None,
        'type': str
    }
    ###############################################################################
    search_engine = {
        'short': '-s',
        'long': '--search_engine',
        'help': 'search engine to use (blast/diamond). Default = blast.',
        'required': False,
        'default': "blast",
        'type': str
    }
    outdir = {
        'short': '-out',
        'long': '--outdir',
        'help': 'name output directory.',
        'required': False,
        'default': "dnds_out",
        'type': str
    }
    control_file = {
        'short': '-ctl',
        'long': '--control_file',
        'help': 'template control file for codeml.',
        'required': False,
        'default': None,
        'type': str
    }
    ref_contigs = {
        'short': '-r',
        'long': '--ref_contigs',
        'help': 'Reference contigs (provide this if you do not have ORFs for the reference genome).',
        'required': False,
        'default': None,
        'type': str
    }
    ref_prots = {
        'short': '-ra',
        'long': '--ref_prots',
        'help': 'Reference ORFs in amino acid format.',
        'required': True,
        'default': None,
        'type': str
    }
    ref_genes = {
        'short': '-rn',
        'long': '--ref_genes',
        'help': 'Reference ORFs in nucleic acid format.',
        'required': True,
        'default': None,
        'type': str
    }
    prots = {
        'short': '-a',
        'long': '--prots',
        'help': 'ORFs in amino acid format.',
        'required': True,
        'default': None,
        'type': str
    }
    genes = {
        'short': '-n',
        'long': '--genes',
        'help': 'ORFs in nucleic acid format.',
        'required': True,
        'default': None,
        'type': str
    }

    if names_only:
        to_remove = ['module', 'kwargs', 'names_only', 'to_remove']
        names = list(locals().keys())
        for i in to_remove:
            names.remove(i)
        return names

    if module == 'annotate':
        usage = '[pseudofinder.py annotate -g GENOME -db DATABASE -op OUTPREFIX]'
        required_args = [genome, database, outprefix]
        optional_args = [threads, intergenic_length, length_pseudo, shared_hits, evalue, distance, hitcap,
                         contig_ends, intergenic_threshold, reference, max_dnds, max_ds, min_ds, diamond, skip]

    elif module == 'reannotate':
        usage = '[pseudofinder.py reannotate -g GENOME -p BLASTP -x BLASTX -hc HITCAP -op OUTPREFIX]'
        required_args = [genome, blastp, blastx, logfile, outprefix]
        optional_args = [length_pseudo, shared_hits, intergenic_threshold, distance, max_dnds, max_ds, min_ds, reference]

    elif module == 'dnds':
        usage = 'pseudofinder.py dnds -a PROTEIN_SEQS -n GENE_SEQS -ra REFERENCE_PROTEINS -rn REFERENCE_GENES -out OUTPUT_DIR'
        required_args = [prots, genes, ref_prots, ref_genes]
        optional_args = [control_file, outdir, search_engine, min_ds, max_ds, max_dnds, threads, ref_contigs]

    elif module == 'genome_map':
        usage = '[pseudofinder.py map -g GENOME -gff GFF -op OUTPREFIX]'
        required_args = [genome, gff, outprefix]
        optional_args = []

    elif module == 'visualize':
        usage = '[pseudofinder.py visualize -g GENOME -op OUTPREFIX -p BLASTP -x BLASTX -log LOGFILE]'
        required_args = [genome, outprefix, blastp, blastx, logfile]
        optional_args = [distance, intergenic_threshold, resolution, keep_files, title]

    elif module == 'test':
        usage = '[pseudofinder.py test -db DATABASE]'
        required_args = [database]
        optional_args = [threads]

    else:
        print("Module not found. Please check your get_args() function call.")
        exit()

    usage = usage + " or [pseudofinder.py %s --help] for more options." % module
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, usage=bold(usage))
    always_required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Adjustable parameters')

    for arg in required_args:
        unpack_arg(always_required, arg)

    for arg in optional_args:
        unpack_arg(optional, arg)

    # parse_known_args will create a tuple of known arguments in the first position and unknown in the second.
    # We only care about the known arguments, so we take [0].
    args = parser.parse_known_args()[0]
    return args


def reconcile_args(command_line_args, logged_args):

    args = command_line_args  # create a namespace with command line args taking priority
    # if argument does not exist or is None, it will be replaced by whatever is found in the log
    for k, v in logged_args.items():
        try:
            if getattr(args, k) is None:
                setattr(args, k, v)
        except AttributeError:
            setattr(args, k, v)
    return args


def file_dict(args, **kwargs):
    outprefix = kwargs.get('outprefix', None)

    if outprefix:
        base_outfile_name = outprefix + "_"
    else:
        base_outfile_name = args.outprefix + "_"

    file_dict = {
        'cds_filename': base_outfile_name + "cds.fasta",
        'ref_cds_filename': base_outfile_name + "ref_cds.fasta",
        'proteome_filename': base_outfile_name + "proteome.faa",
        'ref_proteome_filename': base_outfile_name + "ref_proteome.faa",
        'intergenic_filename': base_outfile_name + "intergenic.fasta",
        'blastp_filename': base_outfile_name + "proteome.faa" + ".blastP_output.tsv",
        'blastx_filename': base_outfile_name + "intergenic.fasta" + ".blastX_output.tsv",
        'pseudos_gff': base_outfile_name + "pseudos.gff",
        'pseudos_fasta': base_outfile_name + "pseudos.fasta",
        'functional_gff': base_outfile_name + "functional.gff",
        'functional_faa': base_outfile_name + "functional.faa",
        'functional_ffn': base_outfile_name + "functional.ffn",
        'chromosome_map': base_outfile_name + "map.pdf",
        'dnds_out': base_outfile_name + "dnds",
        'log': base_outfile_name + "log.txt",
        'ctl': os.path.dirname(os.path.dirname(__file__)) + "/codeml-2.ctl"
    }
    return file_dict


def is_int(x: str):
    try:
        a = float(x)
        b = int(x)
    except ValueError:
        return False
    else:
        return a == b


def is_float(x: str):
    try:
        float(x)
    except ValueError:
        return False
    else:
        return True


def literal_eval(x: str):
    """
    My version of ast.literal_eval() since that function is more complicated than it needs to be for our use case.

    Basically just takes a string and tries to evaluate it to the correct data type in the following hierarchy:
        Bool
        Int
        Float
        Str

    """
    if x == 'True':
        return True
    elif x == 'False':
        return False
    elif is_int(x):
        return int(x)
    elif is_float(x):
        return float(x)
    else:
        return x
#
# print("TRUE TESTS")
# print(literal_eval('True') == True)
# print(literal_eval('False') == False)
# print(literal_eval('0.50') == 0.50)
# print(literal_eval('3') == 3)
#
# print("FALSE TESTS")
# print(literal_eval('4.57') == 0.50)
# print(literal_eval('Truer') == True)
