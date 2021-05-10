#!/usr/bin/env python3

import argparse
from argparse import Namespace
import os
import sys
import re
import warnings
from time import localtime, strftime
from contextlib import contextmanager
from Bio import SeqIO


def bold(x):
    start_bold = '\033[1m'
    end_bold = '\033[0m'
    return start_bold + x + end_bold


def current_time() -> str:
    """Returns the current time when this function was executed."""
    return str(strftime("%Y-%m-%d %H:%M:%S", localtime()))


def print_with_time(x):
    # This will auto indent newlines to accomodate the datetime
    x_indented = '\t' + x.replace('\n', '\n\t\t\t')
    print(bold(current_time()) + x_indented)
    sys.stdout.flush()


def is_int(x: str):
    try:
        a = float(x)
        b = int(x)
    except (ValueError, TypeError):
        return False
    else:
        return a == b


def is_float(x: str):

    if x == 'None' or x is None:
        return False
    elif x == 'True' or x is True:
        return False
    elif x == 'False' or x is False:
        return False
    else:
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
        Nonetype
        Bool
        Int
        Float
        Str
    """

    if x == 'None':
        return None
    elif x == 'True':
        return True
    elif x == 'False':
        return False
    elif is_int(x):
        return int(x)
    elif is_float(x):
        return float(x)
    else:
        return x


@contextmanager
def suppress_output_to_console():
    """Prevents writing to stdout. Use in the following way:

    with suppress_output_to_console():
        action_to_silence()

    outside_of_scope_will_print_again()
    """
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout


def unpack_arg(arg_group, arg: dict, deprecated=False):
    name = (arg['short'], arg['long'])
    del arg['short']
    del arg['long']
    if deprecated:
        arg['help'] = argparse.SUPPRESS
        arg['required'] = False
        del arg['default']
    arg_group.add_argument(*name, **arg)


def arg_string(arg):
    """returns a string for usage message, ie for genome returns -g GENOME."""
    return arg['short'] + ' ' + arg['long'][2:].upper()


def usage_message(module, required_args: list):
    message = '[pseudofinder.py ' + module + ' ' + ' '.join([arg_string(x) for x in required_args]) + ']'
    message = message + " or [pseudofinder.py %s --help] for more options." % module
    return message


def verify_gbk(gbk):
    """
    Ensures genbank files is Genbank/ENA/DDJB compliant.
    """
    compliant = False
    for contig in SeqIO.parse(gbk, "genbank"):
        for feature in contig.features:
            if feature.type == "gene":
                compliant = True
                break

    if compliant:
        return True
    else:
        raise RuntimeError('pseudofinder has detected your GENOME file is not Genbank compliant. '
                           'If generated with Prokka, please ensure the \'--compliant\' flag is used.')


def verify_numeric(args):
    for key, value in vars(args).items():
        if is_float(value):
            if float(value) > 0:
                return True
            else:
                raise RuntimeError('pseudofinder has detected your %s argument is less than or equal to 0 (%s = %s). Please enter a positive value.' % (key, key, value))


# TODO: expand on this if we need to add more deprecated arguments
def verify_deprecated(args, deprecated_args):
    """This is where a real deprecated argument function should go, however right now we will just manually check
    args.distance since it is the only one."""
    if args.distance:
        warnings.warn("--distance parameter is deprecated. Fragment merging algorithm logic improved such that it is no longer needed."
                      "Program will continute to run but please be advised that changing the --distance parameter now has no effect.", Warning)


def verify_args(args, deprecated_args):
    """
    Ensures that the arguments provided are of the correct format beyond data type.
    """
    try:
        verify_deprecated(args, deprecated_args)
        verify_gbk(args.genome)
        verify_numeric(args)
    except AttributeError:
        pass


def get_args(module='None', **kwargs):
    names_only = kwargs.get('names_only', False)
    genome = {
        'short': '-g',
        'long': '--genome',
        'help': 'Please provide your genome file in the genbank format.',
        'required': False if module == 'test' else True
    }
    annotated_genome = {
        'short': '-ag',
        'long': '--annotated_genome',
        'help': 'Testing feature. Genbank file with all extra qualifiers and pseudogene annotations contained.\n'
                'If this is the only argument supplied, the figure will appear in your browser.\n'
                'If you supply an outprefix, figure will be saved as an html file.',
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
        'required': False if module == 'interactive' else True
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

    perc_id = {
        'short': '-id',
        'long': '--perc_id',
        'help': 'the candidate gene/pseudogene needs to be at least this proportion of '
                'reference gene length for frameshift and dN/dS analysis (Sleuth module or Annotate module with -reference). Default is %(default)s.',
        'required': False,
        'default': '0.25',
        'type': str
    }

    perc_cov = {
        'short': '-cov',
        'long': '--perc_cov',
        'help': 'the candidate gene/pseudogene needs to be at least this percent identicial '
                'to reference gene for frameshift and dN/dS analysis (Sleuth module or Annotate module with -reference). Default is %(default)s.',
        'required': False,
        'default': '0.25',
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
        'default': None,
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
    dnds_out = {
        'short': '-do',
        'long': '--dnds_out',
        'help': 'dnds output from previous run. Provide this if you previously ran annotate with the -ref flag',
        'required': False,
        'default': None,
        'type': str
    }
    skip_makedb = {
        'short': '-skpdb',
        'long': '--skip_makedb',
        'help': 'if a dmnd or blast db already exists for the provided database, '
                'please include this flag to skip this step',
        'required': False,
        'default': False,
        'action': 'store_true'
    }
    no_bidirectional_length = {
        'short': '-nbl',
        'long': '--no_bidirectional_length',
        'help': 'Pseudofinder by default will consider the length_pseudo parameter for assessing both truncated and run-on pseudogenes.'
                'By adding this flag, length_pseudo will only be used for truncated pseudogenes, and pseudofinder will not check for run-on pseudogenes.',
        'required': False,
        'default': False,
        'action': 'store_true'
    }
    use_alignment = {
        'short': '-al',
        'long': '--use_alignment',
        'help': 'Feature currently in development, use at own risk. This flag will make pseudofinder consider '
                'alignment quality between query and hits when considering pseudogene flags.',
        'required': False,
        'default': False,
        'action': 'store_true'
    }

    if names_only:
        to_remove = ['module', 'kwargs', 'names_only', 'to_remove']
        names = list(locals().keys())
        for i in to_remove:
            names.remove(i)
        return names

    if module == 'annotate':
        required_args = [genome, database, outprefix]
        optional_args = [threads, intergenic_length, length_pseudo, shared_hits, evalue, hitcap,
                         contig_ends, intergenic_threshold, reference, diamond, skip_makedb,
                         no_bidirectional_length, use_alignment, perc_id, perc_cov, evalue, max_dnds]
        deprecated_args = [distance]

    elif module == 'reannotate':
        required_args = [genome, logfile, outprefix]
        optional_args = [length_pseudo, shared_hits, intergenic_threshold, max_dnds, max_ds, min_ds, dnds_out,
                         no_bidirectional_length, use_alignment, perc_id, perc_cov, evalue, max_dnds]
        deprecated_args = [distance]

    elif module == 'selection':
        required_args = [prots, genes, ref_prots, ref_genes, reference]
        optional_args = [control_file, outdir, search_engine, min_ds, max_ds, max_dnds, threads, ref_contigs, skip]
        deprecated_args = []

    elif module == 'sleuth':
        required_args = [ref_genes, genome, ref_gff]
        optional_args = [outdir, threads, perc_id, perc_cov, evalue, max_dnds, length_pseudo]
        deprecated_args = []

    elif module == 'genome_map':
        required_args = [genome, gff, outprefix]
        optional_args = []
        deprecated_args = []

    elif module == 'visualize':
        required_args = [logfile, outprefix]
        optional_args = [intergenic_threshold, resolution, keep_files]
        deprecated_args = [distance]

    elif module == 'test':
        required_args = [database]
        optional_args = [genome, diamond, threads]
        deprecated_args = []

    elif module == 'interactive':
        required_args = [annotated_genome]
        optional_args = [outprefix]
        deprecated_args = []

    else:
        print("Module not found. Please check your get_args() function call.")
        exit()

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     usage=bold(usage_message(module, required_args)))
    always_required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Adjustable parameters')
    deprecated = parser.add_argument_group('Deprecated parameters')

    for arg in required_args:
        unpack_arg(always_required, arg)

    for arg in optional_args:
        unpack_arg(optional, arg)

    for arg in deprecated_args:
        unpack_arg(deprecated, arg, deprecated=True)

    # parse_known_args will create a tuple of known arguments in the first position and unknown in the second.
    # We only care about the known arguments, so we take [0].
    args = parser.parse_known_args()[0]
    verify_args(args, deprecated_args)
    return args


def parse_log_args(logfile: str):
    """
    Returns a namespace of all args found in a logfile.
    """
    log_args = Namespace()
    with open(logfile, 'r') as log:
        for line in log.readlines():
            line = line.replace('\n', '')
            sep = '\t'
            for arg in get_args(names_only=True):
                if re.match(arg, line, re.IGNORECASE):
                    item = str(line.split(sep=sep)[1])
                    setattr(log_args, arg, literal_eval(item))
    setattr(log_args, 'log_outprefix', log_args.blastp.replace('_proteome.faa.blastP_output.tsv', ''))
    return log_args


def reconcile_args(priority_args, secondary_args):
    """
    Combines two instances of argparse namespace into one.
    First set of args takes priority
    Second set of args will be considered if the attribute is None or does not exist in the first set
    """
    args = priority_args
    for k, v in vars(secondary_args).items():
        try:
            if getattr(args, k) is None:
                setattr(args, k, v)
        except AttributeError:
            setattr(args, k, v)
    return args


def convert_args_to_str(args):
    """
    Converts every attribute in the args namespace to a string for easy printing.
    """
    str_args = Namespace()
    for k, v in vars(args).items():
        setattr(str_args, k, str(v))
    return str_args


def file_dict(args, **kwargs):
    outprefix = kwargs.get('outprefix', None)

    if outprefix:
        base_outfile_name = outprefix
    else:
        base_outfile_name = args.outprefix

    if not base_outfile_name.endswith('_'):
        base_outfile_name += '_'

    file_dict = {
        'base_filename': base_outfile_name,
        'cds_filename': base_outfile_name + "cds.fasta",
        'contigs_filename': base_outfile_name + "contigs.fasta",
        'ref_cds_filename': base_outfile_name + "ref_cds.fasta",
        'proteome_filename': base_outfile_name + "proteome.faa",
        'ref_proteome_filename': base_outfile_name + "ref_proteome.faa",
        'intergenic_filename': base_outfile_name + "intergenic.fasta",
        'input_pseudos_filename': base_outfile_name + "input_pseudos.fasta",
        'blastp_filename': base_outfile_name + "proteome.faa" + ".blastP_output.tsv",
        'blastx_filename': base_outfile_name + "intergenic.fasta" + ".blastX_output.tsv",
        'blastx_pseudos_filename': base_outfile_name + "input_pseudos.fasta" + ".blastX_output.tsv",
        'pseudos_gff': base_outfile_name + "pseudos.gff",
        'new_pseudos_gff': base_outfile_name + "pseudos-2.gff",
        'pseudos_fasta': base_outfile_name + "pseudos.fasta",
        'new_pseudos_fasta': base_outfile_name + "pseudos-2.fasta",
        'intact_gff': base_outfile_name + "intact.gff",
        'new_intact_gff': base_outfile_name + "intact-2.gff",
        'intact_faa': base_outfile_name + "intact.faa",
        'new_intact_faa': base_outfile_name + "intact-2.faa",
        'intact_ffn': base_outfile_name + "intact.ffn",
        'new_intact_ffn': base_outfile_name + "intact-2.ffn",
        'chromosome_map': base_outfile_name + "map.pdf",
        'interactive_bar': base_outfile_name + "interactive_results.html",
        'interactive_map': base_outfile_name + "interactive_map.html",
        'gbk_out': base_outfile_name[:-1] + ".gbk",
        'dnds_out': base_outfile_name + "dnds",
        'log': base_outfile_name + "log.txt",
        'sleuthDir': base_outfile_name + "sleuth",
        'ctl': os.path.dirname(os.path.dirname(__file__)) + "/codeml-2.ctl"
    }
    return file_dict


def write_test_genome_output(file_dict, genome):
    with open(file_dict['base_filename']+"test_genome.gbk", "w") as output_handle:
        SeqIO.write(genome, output_handle, "genbank")
