#!/usr/bin/env python3

from .common import get_args, bold, print_with_time
import os
from os.path import dirname
import re
import subprocess
from collections import OrderedDict
from time import localtime, strftime


def current_time() -> str:
    """Returns the current time when this function was executed."""
    return str(strftime("%Y-%m-%d %H:%M:%S", localtime()))


def manage_folders(path: str):
    """Creates a folder to store test results and returns the name."""

    folder_name = path+strftime("%Y%m%d"+"_results")

    try:  # simplest solution
        os.makedirs(folder_name)

    except FileExistsError:  # unless the folder already exists.
        current_folders = os.listdir(path)
        results_folders = [folder for folder in current_folders if re.match(strftime("%Y%m%d"), folder)]
        folder_numbers = [int(folder[-1]) for folder in results_folders if re.match("[0-9]", folder[-1])]

        if not folder_numbers:  # If no numbered folders exist yet, start with "folder_name_1"
            folder_name = folder_name + "_1"
            os.makedirs(folder_name)

        else:  # increment the number by 1, ie. "folder_name_2"
            biggest_folder_number = sorted(folder_numbers, key=int, reverse=True)[0]
            new_number = biggest_folder_number + 1
            folder_name = folder_name + "_" + str(new_number)
            os.makedirs(folder_name)

    return folder_name


def test_command(command_name: str, full_command: str):
    """
    Tests the given pseudofinder command to make sure that the command runs without an error.
    """
    print_with_time(bold('Testing the %s command.' % command_name))
    print_with_time('Full shell command: %s' % full_command)

    try:
        subprocess.run(full_command, shell=True, check=True)
    except subprocess.CalledProcessError:
        print_with_time('Command failure: ' + command_name)


def main():
    # Inputs for all the test commands
    args = get_args('test')

    path_to_master = dirname(dirname((__file__)))
    path_to_test_data = path_to_master + "/test/"
    path_to_pseudofinder = path_to_master + "/pseudofinder.py"
    folder_name = manage_folders(path_to_test_data)
    genome_name = "candidatus_tremblaya_princeps_PCIT.gbf"
    genome_full_path = path_to_test_data + genome_name
    output_prefix = folder_name + "/test"
    blastp_file = output_prefix + "_proteome.faa.blastP_output.tsv"
    blastx_file = output_prefix + "_intergenic.fasta.blastX_output.tsv"
    log_file = output_prefix + "_log.txt"

    ctl_full_path = path_to_master + "/codeml-2.ctl"
    ref_pep = path_to_test_data + "Mycobacterium_tuberculosis_H37Rv-subset.faa"
    ref_nuc = path_to_test_data + "Mycobacterium_tuberculosis_H37Rv-subset.ffn"
    pep = path_to_test_data + "Mycobacterium_leprae-subset.faa"
    nuc = path_to_test_data + "Mycobacterium_leprae-subset.ffn"
    dndsOutput = "dnds_output"

    # A dictionary to store the names and shell commands for each section of pseudofinder
    command_dict = OrderedDict()
    command_dict['Annotate'] = "python3 %s annotate -g %s -db %s -op %s -t %s --diamond" % (
        path_to_pseudofinder, genome_full_path, args.database, output_prefix, args.threads)
    command_dict['Visualize'] = "python3 %s visualize -g %s -op %s -p %s -x %s -log %s" % (
        path_to_pseudofinder, genome_full_path, output_prefix, blastp_file, blastx_file, log_file)
    command_dict['Reannotate'] = "python3 %s reannotate -g %s -p %s -x %s -log %s -op %s" % (
        path_to_pseudofinder, genome_full_path, blastp_file, blastx_file, log_file, output_prefix)
    # command_dict['dnds'] = "python3 %s dnds -ra %s -rn %s -a %s -n %s -ctl %s -out %s" % (
    #     path_to_pseudofinder, ref_pep, ref_nuc, pep, nuc, ctl_full_path, dndsOutput)

    for command in command_dict:
        test_command(command_name=command, full_command=command_dict[command])

    # os.system("rm -r " + dndsOutput) # TODO: remove?

if __name__ == '__main__':
    main()

