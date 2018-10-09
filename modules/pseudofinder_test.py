#!/usr/bin/env python3

import argparse
import os
import re
import subprocess
from collections import OrderedDict
from time import localtime, strftime


def current_time() -> str:
    """Returns the current time when this function was executed."""
    return str(strftime("%Y-%m-%d %H:%M:%S", localtime()))


def get_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     usage='\033[1m' + "[pseudofinder.py test -db DATABASE] or "
                                                       "[pseudofinder.py test --help] for more options." + '\033[0m')

    # Required argument
    required = parser.add_argument_group('\033[1m' + 'Required arguments' + '\033[0m')
    required.add_argument('-db', '--database',
                          help='Please provide name (if $BLASTB is set on your system)'
                               ' or absolute path of your blast database.',
                          required=True)
    optional = parser.add_argument_group('\033[1m' + 'Optional arguments' + '\033[0m')
    optional.add_argument('-t', '--threads', default=4,
                          help='Please provide total number of threads to use for blast, default is 4.')

    # parse_known_args will create a tuple of known arguments in the first position and unknown in the second.
    # We only care about the known arguments, so we take [0].
    args = parser.parse_known_args()[0]

    return args


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
    print("\033[1m"+"%s\tTesting the %s command.\n\033[0m"
          "%s\tFull shell command: %s" % (current_time(), command_name, current_time(), full_command))

    try:
        subprocess.run(full_command, shell=True, check=True)
    except subprocess.CalledProcessError:
        print("\033[1m"+"%s\tCommand failure: %s\033[0m" % (current_time(), command_name))


def main():
    # Inputs for all the test commands
    args = get_args()
    path_to_test_data = os.path.dirname(__file__) + "/../test/"
    path_to_pseudofinder = os.path.dirname(__file__) + "/../pseudofinder.py"
    folder_name = manage_folders(path_to_test_data)
    genome_name = "candidatus_tremblaya_princeps_PCIT.gbf"
    genome_full_path = path_to_test_data + genome_name
    output_prefix = folder_name+"/test"
    blastp_file = output_prefix+"_proteome.faa.blastP_output.tsv"
    blastx_file = output_prefix+"_intergenic.fasta.blastX_output.tsv"
    log_file = output_prefix+"_log.txt"

    # A dictionary to store the names and shell commands for each section of pseudofinder
    command_dict = OrderedDict()
    command_dict['Annotate'] = "python3 %s annotate -g %s -db %s -op %s -t %s" % (
        path_to_pseudofinder, genome_full_path, args.database, output_prefix, args.threads)
    command_dict['Visualize'] = "python3 %s visualize -g %s -op %s -p %s -x %s -log %s" % (
        path_to_pseudofinder, genome_full_path, output_prefix, blastp_file, blastx_file, log_file)
    command_dict['Reannotate'] = "python3 %s reannotate -g %s -p %s -x %s -log %s -op %s" % (
        path_to_pseudofinder, genome_full_path, blastp_file, blastx_file, log_file, output_prefix)

    for command in command_dict:
        test_command(command_name=command, full_command=command_dict[command])


if __name__ == '__main__':
    main()
