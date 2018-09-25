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
                                     usage='\033[1m' + "[pseudo_finder.py test -db DATABASE] or "
                                                       "[pseudo_finder.py test --help] for more options." + '\033[0m')

    # Required argument
    required = parser.add_argument_group('\033[1m' + 'Required arguments' + '\033[0m')
    required.add_argument('-db', '--database',
                          help='Please provide name (if $BLASTB is set on your system)'
                               ' or absolute path of your blast database.',
                          required=True)

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
    Tests the given pseudofinder command to make sure that:
        1.  The command runs without an error
        2.  The command produces the expected files
    """
    print("%s\tTesting the %s command.\nFull shell command: %s" % (current_time(), command_name, full_command))

    # 1. The command runs without an error
    try:
        subprocess.run(full_command, shell=True, check=True)
    except subprocess.CalledProcessError:
        print("%s\tCommand failure: %s" % (current_time(), command_name))

    # 2. The command produces the expected files
    # TODO: Fill this in

def main():
    # Inputs for all the test commands
    args = get_args()
    path_to_test_data = os.path.dirname(__file__) + "/test/"
    path_to_pseudofinder = os.path.dirname(__file__) + "/pseudo_finder.py"
    folder_name = manage_folders(path_to_test_data)
    genome = path_to_test_data + "candidatus_tremblaya_princeps_PCIT.gbf"
    database = args.database
    output_prefix = folder_name+"_test"
    blastp_file = ""
    blastx_file = ""
    gff_file = ""
    hitcap = ""

    # A dictionary to store the names and shell commands for each section of pseudofinder
    command_dict = OrderedDict()
    command_dict['annotate_command'] = "python3 %s annotate -g %s -db %s -op %s" % (
        path_to_pseudofinder, genome, database, output_prefix)
    command_dict['visualize_command'] = "python3 %s visualize -g %s -op %s -p %s -x %s -hc %s" % (
        path_to_pseudofinder, genome, output_prefix, blastp_file, blastx_file, hitcap)
    command_dict['map_command'] = "pyton3 %s map -g %s -gff %s -op %s" % (
        path_to_pseudofinder, genome, gff_file, output_prefix)

    # TODO: Also check if dependencies are installed.
    for command in command_dict:
        test_command(command_name=command, full_command=command_dict[command])


if __name__ == '__main__':
    main()
