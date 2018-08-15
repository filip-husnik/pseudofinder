#!/usr/bin/env python3

import subprocess
import os
import argparse
import re
from time import localtime, strftime
from collections import OrderedDict


def current_time() -> str:
    """Returns the current time. When this function was executed."""
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


def test_command(command_name: str, full_command: str):
    """
    Tests the given pseudofinder command to make sure that:
        1.  The command runs without an error
        2.  The command produces the expected files
    """
    print("Testing the %s command.\nFull shell command: %s" % (command_name, full_command))
    subprocess.call(full_command, shell=True)


def main():
    # Inputs for all the test commands
    args = get_args()
    path_to_test_data = os.path.dirname(__file__) + "/test/"
    path_to_pseudofinder = os.path.dirname(__file__) + "/pseudo_finder.py"
    folder_name = path_to_test_data+strftime("%Y%m%d"+"_results")
    genome = path_to_test_data + "candidatus_tremblaya_princeps_PCIT.gbf"
    database = args.database
    output_prefix = ""
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

    # Make the directory
    try:
        os.makedirs(folder_name)
    except FileExistsError:
        current_folders = os.listdir(path_to_test_data)
        results_folders = [folder for folder in current_folders if re.match(strftime("%Y%m%d"), folder)]
        folder_numbers = [int(folder[-1]) for folder in results_folders if re.match("[0-9]", folder[-1])]

        if not folder_numbers:
            os.makedirs(folder_name + "_1")
        else:
            biggest_folder_number = sorted(folder_numbers, key=int, reverse=True)[0]
            new_number = biggest_folder_number + 1
            os.makedirs(folder_name + "_" + str(new_number))
    exit()

if __name__ == '__main__':
    main()
