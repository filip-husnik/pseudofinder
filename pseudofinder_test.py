#!/usr/bin/env python3

import subprocess
import os
from collections import OrderedDict


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
    path_to_pseudofinder = os.path.dirname(__file__) + "/pseudo_finder.py"
    genome = ""
    database = ""
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
    print(command_dict['annotate_command'])
    exit()

if __name__ == '__main__':
    main()
