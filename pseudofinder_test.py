#!/usr/bin/env python3
##############  UNDER CONSTRUCTION    ##############
##############  UNDER CONSTRUCTION    ##############
##############  UNDER CONSTRUCTION    ##############
import subprocess, os, argparse
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
    CommandDict = OrderedDict()
    CommandDict['annotate_command'] = "python3 %s annotate -g %s -db %s -op %s" % (path_to_pseudofinder, genome, database, output_prefix)
    CommandDict['visualize_command'] = "python3 %s visualize -g %s -op %s -p %s -x %s -hc %s" % (path_to_pseudofinder, genome, output_prefix, blastp_file, blastx_file, hitcap)
    CommandDict['map_command'] = "pyton3 %s map -g %s -gff %s -op %s" % (path_to_pseudofinder, genome, gff_file, output_prefix)

    print(CommandDict['annotate_command'])
    exit()

main()