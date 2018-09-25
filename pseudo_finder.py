#!/usr/bin/env python3

from sys import argv, stderr
from modules import annotate, visualize, genome_map, pseudofinder_test  # all pseudofinder modules

"""
pseudo_finder.py: A script to find pseudogene candidates in annotated genbank files.
Tested mostly on .gbk files annotated by Prokka with the --compliant flag.

Installation requirements:
    *python3
    *3rd party libraries: biopython, plotly, pandas, numpy, reportlab
    *libraries from the Python standard library: see annotate.py, visualize.py and genome_map.py

 """
__author__ = "Mitch Syberg-Olsen & Filip Husnik"
__version__ = "0.09"
__maintainer__ = "Filip Husnik"
__email__ = "filip.husnik@gmail.com"

errorMessage = "Options: pseudo_finder.py [ annotate | visualize | map | test | help ]\n"

try:
    if argv[1] == "annotate":
        annotate.main()
    elif argv[1] == "visualize":
        visualize.main()
    elif argv[1] == "map":
        genome_map.main()
    elif argv[1] == "test":
        pseudofinder_test.main()
    elif argv[1] == "help":
        stderr.write("\tpseudofinder.py annotate: Flags candidate pseudogenes.\n"
                     "\tpseudofinder.py visualize: Generates a 3D plot to visualize different combinations of "
                     "settings.\n"
                     "\tpseudofinder.py map: Generates a chromosome map to show where pseudogenes have been called.\n"
                     "\tpseudofinder.py test: Runs all commands on a test dataset and checks that the outputs "
                     "are as expected.\n")
        exit()
    else:
        stderr.write(errorMessage)
        exit()
except IndexError:  # catches if the user does not provide any arguments
    stderr.write(errorMessage)
    exit()
