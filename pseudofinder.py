#!/usr/bin/env python3

from sys import argv, stderr
from modules import *

"""
pseudofinder.py: A script to find pseudogene candidates in annotated genbank files.
Tested mostly on .gbk files annotated by Prokka with the --compliant flag.

Installation requirements: Please see github repo for detailed explanation of requirements.
 """
__author__ = "Mitch Syberg-Olsen & Filip Husnik"
__version__ = "0.11"
__maintainer__ = "Filip Husnik"
__email__ = "filip.husnik@gmail.com"

errorMessage = "Options: pseudofinder.py [ annotate | reannotate | visualize | dnds | test | help ]\n"

try:
    argv[1]
except IndexError:
    stderr.write(errorMessage)
    exit()

if argv[1] == "annotate":
    annotate.main()
elif argv[1] == "reannotate":
    reannotate.main()
elif argv[1] == "visualize":
    visualize.main()
elif argv[1] == "dnds":
    dnds.main()
elif argv[1] == "test":
    pseudofinder_test.main()
elif argv[1] == "help":
    stderr.write("\tpseudofinder.py annotate: Flags candidate pseudogenes.\n"
                 "\tpseudofinder.py reannotate: Begins the annotate pipeline post-BLAST.\n"
                 "\tpseudofinder.py visualize: Generates a 3D plot to visualize different combinations of "
                 "settings.\n"
                 "\tpseudofinder.py dnds: dN/dS analysis via pairwise comparison against a reference genome. "
                 "Pseudogenes inferred from relaxed selection.\n"
                 "\tpseudofinder.py test: Runs all commands on a test dataset and checks that the outputs "
                 "are as expected.\n")
    exit()
else:
    stderr.write(errorMessage)
    exit()
