#!/usr/bin/env python3

"""
pseudofinder.py: A script to find pseudogene candidates in annotated genbank files.
Tested mostly on .gbk files annotated by Prokka with the --compliant flag.

Installation requirements: Please see github repo for detailed explanation of requirements.
 """
__author__ = "Mitch Syberg-Olsen, Arkadiy Garber & Filip Husnik"
__version__ = "1.00"
__maintainer__ = "Filip Husnik"
__email__ = "filip.husnik@gmail.com"

try:
    from sys import argv, stderr
    from modules import *

except ModuleNotFoundError as err:
    stderr.write(f"Pseudofinder encountered an error: {str(err)}.\n"
                 f"If you are a conda user, please enable the Pseudofinder environment with the following command:\n"
                 f"conda activate pseudofinder\n"
                 f"If you are not a conda user, please check that you have installed the required dependencies.\n")
    exit()

except ImportError as err:
    stderr.write(f"Pseudofinder encountered an error: {str(err)}.\n"
                 f"Please check that no files are missing from the pseudofinder/modules folder.\n")
    exit()

errorMessage = "Options: " + common.bold("pseudofinder.py [ annotate | reannotate | visualize | sleuth | test | help ]\n")

try:
    module = argv[1]
except IndexError:
    stderr.write(errorMessage)
    exit()

try:
    if module == "annotate":
        annotate.main()
    elif module == "reannotate":
        reannotate.main()
    elif module == "visualize":
        visualize.main()
    elif module == "sleuth":
        sleuth.main()
    elif module == "test":
        pseudofinder_test.main()
    elif module == "interactive":
        interactive.main()
    elif module == "help":
        stderr.write("\tpseudofinder.py annotate: Flags candidate pseudogenes.\n"
                     "\tpseudofinder.py reannotate: Begins the annotate pipeline post-BLAST.\n"
                     "\tpseudofinder.py visualize: Generates a 3D plot to visualize different combinations of "
                     "settings.\n"
                     "\tpseudofinder.py sleuth: pairwise comparison against a reference genome."
                     "Pseudogenes inferred from relaxed selection.\n"
                     "\tpseudofinder.py test: Runs all commands on a test dataset and checks that the outputs "
                     "are as expected.\n")
        exit()
    else:
        stderr.write(errorMessage)
        exit()

except KeyboardInterrupt:
    print("\n")
    common.print_with_time("Pseudofinder process cancelled by user. Exiting now.")
    exit()
