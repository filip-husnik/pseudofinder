#!/usr/bin/env python3
"""
pseudo_finder.py: A script to find pseudogene candidates in annotated genbank files.
Tested mostly on .gbk files annotated by Prokka with the --compliant flag.

Installation requirements:
    python3
    3rd party libraries: biopython
    libraries from the Python standard library: see annotate.py and visualize.py

 """
__author__ = "Mitch Syberg-Olsen & Filip Husnik"
__version__ = "0.08"
__maintainer__ = "Filip Husnik"
__email__ = "filip.husnik@gmail.com"

from sys import argv, stderr
import annotate
import visualize

print("Please note: the pseudo_finder code is currently being reorganized. Up to date documentation will be written shortly.")

errormessage=("Options: [pseudo_finder.py annotate] or [pseudo_finder.py visualize]\n")

try:
    if argv[1] == "annotate":
        annotate.main()
    elif argv[1] == "visualize":
        visualize.main()
    else:
        stderr.write(errormessage)
except IndexError:
    stderr.write(errormessage)