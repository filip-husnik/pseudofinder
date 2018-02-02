# Pseudo finder
Pseudo finder is a python3 script to find pseudogene candidates in annotated genbank files.
It was tested mostly on genbank (.gbf/.gbk) files annotated by Prokka with the --compliant flag.

## Authors
Mitch Syberg-Olsen & Filip Husnik
University of British Columbia, Vancouver, Canada

## Versions and changes


## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

Installation requirements: python3, pip3 (or any other way how to install python libraries), argparse, biopython, ncbi-blast+

Databases: NCBI-NR protein database (or similar such as SwissProt) formatted for blastp/blastx searches

Input files: Genome sequence with annotations in the genbank (.gbf/.gbk) format.


### Installing

A step by step series of commands to install all dependencies:

Installation of python3, pip3, and ncbi-blast+ on Ubuntu:

```
sudo apt-get update
sudo apt-get install python3
sudo apt-get install python3-pip
sudo apt-get install ncbi-blast+
```

Installation of python3 libraries on Ubuntu:

```
sudo pip3 install argparse
sudo pip3 install biopython
```
## Running pseudo finder

```
# Run with test data and 16 processors (for BlastX/BlastP searches)
python3 pseudo_finder.py --genome GENOME.GBF --output PREFIX --database NR --threads 16 
```


## Contributing

We appreciate any critical comments or suggestions for improvements. Please raise issues or submit pull requests.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* This code was inspired by ...

## References
Recognizing the pseudogenes in bacterial genomes: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1142405/

