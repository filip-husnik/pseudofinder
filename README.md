# Pseudo finder
Pseudo finder is a python3 script that detects pseudogene candidates in annotated genbank files of bacterial/archaeal genomes. It was tested mostly on genbank (.gbf/.gbk) files annotated by Prokka [https://github.com/tseemann/prokka] with the --compliant flag.

## Authors
Mitch Syberg-Olsen & Filip Husnik

University of British Columbia, Vancouver, Canada

## Versions and changes


## Getting Started

These instructions will (hopefully) get you the pipeline up and running on your local machine.

### Prerequisites

Installation requirements: python3, pip3 (or any other way how to install python libraries), argparse, biopython, ncbi-blast+

Databases: NCBI-NR protein database (or similar such as SwissProt) formatted for blastp/blastx searches.

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
## Output of pseudo finder

## Contributing

We appreciate any critical comments or suggestions for improvements. Please raise issues or submit pull requests.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* This code was inspired by ...

## References
Recognizing the pseudogenes in bacterial genomes: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1142405/

Taking the pseudo out of pseudogenes: https://www.ncbi.nlm.nih.gov/pubmed/25461580

**Pseudogene detection in a nascent stage of symbiosis (several examples from the Sodalis clade):**

Mobile genetic element proliferation and gene inactivation impact over the genome structure and metabolic capabilities of Sodalis glossinidius, the secondary endosymbiont of tsetse flies: https://www.ncbi.nlm.nih.gov/pubmed/20649993

A Novel Human-Infection-Derived Bacterium Provides Insights into the Evolutionary Origins of Mutualistic Insect–Bacterial Symbioses: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3499248/

Genome Degeneration and Adaptation in a Nascent Stage of Symbiosis: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3914690/

Repeated replacement of an intrabacterial symbiont in the tripartite nested mealybug symbiosis: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5027413/

Large scale and significant expression from pseudogenes in Sodalis glossinidius - a facultative bacterial endosymbiont: https://www.biorxiv.org/content/early/2017/07/23/124388

