# Plant Taxon Name Cleaner

This Python program edits plant taxon names in a way that best aligns with modern databases while adhering to international plant taxonomy standards.

Briefly, it can standardize spacing and capitalization in scientific names and authorities, split taxon names into their parts, replace common misspellings of genera, compare genera with lists of vascular plant, fungus, and bryophyte genera, replace authorities with standard forms, replace special characters, and remove extraneous information about taxonomic alterations/circumscriptions, among other things.

Please read the guide (`Guide_PlantTaxonNameCleaner.pdf`) for much more comprehensive instructions.

## Author

Nicole L. Kinlock (nlkinlock@gmail.com)

## System requirements

Python 3, and Python modules `numpy` and `pandas`

## Installation

Install from GitHub:
   
```bash
git clone https://github.com/nlkinlock/PlantTaxonNameCleaner
```

## Usage

Plant taxon names are provided in a `.csv` input file. The Python program `cleanTaxonNames.py` is executed via a Bash script. To run the program using the provided input file, `Example1.csv`, enter the following from the installed directory:

```bash
./runPlantTaxonNameCleaner.sh Example1.csv -oc
```

To run the program with the provided input file, `Example2.csv`:

```bash
./runPlantTaxonNameCleaner.sh Example2.csv -ac
```

The program will write a cleaned taxon name output `.csv` file as well as a detailed console output `.dat` file. See the guide (`Guide_PlantTaxonNameCleaner.pdf`) for complete usage details and flag meanings.