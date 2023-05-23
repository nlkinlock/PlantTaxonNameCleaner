#!/bin/bash
# runPlantTaxonNameCleaner.sh: run cleanTaxonNames.py, passing CSV file as variable, and write std err to console
usage() {
	# print usage to stderr
	exec 1>&2
	printf "Usage: $0 FilePath -o -a -c\n"
	printf "\tFilePath\tpath to file, e.g., /mnt/d/Nicole/Postdoc/WorldCuP/NewSources/TableToClean.csv or\n\t\t\t/home/nlkinlock/Postdoc/WorldCuP/NewSources/TableToClean.csv (full file path recommended)\n"
	printf "\t-o\tis the file an old (e.g., pre-1900s) text? if this flag is used, authorities are replaced with accepted\n\t\tabbreviations, and (suspected) cultivar names given as infraspecific names are reclassified accordingly\n"
	printf "\t-a\tdo authorities need to be split from taxon names? this only matters if you have a 'TaxonName' column\n"
	printf "\t-c\tshould cultivar names be standardized, e.g., 'Cultivar Name'? if this flag is NOT used,\n\t\tcultivars will be treated as is\n"
	exit 1
}
if [[ "$#" -lt "1" ]] ; then
	usage
fi

FILEPATH=$1; shift
OLDTEXT="False"
AUTHORSPLIT="False"
CV="False"

if [ -n "$1" ]
then :
	while getopts "oac" opt; do
		case $opt in
			o) 
				OLDTEXT="True"
				;;
			a) 
				AUTHORSPLIT="True"
				;;
			c)
				CV="True"
				;;
		esac
	done
	shift $((OPTIND -1))
else :
	printf "no flags, assuming new text, no authorities to split, and no edits to cultivars!\n"
fi

DIR=${FILEPATH%/*}
FILE=${FILEPATH##*/}
BASE=${FILE%.*}

if [[ $FILEPATH != *.csv ]]; then
	echo "file (${FILEPATH}) must be a CSV file!"
	exit
fi

SCRIPTDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

printf "\nin runPlantTaxonNameCleaner.sh (in ${SCRIPTDIR})\n"

printf "\tflags: OLDTEXT = ${OLDTEXT}, AUTHORSPLIT = ${AUTHORSPLIT}, CV = ${CV}\n"	

if [ ! -f "$FILEPATH" ]; then
    echo "file (${FILEPATH}) does not exist!"
	exit
fi

if ! command -v python3 &> /dev/null
then
    echo "python3 is not installed, install and try again"
    exit
fi

DATE=$(date +"%Y%m%d")
TIME=$(date +"%H:%M")
printf "\nrunning cleanTaxonNames.py\n"
printf "\nstarted on ${DATE} at ${TIME}\n"

OUTFILE=${DIR}/${BASE}_Console_${DATE}.dat
printf "\nsaving output to ${OUTFILE}\n"

exec 3> "${OUTFILE}"
python3 ${SCRIPTDIR}/cleanTaxonNames.py ${FILEPATH} ${SCRIPTDIR} ${OLDTEXT} ${AUTHORSPLIT} ${CV} 2>&1 | tr '\r' '\n' >&3
exec 3<&-

echo -en "\007"

DATE=$(date +"%Y%m%d")
TIME=$(date +"%H:%M")
printf "\nfinished on ${DATE} at ${TIME}\n"
