#!/usr/bin/env bash

## CaPReSe 1.0
## by Emma Wahlberg, 2019
## emma.wahlberg@nrm.se
##   CaPReSe has three functions: convert a tab delimited (TSV) file
##   downloaded from BOLD to fasta, prepare a FASTA-file for analysis
##   with the script FACEPAI, and merge two FASTA-files (e.g.
##   from BOLD and GenBank) and prepare the merged file for analysis in
##   FACEPAI. Detailed instructions in using this script and associated
##   script can be found in the readme file.
##
##   Usage: caprese.sh
##   -h  Print this text.
##   -C  Convert a TSV-file from BOLD to FASTA and
##       filter out non-BIN-sequences.
##       Usage: SCRIPTNAME -C INPUTFILE OUTPUTFILE
##   -P  Prepare a single file for analysis.
##       Usage: SCRIPTNAME -P NAME_OF_SOURCE INPUTFILE
##       NAME_OF_SOURCE = e.g. BOLD or GenBank.
##   -M  Merge two FASTA-files, e.g. from BOLD and GenBank,
##       and prepare for analysis.
##       Usage: SCRIPTNAME -M NAME_OF_SOURCE1 NAME_OF_SOURCE2 INPUTFILE1 INPUTFILE2
##       NAME_OF_SOURCE = e.g. BOLD or GenBank.
##
##   Software prerequisites: vsearch (see readme for version and details).


# Function for spinner to keep user from becoming bored.
show_spinner()
{
  local -r pid="${1}"
  local -r delay='0.75'
  local spinstr='\|/-'
  local temp
  while ps a | awk '{print $1}' | grep -q "${pid}"; do
    temp="${spinstr#?}"
    printf " [%c]  " "${spinstr}"
    spinstr=${temp}${spinstr%"${temp}"}
    sleep "${delay}"
    printf "\b\b\b\b\b\b"
  done
  printf "    \b\b\b\b"
}

# Function for preparing a FASTA-file for analysis..
prepare_file()
{
  echo 'Prepare '"$1"'...'

  # Change - to N in sequences, change space to underscore, add BOLD identifier, keep only ID and name:
  while true;do show_spinner;sleep 1;done &
    LC_CTYPE=C && LANG=C cat "$2" | sed '/^[^>]/s/-/N/g' | sed 's/ /_/g' | sed "s/>/>$1|/g" > "$1"_prepared.fasta
    kill $!; trap 'kill $!' SIGTERM
    echo "$1"' prepared.'
}

# Function for merging and dereplication two FASTA-files.
merge_files()
{
  echo 'Merging '"$1"' and '"$2"'...'

  TMP_CONC=$(mktemp)
  while true;do show_spinner;sleep 1;done &
    cat "$1"_prepared.fasta "$2"_prepared.fasta > "${TMP_CONC}"
    kill $!; trap 'kill $!' SIGTERM
    echo "$1"' and '"$2"' merged.'

  echo 'Dereplicate sequences using vsearch, write final file...'
  vsearch -derep_fulllength "${TMP_CONC}" --output "$1"_"$2"_merged.fasta
  echo -e '\nMerged files dereplicated and written to: '"$1"'_'"$2"'_merged.fasta'

  rm -f "${TMP_CONC}"
}

# Function for converting a TSV-file from BOLD to FASTA-format.
convert_TSV()
{
  echo 'Converting BOLD TSV-file to FASTA-file...'

  TMP_CONV=$(mktemp)
  while true;do show_spinner;sleep 1;done &
    awk -F "\"*\t\"*" '$8 != "" {print $1,"|",$8,"|",$22,"|",$71,"|",$55,"|",$10,"|",$12,"|",$14,"|",$16,"|",$18,"|",$20,"|",$72}' < "$1" > "${TMP_CONV}"
    sed -i 1d "${TMP_CONV}"
    sed -i 's/^/>/' "${TMP_CONV}"
    sed -i 's/\(.*\) | /\1\n/' "${TMP_CONV}"
    sed -i 's/ | /|/g' "${TMP_CONV}"
    cat "${TMP_CONV}" > "$2"
    kill $!; trap 'kill $!' SIGTERM
    echo "$1"' converted and written to '"$2"'.'

  rm -f "${TMP_CONV}"
}

# Initate script, define variable and read what the user wants.

OPTION=$1

if [ $OPTION == '-M' ]
then
  echo -en "\nChecking dependencies... "
  for name in vsearch
    do
      [[ $(which $name 2>/dev/null) ]] || { echo -en "\n$name needs to be installed.";deps=1; }
    done
  [[ $deps -ne 1 ]] && echo -en "OK\n\n" || { echo -en "\nInstall the above and rerun this script.\n";exit 1; }
  prepare_file $2 $4
  prepare_file $3 $5
  merge_files $2 $3
elif [ $OPTION == '-P' ]
then
  prepare_file $2 $3
elif [ $OPTION == '-C' ]
then
  convert_TSV $2 $3
elif [ $OPTION == '-h' ]
then
  abspath=$(cd ${0%/*} && echo $PWD/${0##*/})
  grep '^\##.*$' $abspath
  exit 1;
else
  echo -e '\nInvalid option! Use -h to show help text, or consult the documentation.\n'
  exit 1;
fi

# End script.

echo 'Finished!'
