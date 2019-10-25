#!/usr/bin/env bash

## FACEPAI 1.0
## by Emma Wahlberg, 2019
## emma.wahlberg@nrm.se
##   FACEPAI will take eDNA sequences in demultiplexed FASTQ-format and pool several
##   replicates, filter, merge reads, dereplicate, cluster and finally perform
##   taxonomic identification. FACEPAI is constructed to be fast and to be
##   used with minimum effort and hands on time. Results in table format can
##   thereafter be manually controlled, checked and validated. For detailed
##   information and instructions, consult the readme file.
##
##   Usage: FACEPAI SAMPLE_NAME FORWARD_IDENTIFIER REVERSE_IDENTIFIER
##           SAMPLE_NAME = the name of the sample.
##           FORWARD_IDENTIFIER = string that identifies the file or
##           files containing the forward reads in FASTQ format,
##           this may be one or multiple files, as long as the identifier
##           string is the same.
##           REVERSE_IDENTIFIER = same as above but for reverse reads.
##   -h  Print this text.
##
##   Preferences for primer sequences, minminum length, reference database,
##   and header variables are set in the file "options.config".
##   Software prerequisites: fastp, vsearch, cutadapt, swarm, blastn,
##   and makdeblastdb (see readme for versions and details).

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

# Display help if requested.

if [ "$1" == '-h' ]
then
  abspath=$(cd ${0%/*} && echo $PWD/${0##*/})
  grep '^\##.*$' $abspath
  exit 1;
fi

# Read options file

ABSPATH=$(readlink -f $0)
ABSDIR=$(dirname $ABSPATH)
. "$ABSDIR/options.config"

# Define variables

VSEARCH=$(which vsearch)
THREADS=4
ENCODING=33

INPUT=$1
INPUTENDINGF=$2
INPUTENDINGR=$3
INPUTS=1
OUTPUT=${INPUT}
REPFILEENDINGF="_pooled_F.fastq"
REPFILEENDINGR="_pooled_R.fastq"
POOLEDF="$INPUT$REPFILEENDINGF"
POOLEDR="$INPUT$REPFILEENDINGR"
FILTEREDENDINGF="_pooled_F_filtered.fastq"
FILTEREDENDINGR="_pooled_R_filtered.fastq"
FILTEREDNAMEF="$INPUT$FILTEREDENDINGF"
FILTEREDNAMER="$INPUT$FILTEREDENDINGR"
BLASTENDING="_BLAST_results.tab"
BLASTRESULTS="$INPUT$BLASTENDING"

# Setup log file

LOGFILENAME="${INPUT}_analysis.log"
exec &> >(tee -a "$LOGFILENAME")

# Log time

echo -en "\nAnalysis started at:\n"
date +"Date: %d/%m/%Y Time: %H:%M:%S"

# Check for required dependencies

echo -en "\nChecking dependencies... "
for name in fastp vsearch cutadapt swarm blastn
  do
    [[ $(which $name 2>/dev/null) ]] || { echo -en "\n$name needs to be installed.";deps=1; }
  done
[[ $deps -ne 1 ]] && echo -en "OK\n" || { echo -en "\nInstall the above and rerun this script.\n";exit 1; }

# Check input

echo -en "\nChecking input... "
if [ -z "$INPUT" ]; then
  echo -en "\nA sample name needs to be provided.";INPUTS=$((INPUTS+1));
fi
if [ -z "$INPUTENDINGF" ]; then
  echo -en "\nFile ending of forward reads needs to be provided.";INPUTS=$((INPUTS+1));
fi
if [ -z "$INPUTENDINGR" ]; then
  echo -en "\nFile ending of reverse reads needs to be provided.";INPUTS=$((INPUTS+1));
fi
if [ $INPUTS -ne 1 ]; then
  echo -en "\nCheck above inputs and rerun this script.\nFor help, run script with argument -h or consult the documentation.\n";exit 1;
else
  echo -en "OK\n";
fi

# Merge, filter and pool samples

read -p $'\n'$"Merging reads, filtering and pooling samples..."$'\n' -t 1

RAW_FILES=*${INPUTENDINGF}
for f in $RAW_FILES
do
  FORWARDREAD=$f
  REVERSEREAD=${f/$INPUTENDINGF/$INPUTENDINGR}
  REPMERGED=${f/$INPUTENDINGF/_merged.fastq}
  read -p $'\n'$"Processing ${FORWARDREAD} and ${REVERSEREAD}."$'\n' -t 0.5
  read -p $'\n'$""$'\n' -t 0.5
  fastp -m -u 100 -q 15 -c -L --overlap_len_require 10 --overlap_diff_percent_limit 40 -i "${FORWARDREAD}" -I "${REVERSEREAD}" --merged_out "${REPMERGED}"
done

cat *_merged.fastq > "${INPUT}"

# Remove primers

read -p $'\n'$"Primers and minimum length defined:"$'\n' -t 0.5

echo Forward primer: ${PRIMER_F}
echo Reverse complemented primer: ${PRIMER_R}
echo Minimum length: ${MIN_LENGTH}

read -p $'\n'$"Removing primers..."$'\n' -t 1

MIN_F=$(( ${#PRIMER_F} * 2 / 3 ))
MIN_R=$(( ${#PRIMER_R} * 2 / 3 ))


# Define binaries, temporary files and output files

CUTADAPT="$(which cutadapt) --discard-untrimmed --minimum-length ${MIN_LENGTH}"
VSEARCH=$(which vsearch)
INPUT_REVCOMP=$(mktemp)
TMP_FASTQ=$(mktemp)
TMP_FASTQ2=$(mktemp)
TMP_FASTA=$(mktemp)
OUTPUT=$(mktemp)
QUALITY_FILE="${INPUT/.fastq/.qual}"

# Reverse complement fastq file

"${VSEARCH}" --quiet \
             --fastx_revcomp "${INPUT}" \
             --fastqout "${INPUT_REVCOMP}"

    LOG="${INPUT}.log"
    FINAL_FASTA="${INPUT}.fas"

# Trim forward & reverse primers (search normal and antisens)
cat "${INPUT}" "${INPUT_REVCOMP}" | \
      ${CUTADAPT} -g "${PRIMER_F}" -O "${MIN_F}" - 2>> "${LOG}" | \
      ${CUTADAPT} -a "${PRIMER_R}" -O "${MIN_R}" - 2>> "${LOG}" > "${TMP_FASTQ}"

read -p $'\n'$"Additional filtering using VSEARCH..."$'\n' -t 0.5
read -p $'\n'$"Discarding sequences with N:s, converting to FASTA and dereplicate."$'\n' -t 1

# Discard sequences containing Ns, add expected error rates
"${VSEARCH}" \
  --quiet \
  --fastq_filter "${TMP_FASTQ}" \
  --fastq_maxns 0 \
  --relabel_sha1 \
  --eeout \
  --fastqout "${TMP_FASTQ2}" 2>> "${LOG}"

# Discard sequences containing Ns, convert to fasta
"${VSEARCH}" \
  --quiet \
  --fastq_filter "${TMP_FASTQ}" \
  --fastq_maxns 0 \
  --fastaout "${TMP_FASTA}" 2>> "${LOG}"

# Dereplicate at the study level
"${VSEARCH}" \
  --quiet \
  --derep_fulllength "${TMP_FASTA}" \
  --minuniquesize 10 \
  --sizeout \
  --fasta_width 0 \
  --relabel_sha1 \
  --output "${FINAL_FASTA}" 2>> "${LOG}"

# Discard quality lines, extract hash, expected error rates and read length
sed 'n;n;N;d' "${TMP_FASTQ2}" | \
  awk 'BEGIN {FS = "[;=]"}
  {if (/^@/) {printf "%s\t%s\t", $1, $3} else {print length($1)}}' | \
  tr -d "@" >> "${OUTPUT}"

# Produce the final quality file

read -p $'\n'$"Cleaning up..."$'\n' -t 1

sort -k3,3n -k1,1d -k2,2n "${OUTPUT}" | \
    uniq --check-chars=40 > "${QUALITY_FILE}"

# Clean

rm -f "${INPUT_REVCOMP}" "${TMP_FASTQ}" "${TMP_FASTA}" "${TMP_FASTQ2}" "${OUTPUT}"

read -p $'\n'$"Finished preparing files."$'\n' -t 0.5
read -p $'\n'$"Moving on to swarming and BLASTing..."$'\n' -t 1

SWARM=$(which swarm)
TMP_FASTA=$(mktemp --tmpdir=".")
FE=".fasta"
FINAL_FASTA="$INPUT$FE"

# Read sequences

cat *.fas > "${TMP_FASTA}"

# Dereplicate (vsearch)

read -p $'\n'$"Dereplicating..."$'\n' -t 1

"${VSEARCH}" --derep_fulllength "${TMP_FASTA}" \
             --sizein \
             --sizeout \
             --fasta_width 0 \
             --output "${FINAL_FASTA}" > /dev/null

rm -f "${TMP_FASTA}"

# Clustering

read -p $'\n'$"Clustering..."$'\n' -t 1

THREADS=8
TMP_REPRESENTATIVES=$(mktemp --tmpdir=".")
"${SWARM}" \
    -d 1 -f -t ${THREADS} -z \
    -i ${FINAL_FASTA/.fas/_1f.struct} \
    -s ${FINAL_FASTA/.fas/_1f.stats} \
    -w ${TMP_REPRESENTATIVES} \
    -o ${FINAL_FASTA/.fas/_1f.swarms} < ${FINAL_FASTA}

# Sort representatives

read -p $'\n'$"Sorting representatives..."$'\n' -t 1

"${VSEARCH}" --fasta_width 0 \
             --sortbysize ${TMP_REPRESENTATIVES} \
             --output ${FINAL_FASTA/.fas/_1f_representatives.fas}
rm ${TMP_REPRESENTATIVES}

# Chimera checking

read -p $'\n'$"Checking for chimeras..."$'\n' -t 1

REPRESENTATIVES=${FINAL_FASTA/.fas/_1f_representatives.fas}
UCHIME=${REPRESENTATIVES/.fas/.uchime}
"${VSEARCH}" --uchime_denovo "${REPRESENTATIVES}" \
             --uchimeout "${UCHIME}"

# Perform local BLAST

read -p $'\n'$"Completed swarming and chimera checking. Prepare for BLAST!"$'\n' -t 0.5

echo -en "\nBLAST started at:\n"
date +"Date: %d/%m/%Y Time: %H:%M:%S"

read -p $'\n'$"Be aware, this may take some time..."$'\n' -t 0.5

TMP_BLAST=$(mktemp)
TMP_BLAST2=$(mktemp)
TMP_BLAST3=$(mktemp)
TMP_BLAST4=$(mktemp)
TMP_BLAST5=$(mktemp)

sed -i 's/^\([^_]*\(_[^_]*\)\{2\}\).*/\1/' *_1f_representatives.fasta

while true;do show_spinner;sleep 1;done &
  LC_CTYPE=C && LANG=C ionice -c2 -n7 blastn -query *_1f_representatives.fasta -db "${BLASTDB}" -max_target_seqs 10 -outfmt '6 qseqid pident evalue qcovs sseqid' > "${TMP_BLAST}"
  kill $!; trap 'kill $!' SIGTERM

read -p $'\n'$"BLAST complete, preparing results..."$'\n' -t 1

tr '|' \\t < "${TMP_BLAST}" > "${TMP_BLAST2}"

tr '_' ' ' < "${TMP_BLAST2}" > "${TMP_BLAST3}"

awk -vRS=";size=" -vORS="\t" '1' "${TMP_BLAST3}" > "${TMP_BLAST4}"

tr -d ';' < "${TMP_BLAST4}" > "${TMP_BLAST5}"

echo -e "${TABHEADER//,/\\t}" | cat - "${TMP_BLAST5}" > "${BLASTRESULTS}"

rm -f "${TMP_BLAST}" "${TMP_BLAST2}" "${TMP_BLAST3}" "${TMP_BLAST4}" "${TMP_BLAST5}" "${TMP_DB}"

# Log time

echo -en "\nAnalysis finished at:\n"
date +"Date: %d/%m/%Y Time: %H:%M:%S"

read -p $'\n'$"Results complete. Carry on, captain!"$'\n' -t 1
