#! /bin/bash


### Pipefail
set -o pipefail



### Argument processing
INPUT=""
OUTPUT="output_file"
VERBOSE=3
CLEAN=0
KEEP=0
START=0
END=100
CURR=0
THREADS=1
PAIRS=0

while getopts 'f:v:cks:e:pt:o:' flag ; do
	case $flag in
		f)
			INPUT=$OPTARG
			;;

		v)
			VERBOSE=$OPTARG
			;;

		c)
			CLEAN=1
			;;

		k)
			KEEP=1
			;;

		s)
			START=${OPTARG}
			;;

		e)
			END=${OPTARG}
			;;

		p)
			PAIRS=1
			;;

		t)
			THREADS=${OPTARG}
			;;

		o)
			OUTPUT=${OPTARG}
			;;

		\?)
			err "Invalid flags. Exiting ... "
			exit 12
			;;

		:)
			err "Flag requires an argument. Exiting ... "
			exit 13
			;;

		*)
			err "Error with flags. Exiting ... "
			exit 14
			;;
	esac
done



MERGED_INPUT_BAM=${OUTPUT}_merged.bam
INDEX_INPUT_BAM=${OUTPUT}_merged.ind.bam

SORTED_INPUT_BAM=${OUTPUT}_merged.ind.sort.bam
SORTED_INPUT_HEAD_TXT=${OUTPUT}_merged.ind.sort.head.txt
SORTED_INPUT_TAIL_TXT=${OUTPUT}_merged.ind.sort.tail.txt
SORTED_INPUT_REMOVE_TXT=${OUTPUT}_merged.ind.sort.remove.txt

REVERT_INPUT_BAM=${OUTPUT}_merged.ind.sort.revert.bam
REVERT_INPUT_BAD_FIRST=${OUTPUT}_merged.ind.sort.revert.first.txt
REVERT_INPUT_BAD_SECOND=${OUTPUT}_merged.ind.sort.revert.second.txt
REVERT_INPUT_BAD_REMOVE=${OUTPUT}_merged.ind.sort.revert.remove.txt

REMOVE_READS_TXT=${OUTPUT}_merged.removeReads.txt
READS_INPUT_BAM=${OUTPUT}_merged.ind.sort.revert.reads.bam

MATE_INPUT_BAM=${OUTPUT}_merged.ind.sort.revert.reads.mate.bam
RG_INPUT_BAM=${OUTPUT}_merged.ind.sort.revert.reads.mate.rg.bam
VALIDATE_INPUT_TXT=${OUTPUT}_merged.ind.sort.revert.reads.mate.rg.validate.txt

READY_BAM=${OUTPUT}_ready.bam

MERGED_BAM=${OUTPUT}_shuffled.bam
CLEAN_BAM=${OUTPUT}_shuffled.clean.bam



TEMP_DIR="${PWD}/tmp"
LOG_DIR="${PWD}/log"
RESULTS_DIR="${PWD}/results"

MAIN_LOG_FILE=${LOG_DIR}/main.log


SCRIPT_DIR="/home/shared/scripts"
CURR_SCRIPT_DIR="${SCRIPT_DIR}/BAMtoFASTQ"


TOOL_DIR="/home/shared/tools"
PICARD_FILE="${TOOL_DIR}/picard.jar"


### Text variables
PASS_TEST_LIGHT="[\e[102mPASSED\e[0m]"
PASS_TEST="[\e[42mPASSED\e[0m]"
FAIL_TEST_LIGHT="[\e[101mFAILED\e[0m]"
FAIL_TEST="[\e[41mFAILED\e[0m]"


# Create a logger for messages. 
DEBUG="\e[94mDEBUG\e[0m"
INFO="INFO \e[0m"
WARN="\e[1m\e[93mWARN \e[0m"
ERROR="\e[1m\e[91mERROR\e[0m"

log() {
	DATE=`date +"%H:%M:%S %Y/%m/%d"`

	case $2 in
		4)
			MSG=${DEBUG}
			;;

		3)
			MSG=${INFO}
			;;

		2)
			MSG=${WARN}
			;;

		1)
			MSG=${ERROR}
			;;

		*)
			MSG=${INFO}
			;;
	esac

	if [[ $2 -le ${VERBOSE}  ]]
	then
		MESSAGE="${MSG}\t${DATE}\t$1"
		echo -e "${MESSAGE}" >> ${MAIN_LOG_FILE}
		echo -e "${MESSAGE}"
	fi
}



if [[ ! -d ${LOG_DIR} ]]
then
	mkdir ${LOG_DIR}
fi


if [ -f ${MAIN_LOG_FILE} ]
then
	mv ${MAIN_LOG_FILE} ${MAIN_LOG_FILE}.$(date +%Y-%m-%d_%H.%M.%S)
	touch ${MAIN_LOG_FILE}
	log "Renaming previous log file" 4
else
	touch ${MAIN_LOG_FILE}
	log "Creating the main log file: ${MAIN_LOG_FILE}" 4
fi



testResult() {
	if [ $1 != "0" ]
	then
		log "${FAIL_TEST}" 2
		exit 1
	else
		log "${PASS_TEST}" 3
	fi
	
	log " " 3
}


ulimit -n $(ulimit -Hn) 2> /dev/null
log "Note: ulimit is $(ulimit -Sn), max limit is $(ulimit -Hn)" 4


### Heading
log "$(printf '#%.0s' $(seq 1 $(($(tput cols)-35)) ))" 3
log " " 3 
log " " 3 
log "Convert BAM to FASTQ with shuffling" 3
log " " 3
log " " 3


### Initial checks
if [[ -z ${INPUT} ]]
then
	log "Error: no input specified with -f. Exiting ... " 1
	exit 10
fi



if [[ ${CLEAN} -ne 0 ]]
then 
	if [[ ${START} -ne 0 || ${END} -ne 100  ]]
	then
		log "Cleaning flag passed, but start/end specified. " 2
	fi


	if [[ -f ${INPUT}.bai ]]
	then
		rm -f ${INPUT}.bai
		log "Previous index file removed. " 3
	fi

	rm -rf ${TEMP_DIR}
	log "Previous tmp dir removed. " 3

	rm -rf ${LOG_DIR}
	log "Previous log dir removed. " 3

	rm -rf ${RESULTS_DIR}
	log "Previous results file removed" 3
fi


if [[ ! -d ${TEMP_DIR} ]]
then
	mkdir ${TEMP_DIR}
fi


if [[ ! -d ${TEMP_DIR}/SPLIT ]]
then
	mkdir ${TEMP_DIR}/SPLIT
fi


if [[ ! -d ${TEMP_DIR}/SHUF ]]
then
	mkdir ${TEMP_DIR}/SHUF
fi


if [[ ! -d ${RESULTS_DIR} ]]
then
	mkdir ${RESULTS_DIR}
fi




### MAIN SECTION
log $(printf '#%.0s' $(seq 1 $(($(tput cols)-35)) ) ) 3 2> /dev/null
log " " 3 


### Merge all input BAM parts
# note: assumes parts are in form ${OUTPUT}_part1.bam
CURR=$(( CURR + 1 ))

if [[ ${START} -le ${CURR} && ${END} -ge ${CURR} ]]
then
	NAME="01_mergeBAM"
	log ${NAME} 3
	(. ${CURR_SCRIPT_DIR}/${NAME}.sh ${NAME})
	testResult $?
fi



### Fix indexing errors by re-converting SAM to BAM
CURR=$(( CURR + 1 ))

if [[ ${START} -le ${CURR} && ${END} -ge ${CURR} ]]
then
	NAME="02_fixIndex"
	log ${NAME} 3
	(. ${CURR_SCRIPT_DIR}/${NAME}.sh ${NAME})
	testResult $?
fi



### Sort reads by name
CURR=$(( CURR + 1 ))

if [[ ${START} -le ${CURR} && ${END} -ge ${CURR} ]]
then
	NAME="03_sort"
	log ${NAME} 3
	(. ${CURR_SCRIPT_DIR}/${NAME}.sh ${NAME})
	testResult $?
fi



### Get the head and tail, to check for malformed reads
### Note, checking is manual
CURR=$(( CURR + 1 ))

if [[ ${START} -le ${CURR} && ${END} -ge ${CURR} ]]
then
	NAME="04_ends"
	log ${NAME} 3
	(. ${CURR_SCRIPT_DIR}/${NAME}.sh ${NAME})
	testResult $?
fi



### Strip out all alignment information
CURR=$(( CURR + 1 ))

if [[ ${START} -le ${CURR} && ${END} -ge ${CURR} ]]
then
	NAME="05_revert"
	log ${NAME} 3
	(. ${CURR_SCRIPT_DIR}/${NAME}.sh ${NAME})
	testResult $?
fi



### Detect reads which have the same name, are paired, are not secondary 
### alignment and are both first or second reads. This process takes
### a long time, so is best avoided if it's not needed. Run with the -p flag
CURR=$(( CURR + 1 ))

if [[ ${START} -le ${CURR} && ${END} -ge ${CURR} ]]
then
	NAME="06_badPairs"
	log ${NAME} 3
	(. ${CURR_SCRIPT_DIR}/${NAME}.sh ${NAME})
	testResult $?
fi



### Remove reads which have been flagged as bad 
CURR=$(( CURR + 1 ))

if [[ ${START} -le ${CURR} && ${END} -ge ${CURR} ]]
then
	NAME="07_removeReads"
	log ${NAME} 3
	(. ${CURR_SCRIPT_DIR}/${NAME}.sh ${NAME})
	testResult $?
fi



### Fix mate information
CURR=$(( CURR + 1 ))

if [[ ${START} -le ${CURR} && ${END} -ge ${CURR} ]]
then
	NAME="08_fixMate"
	log ${NAME} 3
	(. ${CURR_SCRIPT_DIR}/${NAME}.sh ${NAME})
	testResult $?
fi



### Add a read group
CURR=$(( CURR + 1 ))

if [[ ${START} -le ${CURR} && ${END} -ge ${CURR} ]]
then
	NAME="09_rg"
	log ${NAME} 3
	(. ${CURR_SCRIPT_DIR}/${NAME}.sh ${NAME})
	testResult $?
fi



### Validate and check for errors. There may be some left over, but most
### should have been removed by previous steps
CURR=$(( CURR + 1 ))

if [[ ${START} -le ${CURR} && ${END} -ge ${CURR} ]]
then
	NAME="10_validate"
	log ${NAME} 3
	(. ${CURR_SCRIPT_DIR}/${NAME}.sh ${NAME})
	testResult $?
fi



### Split the unmapped file into 10 Mbp chunks
CURR=$(( CURR + 1 ))

if [[ ${START} -le ${CURR} && ${END} -ge ${CURR} ]]
then
	NAME="11_split"
	log ${NAME} 3
	(. ${CURR_SCRIPT_DIR}/${NAME}.sh ${NAME})
	testResult $?
fi



### Shuffle the split files
CURR=$(( CURR + 1 ))

if [[ ${START} -le ${CURR} && ${END} -ge ${CURR} ]]
then
	NAME="12_shuf"
	log ${NAME} 3
	(. ${CURR_SCRIPT_DIR}/${NAME}.sh ${NAME})
	testResult $?
fi



### Merge the shuffled files
CURR=$(( CURR + 1 ))

if [[ ${START} -le ${CURR} && ${END} -ge ${CURR} ]]
then
	NAME="13_merge"
	log ${NAME} 3
	(. ${CURR_SCRIPT_DIR}/${NAME}.sh ${NAME})
	testResult $?
fi



### Convert BAM to FASTQ
CURR=$(( CURR + 1 ))

if [[ ${START} -le ${CURR} && ${END} -ge ${CURR} ]]
then
	NAME="14_convert"
	log ${NAME} 3
	(. ${CURR_SCRIPT_DIR}/${NAME}.sh ${NAME})
	testResult $?
fi


### Compress FASTQ files
CURR=$(( CURR + 1 ))

if [[ ${START} -le ${CURR} && ${END} -ge ${CURR} ]]
then
	NAME="15_compress"
	log ${NAME} 3
	(. ${CURR_SCRIPT_DIR}/${NAME}.sh ${NAME})
	testResult $?
fi


log $(printf '#%.0s' $(seq 1 $(($(tput cols)-35)) ) ) 3 2> /dev/null

# End of script
