#! /bin/bash


### Pipefail
set -o pipefail



### Argument processing
INPUT=""
VERBOSE=3
CLEAN=0
KEEP=0
START=0
END=100
CURR=0
PAUSE=0
THREADS=1

while getopts 'f:v:cks:e:pt:' flag ; do
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
			PAUSE=1
			;;

		t)
			THREADS=${OPTARG}
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



PREFIX=${INPUT%.bam}
MERGED_INPUT_BAM=${PREFIX}_merged.bam
INDEX_INPUT_BAM=${PREFIX}_merged.ind.bam
CLEAN_INPUT_BAM=${PREFIX}_merged.ind.clean.bam
MATE_INPUT_BAM=${PREFIX}_merged.ind.clean.mate.bam
RG_INPUT_BAM=${PREFIX}_merged.ind.clean.mate.rg.bam
VALIDATE_INPUT_TXT=${PREFIX}_merged.ind.clean.mate.rg.validate.txt

READY_BAM=${PREFIX}.bam

MERGED_BAM=${PREFIX}_shuffled.bam
SORTED_BAM=${PREFIX}_shuffled.sort.bam
CLEAN_BAM=${PREFIX}_shuffled.sort.clean.bam


TEMP_DIR="${PWD}/tmp"
LOG_DIR="${PWD}/log"
RESULTS_DIR="${PWD}/results"

SCRIPT_DIR="/home/shared/scripts/"
CURR_SCRIPT_DIR="${SCRIPT_DIR}BAMtoFASTQ"


MAIN_LOG_FILE=${LOG_DIR}/main.log


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


ulimit -n unlimited
log "Note: ulimit is $(ulimit -n), max limit is $(ulimit -Hn)" 4


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
	mkdir ${TEMP_DIR}/CHR
	mkdir ${TEMP_DIR}/SPLIT
	mkdir ${TEMP_DIR}/SHUF
fi


if [[ ! -d ${LOG_DIR} ]]
then
	mkdir ${LOG_DIR}
fi


if [[ ! -d ${RESULTS_DIR} ]]
then
	mkdir ${RESULTS_DIR}
fi




### MAIN SECTION
log $(printf '#%.0s' $(seq 1 $(($(tput cols)-35)) ) ) 3 2> /dev/null
log " " 3 


### Merge all input BAM parts
# note: assumes parts are in form ${PREFIX}_part1.bam
CURR=$(( CURR + 1 ))

if [[ ${START} -le ${CURR} && ${END} -ge ${CURR} ]]
then
	NAME="01_mergeBAM"
	log ${NAME} 3
	(. ${CURR_SCRIPT_DIR}/${NAME}.sh ${NAME})
	testResult $?
	
	if [[ ${PAUSE} -ne 0 ]]
	then
		read -n1 -r -p "Press any key to continue" key
	fi
fi



### Add read group to merged input BAM
CURR=$(( CURR + 1 ))

if [[ ${START} -le ${CURR} && ${END} -ge ${CURR} ]]
then
	NAME="02_fixIndex"
	log ${NAME} 3
	(. ${CURR_SCRIPT_DIR}/${NAME}.sh ${NAME})
	testResult $?
	
	if [[ ${PAUSE} -ne 0 ]]
	then
		read -n1 -r -p "Press any key to continue" key
	fi
fi



### Fix indexing errors by re-converting SAM to BAM
CURR=$(( CURR + 1 ))

if [[ ${START} -le ${CURR} && ${END} -ge ${CURR} ]]
then
	NAME="03_clean"
	log ${NAME} 3
	(. ${CURR_SCRIPT_DIR}/${NAME}.sh ${NAME})
	testResult $?
	
	if [[ ${PAUSE} -ne 0 ]]
	then
		read -n1 -r -p "Press any key to continue" key
	fi
fi



### Fix clipping errors with CleanSam
CURR=$(( CURR + 1 ))

if [[ ${START} -le ${CURR} && ${END} -ge ${CURR} ]]
then
	NAME="04_fixMate"
	log ${NAME} 3
	(. ${CURR_SCRIPT_DIR}/${NAME}.sh ${NAME})
	testResult $?
	
	if [[ ${PAUSE} -ne 0 ]]
	then
		read -n1 -r -p "Press any key to continue" key
	fi
fi



### Fix mate errors
CURR=$(( CURR + 1 ))

if [[ ${START} -le ${CURR} && ${END} -ge ${CURR} ]]
then
	NAME="05_rg"
	log ${NAME} 3
	(. ${CURR_SCRIPT_DIR}/${NAME}.sh ${NAME})
	testResult $?
	
	if [[ ${PAUSE} -ne 0 ]]
	then
		read -n1 -r -p "Press any key to continue" key
	fi
fi



### Validate the final, cleaned BAM
CURR=$(( CURR + 1 ))

if [[ ${START} -le ${CURR} && ${END} -ge ${CURR} ]]
then
	NAME="06_validate"
	log ${NAME} 3
	(. ${CURR_SCRIPT_DIR}/${NAME}.sh ${NAME})
	testResult $?
	
	if [[ ${PAUSE} -ne 0 ]]
	then
		read -n1 -r -p "Press any key to continue" key
	fi
fi



### Index the file and get the contigs
CURR=$(( CURR + 1 ))

if [[ ${START} -le ${CURR} && ${END} -ge ${CURR} ]]
then
	NAME="07_index"
	log ${NAME} 3
	(. ${CURR_SCRIPT_DIR}/${NAME}.sh ${NAME})
	testResult $?
	
	if [[ ${PAUSE} -ne 0 ]]
	then
		read -n1 -r -p "Press any key to continue" key
	fi
fi



### Split by chromosome
CURR=$(( CURR + 1 ))

if [[ ${START} -le ${CURR} && ${END} -ge ${CURR} ]]
then
	NAME="08_chr"
	log ${NAME} 3
	(. ${CURR_SCRIPT_DIR}/${NAME}.sh ${NAME})
	testResult $?
	
	if [[ ${PAUSE} -ne 0 ]]
	then
		read -n1 -r -p "Press any key to continue" key
	fi
fi



### Split by 10Mbp, and delete CHR BAM
CURR=$(( CURR + 1 ))

if [[ ${START} -le ${CURR} && ${END} -ge ${CURR} ]]
then
	NAME="09_split"
	log ${NAME} 3
	(. ${CURR_SCRIPT_DIR}/${NAME}.sh ${NAME})
	testResult $?
	
	if [[ ${PAUSE} -ne 0 ]]
	then
		read -n1 -r -p "Press any key to continue" key
	fi
fi



### Shuffle individual files
CURR=$(( CURR + 1 ))

if [[ ${START} -le ${CURR} && ${END} -ge ${CURR} ]]
then
	NAME="10_shuf"
	log ${NAME} 3
	(. ${CURR_SCRIPT_DIR}/${NAME}.sh ${NAME})
	testResult $?
	
	if [[ ${PAUSE} -ne 0 ]]
	then
		read -n1 -r -p "Press any key to continue" key
	fi
fi



### Merge shuffled files
CURR=$(( CURR + 1 ))

if [[ ${START} -le ${CURR} && ${END} -ge ${CURR} ]]
then
	NAME="11_merge"
	log ${NAME} 3
	(. ${CURR_SCRIPT_DIR}/${NAME}.sh ${NAME})
	testResult $?
	
	if [[ ${PAUSE} -ne 0 ]]
	then
		read -n1 -r -p "Press any key to continue" key
	fi
fi



### Set MAPQ to 0 if the read is unmapped
CURR=$(( CURR + 1 ))

if [[ ${START} -le ${CURR} && ${END} -ge ${CURR} ]]
then
	NAME="12_clean"
	log ${NAME} 3
	(. ${CURR_SCRIPT_DIR}/${NAME}.sh ${NAME})
	testResult $?
	
	if [[ ${PAUSE} -ne 0 ]]
	then
		read -n1 -r -p "Press any key to continue" key
	fi
fi



### Convert BAM to FASTQ
CURR=$(( CURR + 1 ))

if [[ ${START} -le ${CURR} && ${END} -ge ${CURR} ]]
then
	NAME="13_convert"
	log ${NAME} 3
	(. ${CURR_SCRIPT_DIR}/${NAME}.sh ${NAME})
	testResult $?
	
	if [[ ${PAUSE} -ne 0 ]]
	then
		read -n1 -r -p "Press any key to continue" key
	fi
fi



### Compress FASTQ files
CURR=$(( CURR + 1 ))

if [[ ${START} -le ${CURR} && ${END} -ge ${CURR} ]]
then
	NAME="14_compress"
	log ${NAME} 3
	(. ${CURR_SCRIPT_DIR}/${NAME}.sh ${NAME})
	testResult $?
	
	if [[ ${PAUSE} -ne 0 ]]
	then
		read -n1 -r -p "Press any key to continue" key
	fi
fi




log $(printf '#%.0s' $(seq 1 $(($(tput cols)-35)) ) ) 3 2> /dev/null

# End of script
