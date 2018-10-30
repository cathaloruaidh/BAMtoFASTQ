#! /bin/bash


java -Djava.io.temp=${TEMP_DIR} -jar /home/shared/tools/picard.jar CleanSam \
	I=${RESULTS_DIR}/${INDEX_INPUT_BAM} \
	O=${RESULTS_DIR}/${CLEAN_INPUT_BAM} \
	TMP_DIR=${TEMP_DIR} \
	> >(tee ${LOG_DIR}/${PREFIX}.${1}.out.txt) \
	2> >(tee ${LOG_DIR}/${PREFIX}.${1}.err.txt >&2)

CLEAN_RET=$?


if [[ ${CLEAN_RET} -ne 0 ]]
then
	log "Error with CleanSam. Exiting ... " 1
	exit 10
fi


log "CleanSam complete. " 4


if [[ ${KEEP} -eq 0 ]]
then
	rm -f ${RESULTS_DIR}/${INDEX_INPUT_BAM}
fi
