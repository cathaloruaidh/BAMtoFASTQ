#! /bin/bash


java -Xmx6G -Xms6G -Djava.io.temp=${TEMP_DIR} -jar ${PICARD_FILE} RevertSam \
	I=${RESULTS_DIR}/${SORTED_INPUT_BAM} \
	O=${RESULTS_DIR}/${REVERT_INPUT_BAM} \
	VALIDATION_STRINGENCY=SILENT \
	TMP_DIR=${TEMP_DIR} \
	MAX_RECORDS_IN_RAM=1000000 \
	2> >(tee ${LOG_DIR}/${OUTPUT}.${1}.err.txt >&2)

REV_RET=$?


if [[ ${REV_RET} -ne 0 ]]
then
	log "Error with RevertSam. Exiting ... " 1
	exit 10
fi


log "RevertSam complete. " 4



if [[ ${KEEP} -eq 0 ]]
then
	rm ${RESULTS_DIR}/${SORTED_INPUT_BAM}
fi
