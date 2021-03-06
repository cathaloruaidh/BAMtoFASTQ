#! /bin/bash


java -Xmx6G -Xms6G -Djava.io.temp=${TEMP_DIR} -jar ${PICARD_FILE} ValidateSamFile \
	I=${RESULTS_DIR}/${RG_INPUT_BAM} \
	O=${RESULTS_DIR}/${VALIDATE_INPUT_TXT} \
	MODE=SUMMARY \
	MAX_OUTPUT=null \
	TMP_DIR=${TEMP_DIR} \
	> >(tee ${LOG_DIR}/${OUTPUT}.${1}.validate.out.txt) \
	2> >(tee ${LOG_DIR}/${OUTPUT}.${1}.validate.err.txt >&2)

VAL_RET=$?


if [[ ${FIX_RET} -ne 0 ]]
then
	log "Error with ValidateSamFile. Exiting ... " 1
	exit 10
fi


log "ValidateSamFile complete. " 4


cp ${RESULTS_DIR}/${RG_INPUT_BAM} ${RESULTS_DIR}/${READY_BAM}
samtools index ${RESULTS_DIR}/${READY_BAM}
