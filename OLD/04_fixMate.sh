#! /bin/bash


java -Djava.io.temp=${TEMP_DIR} -Xmx12G -Xms12G -jar /home/shared/tools/picard.jar FixMateInformation \
	I=${RESULTS_DIR}/${CLEAN_INPUT_BAM} \
	O=${RESULTS_DIR}/${MATE_INPUT_BAM} \
	TMP_DIR=${TEMP_DIR} \
	VALIDATION_STRINGENCY=LENIENT \
	MAX_RECORDS_IN_RAM=800000 \
	> >(tee ${LOG_DIR}/${PREFIX}.${1}.out.txt) \
	2> >(tee ${LOG_DIR}/${PREFIX}.${1}.err.txt >&2)

FIX_RET=$?


if [[ ${FIX_RET} -ne 0 ]]
then
	log "Error with FixMateInformation. Exiting ... " 1
	exit 10
fi


log "FixMateInformation complete. " 4


if [[ ${KEEP} -eq 0 ]]
then
	rm -f ${RESULTS_DIR}/${CLEAN_INPUT_BAM}
fi
