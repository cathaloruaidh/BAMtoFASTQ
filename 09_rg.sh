#! /bin/bash


java -Xmx6G -Xms6G -Djava.io.temp=${TEMP_DIR} -jar ${PICARD_FILE} AddOrReplaceReadGroups \
	I=${RESULTS_DIR}/${MATE_INPUT_BAM} \
	O=${RESULTS_DIR}/${RG_INPUT_BAM} \
	RGID=${OUTPUT} \
	RGLB=lib \
	RGPL=Illumina \
	RGPU=${OUTPUT} \
	RGSM=${OUTPUT} \
	VALIDATION_STRINGENCY=LENIENT \
	TMP_DIR=${TEMP_DIR} \
	> >(tee ${LOG_DIR}/${OUTPUT}.${1}.out.txt) \
	2> >(tee ${LOG_DIR}/${OUTPUT}.${1}.err.txt >&2)

RG_RET=$?


if [[ ${RG_RET} -ne 0 ]]
then
	log "Error with creating read group. Exiting ... " 1
	exit 10
fi


log "Created read group" 4


if [[ ${KEEP} -eq 0 ]]
then
	rm -f ${RESULTS_DIR}/${MATE_INPUT_BAM}
fi


