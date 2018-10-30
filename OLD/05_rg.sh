#! /bin/bash


java -Djava.io.temp=${TEMP_DIR} -Xmx6G -Xms6G -jar /home/shared/tools/picard.jar AddOrReplaceReadGroups \
	I=${RESULTS_DIR}/${MATE_INPUT_BAM} \
	O=${RESULTS_DIR}/${RG_INPUT_BAM} \
	RGID=${PREFIX} \
	RGLB=lib \
	RGPL=Illumina \
	RGPU=${PREFIX} \
	RGSM=${PREFIX} \
	VALIDATION_STRINGENCY=LENIENT \
	TMP_DIR=${TEMP_DIR} \
	> >(tee ${LOG_DIR}/${PREFIX}.${1}.out.txt) \
	2> >(tee ${LOG_DIR}/${PREFIX}.${1}.err.txt >&2)

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


