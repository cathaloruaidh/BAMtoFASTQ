#! /bin/bash


java -Djava.io.temp=${TEMP_DIR} -Xmx6G -Xms6G -jar /home/shared/tools/picard.jar CleanSam \
	I=${RESULTS_DIR}/${MERGED_BAM} \
	O=${RESULTS_DIR}/${CLEAN_BAM} \
	VALIDATION_STRINGENCY=SILENT \
	1> >(tee ${LOG_DIR}/${PREFIX}.${1}.out.txt) \
	2> >(tee ${LOG_DIR}/${PREFIX}.${1}.stderr.txt >&2)


CLEAN_RET=$?

if [ ${CLEAN_RET} -ne 0 ]
then
	log "Error with CleanSam. Exiting ... " 1
	exit 10
fi


if [[ ${KEEP} -eq 0 ]]
then
	rm -f ${RESULTS_DIR}/${PREFIX}.shuf.merged.clean.bam
fi
