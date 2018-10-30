#! /bin/bash




samtools merge -f -@ ${THREADS} ${RESULTS_DIR}/${MERGED_BAM} $(ls ${TEMP_DIR}/SHUF/${OUTPUT}*shuf.bam | shuf) \
	2> >(tee ${LOG_DIR}/${OUTPUT}.${1}.err.txt >&2)

MERGE_RET=$?

if [[ ${MERGE_RET} -ne 0 ]]
then
	log "Error with merge. Exiting ... " 1
	exit 10
fi


log "File merged" 4





if [[ ${KEEP} -eq 0 ]]
then
	rm -f ${TEMP_DIR}/SHUF/${OUTPUT}*bam
fi


