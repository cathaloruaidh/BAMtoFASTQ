#! /bin/bash


samtools merge -f -@ ${THREADS} ${RESULTS_DIR}/${MERGED_INPUT_BAM} $(ls ${INPUT}_part*.bam) \
	2> >(tee ${LOG_DIR}/${OUTPUT}.${1}.err.txt >&2)

MERGE_RET=$?


if [[ ${MERGE_RET} -ne 0 ]]
then
	log "Error with merging BAM files. Exiting ... " 1
	exit 10
fi


log "File merged" 4

