#! /bin/bash


samtools merge -f -@ ${THREADS} ${RESULTS_DIR}/${MERGED_INPUT_BAM} $(ls ${PWD}/${PREFIX}_part*.bam) \
	> >(tee ${LOG_DIR}/${PREFIX}.${1}.out.txt) \
	2> >(tee ${LOG_DIR}/${PREFIX}.${1}.err.txt >&2)

MERGE_RET=$?


if [[ ${MERGE_RET} -ne 0 ]]
then
	log "Error with merging BAM files. Exiting ... " 1
	exit 10
fi


log "File merged" 4

