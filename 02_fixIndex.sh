#! /bin/bash



# Fix indexing errors by converting BAM -> SAM -> BAM
samtools view -h ${RESULTS_DIR}/${MERGED_INPUT_BAM} \
	| samtools view -@ ${THREADS} -h -Sb - > ${RESULTS_DIR}/${INDEX_INPUT_BAM} \
	2> >(tee ${LOG_DIR}/${OUTPUT}.${1}.view.err.txt >&2) 
BAM_RT=$?


# Index the file
samtools index ${RESULTS_DIR}/${MERGED_INPUT_BAM}


if [[ ${BAM_RT} -ne 0 ]]
then
	log "Error with reading merged BAM file. Exiting ... " 1
	exit 10
fi


log "Fixed index errors" 4


if [[ ${KEEP} -eq 0 ]]
then
	rm ${RESULTS_DIR}/${MERGED_INPUT_BAM}
fi
