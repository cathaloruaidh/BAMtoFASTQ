#! /bin/bash



# Fix errors with pair membership for reads with the same read name
# Sort by name and run "resolvepair.py", convert to BAM
#samtools sort -@${THREADS} -T ${TEMP_DIR} -n ${RESULTS_DIR}/${MERGED_INPUT_BAM} | ${SCRIPT_DIR}/resolvepair.py | samtools view -@ ${THREADS} -h -Sb - > ${RESULTS_DIR}/${MERGED_INPUT_BAM}.resolvepair
#RES_RT=$?
#log "Pairs resolved. " 4


# Sort the File by read name
samtools sort -n -@ ${THREADS} -T ${TEMP_DIR} ${RESULTS_DIR}/${MERGED_INPUT_BAM} -o ${RESULTS_DIR}/${MERGED_INPUT_BAM}.sorted \
	> >(tee ${LOG_DIR}/${PREFIX}.${1}.sort.out.txt) \
	2> >(tee ${LOG_DIR}/${PREFIX}.${1}.sort.err.txt >&2) 
SRT_RT=$?
log "Sorted file. " 4


# Fix indexing errors by converting BAM -> SAM -> BAM
samtools view -h ${RESULTS_DIR}/${MERGED_INPUT_BAM}.sorted \
	| samtools view -@ ${THREADS} -h -Sb - > ${RESULTS_DIR}/${INDEX_INPUT_BAM} \
	> >(tee ${LOG_DIR}/${PREFIX}.${1}.view.out.txt) \
	2> >(tee ${LOG_DIR}/${PREFIX}.${1}.view.err.txt >&2) 
BAM_RT=$?


# Index the file
samtools index ${RESULTS_DIR}/${INDEX_INPUT_BAM}


#if [[ ${BAM_RT} -ne 0 || ${RES_RT} -ne 0 || ${SRT_RT} -ne 0 ]]
if [[ ${BAM_RT} -ne 0 || ${SRT_RT} -ne 0 ]]
then
	log "Error with reading merged BAM file. Exiting ... " 1
	exit 10
fi


log "Fixed index errors" 4


if [[ ${KEEP} -eq 0 ]]
then
	rm ${RESULTS_DIR}/${MERGED_INPUT_BAM}
	rm ${RESULTS_DIR}/${MERGED_INPUT_BAM}.sorted
fi
