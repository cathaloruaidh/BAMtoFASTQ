#! /bin/bash



# Sort the File by read name
java -Xmx6G -Xms6G -Djava.io.temp=${TEMP_DIR} -jar ${PICARD_FILE} SortSam \
	I=${RESULTS_DIR}/${INDEX_INPUT_BAM} \
	O=${RESULTS_DIR}/${SORTED_INPUT_BAM} \
	SORT_ORDER=queryname \
	VALIDATION_STRINGENCY=SILENT \
	MAX_RECORDS_IN_RAM=1000000 \
	TMP_DIR=${TEMP_DIR} \
	2> >(tee ${LOG_DIR}/${OUTPUT}.${1}.err.txt >&2) 
SRT_RT=$?
log "Sorted file. " 4


# Index the file
samtools index ${RESULTS_DIR}/${SORTED_INPUT_BAM}


if [[ ${SRT_RT} -ne 0 ]]
then
	log "Error with sorting BAM file by read name. Exiting ... " 1
	exit 10
fi


log "Sorted BAM file by read name" 4


if [[ ${KEEP} -eq 0 ]]
then
	rm ${RESULTS_DIR}/${INDEX_INPUT_BAM}
fi
