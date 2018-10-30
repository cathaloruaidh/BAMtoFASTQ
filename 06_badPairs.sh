#! /bin/bash

rm -f ${RESULTS_DIR}/${REVERT_INPUT_BAD_FIRST}
touch ${RESULTS_DIR}/${REVERT_INPUT_BAD_FIRST}
rm -f ${RESULTS_DIR}/${REVERT_INPUT_BAD_SECOND}
touch ${RESULTS_DIR}/${REVERT_INPUT_BAD_SECOND}


if [[ ${PAIRS} -ne 0 ]]
then
	# Get first read of all paired, primary-mapped reads. 
	# Extract all non-unique reads names, to be removed 
	samtools view -@ ${THREADS} -f 77 -F 2048 ${RESULTS_DIR}/${REVERT_INPUT_BAM} \
		| cut -f 1 | sort -T ${TEMP_DIR} | uniq -D > ${RESULTS_DIR}/${REVERT_INPUT_BAD_FIRST} \
		2> >(tee ${LOG_DIR}/${OUTPUT}.${1}.first.txt >&2)

	FIR_RET=$?

	if [[ ${FIR_RET} -ne 0 ]]
	then
		log "Error with getting bad first reads. Exiting ... " 1
		exit 10
	fi
	log "Found $( wc -l ${RESULTS_DIR}/${REVERT_INPUT_BAD_FIRST} | cut -f 1 -d ' ') first reads to be removed." 4



	# Get second read of all paired, primary-mapped reads. 
	# Extract all non-unique reads names, to be removed 
	samtools view -@ ${THREADS} -f 141 -F 2048 ${RESULTS_DIR}/${REVERT_INPUT_BAM} \
		| cut -f 1 | sort -T ${TEMP_DIR} | uniq -D > ${RESULTS_DIR}/${REVERT_INPUT_BAD_SECOND} \
		2> >(tee ${LOG_DIR}/${OUTPUT}.${1}.second.txt >&2)

	SEC_RET=$?

	if [[ ${SEC_RET} -ne 0 ]]
	then
		log "Error with bad second reads. Exiting ... " 1
		exit 10
	fi
	log "Found $( wc -l ${RESULTS_DIR}/${REVERT_INPUT_BAD_SECOND} | cut -f 1 -d ' ') second reads to be removed." 4



	# Combine all reads for removal 
	cat ${RESULTS_DIR}/${REVERT_INPUT_BAD_FIRST} ${RESULTS_DIR}/${REVERT_INPUT_BAD_SECOND} \
		| sort -T ${TEMP_DIR} | uniq > ${RESULTS_DIR}/${REVERT_INPUT_BAD_REMOVE}

	COM_RET=$?

	if [[ ${COM_RET} -ne 0 ]]
	then
		log "Error with combining first and second reads. Exiting ... " 1
		exit 10
	fi
	log "Selected all poorly paired reads (duplicate read name per alignment and pair position)" 4
fi

