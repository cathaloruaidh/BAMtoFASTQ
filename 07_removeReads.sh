#! /bin/bash



# send an email to inform about removal of bad reads
if [ -z "$STY" ]
then
	LOC=$(tty)
else
	LOC=$(tty)", through screen $STY"
fi


echo -e "You need to examine the head and tail of the BAM file in ${RESULTS_DIR}. Check the following files: \n\n \
${RESULTS_DIR}/${SORTED_INPUT_HEAD_TXT} \n \
${RESULTS_DIR}/${SORTED_INPUT_TAIL_TXT} \n\n \
Add any reads to be removed to ${RESULTS_DIR}/${SORTED_INPUT_REMOVE_TXT} \n\n \
Run from: $LOC\n\n" | mailx -s "Job on Vishnu: Action required" caormond@tcd.ie


log "Note! Check ${RESULTS_DIR}/${SORTED_INPUT_HEAD_TXT} and ${RESULTS_DIR}/${SORTED_INPUT_TAIL_TXT}, and add any files to be removed to ${RESULTS_DIR}/${SORTED_INPUT_REMOVE_TXT}" 2
read -e -t2
read -n 1 -r -p "Press any key to continue once bad reads have been added to file."


if [ ! -f ${RESULTS_DIR}/${REVERT_INPUT_BAD_REMOVE} ]
then
	touch ${RESULTS_DIR}/${REVERT_INPUT_BAD_REMOVE}
fi


cat ${RESULTS_DIR}/${REVERT_INPUT_BAD_REMOVE} ${RESULTS_DIR}/${SORTED_INPUT_REMOVE_TXT} > ${RESULTS_DIR}/${REMOVE_READS_TXT}


if [[ $(wc -l < ${RESULTS_DIR}/${REMOVE_READS_TXT}) -ne 0 ]]
then
	# Remove all badly-paired reads
	java -Djava.io.temp=${TEMP_DIR} -Xmx6G -Xms6G -jar ${PICARD_FILE} FilterSamReads \
		I=${RESULTS_DIR}/${REVERT_INPUT_BAM} \
		O=${RESULTS_DIR}/${READS_INPUT_BAM} \
		FILTER=excludeReadList \
		READ_LIST_FILE=${RESULTS_DIR}/${REMOVE_READS_TXT} \
		SORT_ORDER=queryname \
		TMP_DIR=${TEMP_DIR} \
		2> >(tee ${LOG_DIR}/${OUTPUT}.${1}.filter.err.txt >&2)

	REM_RET=$?

	if [[ ${REM_RET} -ne 0 ]]
	then
		log "Error with FilterSamReads. Exiting ... " 1
		exit 10
	fi

	log "Removed all poorly paired reads (duplicate read name per alignment and pair position)" 4
else
	cp ${RESULTS_DIR}/${REVERT_INPUT_BAM} ${RESULTS_DIR}/${READS_INPUT_BAM}
fi


if [[ ${KEEP} -eq 0 ]]
then
	rm ${RESULTS_DIR}/${REVERT_INPUT_BAM}
fi
