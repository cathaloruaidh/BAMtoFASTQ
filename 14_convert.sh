#! /bin/bash


java -Xmx16G -Xms16G -Djava.io.temp=${TEMP_DIR} -jar ${PICARD_FILE} SamToFastq \
	I=${RESULTS_DIR}/${MERGED_BAM} \
	FASTQ=${RESULTS_DIR}/${OUTPUT}_R1.fastq \
	SECOND_END_FASTQ=${RESULTS_DIR}/${OUTPUT}_R2.fastq \
	UNPAIRED_FASTQ=${RESULTS_DIR}/${OUTPUT}_unpaired.fastq \
	INCLUDE_NON_PF_READS=true \
	INCLUDE_NON_PRIMARY_ALIGNMENTS=true \
	VALIDATION_STRINGENCY=SILENT \
	TMP_DIR=${TEMP_DIR} \
	2> >(tee ${LOG_DIR}/${OUTPUT}.${1}.stderr.txt >&2)

CONVERT_RET=$?

if [ ${CONVERT_RET} -ne 0 ]
then
	log "Error with SamToFastq. Exiting ... " 1
	exit 10
fi
