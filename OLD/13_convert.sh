#! /bin/bash


java -Djava.io.temp=${TEMP_DIR} -Xmx16G -Xms16G -jar ~/tools/picard.jar SamToFastq \
	I=${RESULTS_DIR}/${CLEAN_BAM} \
	FASTQ=${RESULTS_DIR}/${PREFIX}_R1.fastq \
	SECOND_END_FASTQ=${RESULTS_DIR}/${PREFIX}_R2.fastq \
	UNPAIRED_FASTQ=${RESULTS_DIR}/${PREFIX}_unpaired.fastq \
	INCLUDE_NON_PF_READS=true \
	INCLUDE_NON_PRIMARY_ALIGNMENTS=true \
	VALIDATION_STRINGENCY=SILENT \
	1> >(tee ${LOG_DIR}/${PREFIX}.${1}.stdout.txt) \
	2> >(tee ${LOG_DIR}/${PREFIX}.${1}.stderr.txt >&2)

CONVERT_RET=$?

if [ ${CONVERT_RET} -ne 0 ]
then
	log "Error with SamToFastq. Exiting ... " 1
	exit 10
fi
