#! /bin/bash


log "Read 1" 3
gzip ${RESULTS_DIR}/${OUTPUT}_R1.fastq \
	1> >(tee ${LOG_DIR}/${OUTPUT}.${1}.1.out.txt)  \
	2> >(tee ${LOG_DIR}/${OUTPUT}.${1}.1.err.txt >&2)
R1_RET=$?

log "Read 2" 3
gzip ${RESULTS_DIR}/${OUTPUT}_R2.fastq \
	1> >(tee ${LOG_DIR}/${OUTPUT}.${1}.1.out.txt)  \
	2> >(tee ${LOG_DIR}/${OUTPUT}.${1}.1.err.txt >&2)
R2_RET=$?

log "Unpaired Reads" 3
gzip ${RESULTS_DIR}/${OUTPUT}_unpaired.fastq \
	1> >(tee ${LOG_DIR}/${OUTPUT}.${1}.1.out.txt)  \
	2> >(tee ${LOG_DIR}/${OUTPUT}.${1}.1.err.txt >&2)
RU_RET=$?


if [ ${R1_RET} -ne 0 ]
then
	log "Error with gzip on read 1. Exiting ... " 1
	exit 1
fi


if [ ${R2_RET} -ne 0 ]
then
	log "Error with gzip on read 2. Exiting ... " 1
	exit 1
fi


if [ ${RU_RET} -ne 0 ]
then
	log "Error with gzip on unpaired read. Exiting ... " 1
	exit 1
fi
