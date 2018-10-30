#! /bin/bash


for CHROMOSOME in `cut -f 1 ${TEMP_DIR}/std_contigs.bed` 
do
	echo ${CHROMOSOME} >> ${LOG_DIR}/${PREFIX}.${1}.stderr.txt

	samtools view -h ${RESULTS_DIR}/${READY_BAM} ${CHROMOSOME} | samtools view -@ ${THREADS} -Sb -h - > ${TEMP_DIR}/CHR/${PREFIX}.${CHROMOSOME}.bam \
		2>> ${LOG_DIR}/${PREFIX}.${1}.stderr.txt
	SPLIT_RET=$?

	if [[ ${SPLIT_RET} -ne 0 ]]
	then
		log "Error with splitting chromosome ${CHROMOSOME}. Exiting ... " 1
		exit 10
	fi

	log "Split contig: ${CHROMOSOME}" 4

	echo >> ${LOG_DIR}/${PREFIX}.${1}.stderr.txt
done



## get other contigs 
echo "Other" >> ${LOG_DIR}/${PREFIX}.${1}.stderr.txt

samtools view -h ${RESULTS_DIR}/${READY_BAM} -L ${TEMP_DIR}/oth_contigs.bed | samtools view -@ ${THREADS} -Sb -h - > ${TEMP_DIR}/CHR/${PREFIX}.other.bam \
	2>> ${LOG_DIR}/${PREFIX}.${1}.stderr.txt
SPLIT_RET=$?


if [[ ${SPLIT_RET} -ne 0 ]]
then
	log "Error with splitting chromosome ${CHROMOSOME}. Exiting ... " 1
	exit 10
fi


log "Split contig: others" 4

echo >> ${LOG_DIR}/${PREFIX}.${1}.stderr.txt



## get unmapped reads
echo "Unmapped" >> ${LOG_DIR}/${PREFIX}.chr.stderr.txt

samtools view -h ${RESULTS_DIR}/${READY_BAM} '*' | samtools view -@ ${THREADS} -Sb -h - > ${TEMP_DIR}/CHR/${PREFIX}.unmapped.bam \
	2>> ${LOG_DIR}/${PREFIX}.${1}.stderr.txt
SPLIT_RET=$?

if [[ ${SPLIT_RET} -ne 0 ]]
then
	log "Error with splitting chromosome ${CHROMOSOME}. Exiting ... " 1
	exit 10
fi

log "Split contig: unmapped" 4

echo >> ${LOG_DIR}/${PREFIX}.${1}.stderr.txt


