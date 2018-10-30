#! /bin/bash

### Index the file and get the contigs


log "Creating index BAI file. " 3

echo ${READY_BAM} >> ${LOG_DIR}/${PREFIX}.${1}.stdout.txt 
echo ${READY_BAM} >> ${LOG_DIR}/${PREFIX}.${1}.stderr.txt

samtools index -@ ${THREADS} ${RESULTS_DIR}/${READY_BAM} \
	1>> ${LOG_DIR}/${PREFIX}.${1}.stdout.txt 2>> ${LOG_DIR}/${PREFIX}.${1}.stderr.txt
IDX_RET=$?

log "Index BAI file created." 4

echo >> ${LOG_DIR}/${PREFIX}.${1}.stdout.txt 
echo >> ${LOG_DIR}/${PREFIX}.${1}.stderr.txt


log "Creating contig text files. " 4
samtools idxstats ${RESULTS_DIR}/${READY_BAM} 2>> ${LOG_DIR}/${PREFIX}.${1}.stderr.txt | awk 'BEGIN{OFS="\t";} {if($3+$4 > 0 && $1 != "*") print $1,1,$2}' > ${TEMP_DIR}/all_contigs.bed

grep -E '_|[*]|EBV' ${TEMP_DIR}/all_contigs.bed > ${TEMP_DIR}/oth_contigs.bed
grep -v -E '_|[*]|EBV' ${TEMP_DIR}/all_contigs.bed > ${TEMP_DIR}/std_contigs.bed



if [[ ${IDX_RET} -ne 0 ]]
then
	exit 10
fi
