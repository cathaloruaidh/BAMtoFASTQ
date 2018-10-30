#! /bin/bash


samtools idxstats ${RESULTS_DIR}/${READY_BAM} 2>> ${LOG_DIR}/${OUTPUT}.${1}.stats.err.txt | awk 'BEGIN{OFS="\t";} {if($3+$4 > 0 && $1 != "*") print $1,$2}' > ${TEMP_DIR}/all_contigs.bed



samtools view -@ ${THREADS} ${RESULTS_DIR}/${READY_BAM} | split -d -l 10000000 --verbose --additional-suffix=.sam - ${TEMP_DIR}/SPLIT/${READY_BAM%.bam}.split_ \
	2> >(tee ${LOG_DIR}/${OUTPUT}.${1}.err.txt >&2)

SPLIT_RET=$?
if [[ ${SPLIT_RET} -ne 0 ]]
then
	log "Error with splitting file. Exiting ... " 1
	exit 10
fi

log "File split." 4



for FILE in `ls ${TEMP_DIR}/SPLIT/${OUTPUT}*.split_*.sam`
do
	echo ${FILE} >> ${LOG_DIR}/${OUTPUT}.${1}.view.stdout.txt 
	echo ${FILE} >> ${LOG_DIR}/${OUTPUT}.${1}.view.stderr.txt

	#cat <(samtools view -H ${RESULTS_DIR}/${READY_BAM}) <(samtools view ${FILE}) | samtools view -@ ${THREADS} -Sb -h - > ${FILE%sam}bam 2>> ${LOG_DIR}/${OUTPUT}.${1}.view.stderr.txt
	#cat <(samtools view -H ${RESULTS_DIR}/${READY_BAM}) ${FILE} | samtools view -@ ${THREADS} -Sb -h - > ${FILE%sam}bam 2>> ${LOG_DIR}/${OUTPUT}.${1}.view.stderr.txt
	cat <(samtools view -H ${RESULTS_DIR}/${READY_BAM}) <(samtools view -t ${TEMP_DIR}/all_contigs.bed ${FILE}) | samtools view -@ ${THREADS} -Sb -h - > ${FILE%sam}bam 2>> ${LOG_DIR}/${OUTPUT}.${1}.view.err.txt

	VIEW_RET=$?


	if [ ${VIEW_RET} -ne 0 ]
	then
		log "Error with creating BAM file for ${FILE}. Exiting ... " 1
		exit 10
	fi

	rm -f ${FILE}
	log "Created: ${FILE%sam}bam" 4
	
	echo >> ${LOG_DIR}/${OUTPUT}.${1}.view.stdout.txt 
	echo >> ${LOG_DIR}/${OUTPUT}.${1}.view.stderr.txt
done

log " " 4


if [[ ${KEEP} -eq 0 ]]
then
	rm -f ${TEMP_DIR}/CHR/${OUTPUT}*bam
fi

