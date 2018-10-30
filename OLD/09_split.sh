#! /bin/bash


for CHR_BAM in `ls ${TEMP_DIR}/CHR/${PREFIX}.*bam` 
do
	echo ${CHR_BAM} >> ${LOG_DIR}/${PREFIX}.${1}.stdout.txt 
	echo ${CHR_BAM} >> ${LOG_DIR}/${PREFIX}.${1}.stderr.txt

	samtools view -@ ${THREADS} ${CHR_BAM} | split -d -l 10000000 --verbose --additional-suffix=.sam - ${CHR_BAM%.bam}.split_ \
		1>> ${LOG_DIR}/${PREFIX}.${1}.stdout.txt 2>> ${LOG_DIR}/${PREFIX}.${1}.stderr.txt

	SPLIT_RET=$?
	if [[ ${SPLIT_RET} -ne 0 ]]
	then
		log "Error with splitting chromosome ${CHROMOSOME}. Exiting ... " 1
		exit 10
	fi

	mv ${TEMP_DIR}/CHR/${PREFIX}.*.sam ${TEMP_DIR}/SPLIT
	log "Split: ${CHR_BAM}" 4

	echo >> ${LOG_DIR}/${PREFIX}.${1}.stdout.txt 
	echo >> ${LOG_DIR}/${PREFIX}.${1}.stderr.txt
done


for FILE in `ls ${TEMP_DIR}/SPLIT/${PREFIX}.*.split_*.sam`
do
	echo ${FILE} >> ${LOG_DIR}/${PREFIX}.${1}.view.stdout.txt 
	echo ${FILE} >> ${LOG_DIR}/${PREFIX}.${1}.view.stderr.txt

	cat <(samtools view -H ${RESULTS_DIR}/${READY_BAM}) <(samtools view -t <( awk '{print $1,$3}' OFS="\t" ${TEMP_DIR}/all_contigs.bed) ${FILE}) | samtools view -@ ${THREADS} -Sb -h - > ${FILE%sam}bam 2>> ${LOG_DIR}/${PREFIX}.${1}.view.stderr.txt

	VIEW_RET=$?


	if [ ${VIEW_RET} -ne 0 ]
	then
		log "Error with creating BAM file for ${FILE}. Exiting ... " 1
		exit 10
	fi

	log "Created: ${FILE%sam}bam" 4
	
	echo >> ${LOG_DIR}/${PREFIX}.${1}.view.stdout.txt 
	echo >> ${LOG_DIR}/${PREFIX}.${1}.view.stderr.txt
done

log " " 4


if [[ ${KEEP} -eq 0 ]]
then
	rm -f ${TEMP_DIR}/CHR/${PREFIX}*bam
fi

