#! /bin/bash


for FILE in `ls ${TEMP_DIR}/SPLIT/${OUTPUT}*.split_*bam`
do
	echo ${FILE} >> ${LOG_DIR}/${OUTPUT}.${1}.stdout.txt 
	echo ${FILE} >> ${LOG_DIR}/${OUTPUT}.${1}.stderr.txt

	cat <(samtools view -H ${FILE}) <(samtools view -@4 ${FILE} | shuf 2>> ${LOG_DIR}/${OUTPUT}.${1}.stderr.txt) | samtools view -@ ${THREADS} -Sb -h > ${FILE%.bam}.${1}.bam \
		2>> ${LOG_DIR}/${OUTPUT}.${1}.stderr.txt
	
	SHUF_RET=$?

	if [ ${SHUF_RET} -ne 0 ]
	then
		log "Error with shuffling chromosome ${CHROMOSOME}. Exiting ... " 1
		exit 10
	fi


	log "${FILE} shuffled" 4


	echo >> ${LOG_DIR}/${OUTPUT}.${1}.stdout.txt 
	echo >> ${LOG_DIR}/${OUTPUT}.${1}.stderr.txt
done

mv ${TEMP_DIR}/SPLIT/${OUTPUT}*shuf.bam ${TEMP_DIR}/SHUF

if [[ ${KEEP} -eq 0 ]]
then
	rm -f ${TEMP_DIR}/SPLIT/${OUTPUT}*split*bam
	rm -f ${TEMP_DIR}/SPLIT/${OUTPUT}*split*sam
fi


