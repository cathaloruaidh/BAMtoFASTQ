#! /bin/bash



# Get head of BAM file
samtools view ${RESULTS_DIR}/${SORTED_INPUT_BAM} | head -n 1000 > ${RESULTS_DIR}/${SORTED_INPUT_HEAD_TXT}
	2> >(tee ${LOG_DIR}/${OUTPUT}.${1}.head.err.txt >&2) 
HEAD_RT=$?
log "Head of file. " 4


# Get tail of BAM file
samtools view ${RESULTS_DIR}/${SORTED_INPUT_BAM} | tail -n 1000 > ${RESULTS_DIR}/${SORTED_INPUT_TAIL_TXT}
	2> >(tee ${LOG_DIR}/${OUTPUT}.${1}.tail.err.txt >&2) 
TAIL_RT=$?
log "Tail of file. " 4



rm -f ${RESULTS_DIR}/${SORTED_INPUT_REMOVE_TXT}
touch ${RESULTS_DIR}/${SORTED_INPUT_REMOVE_TXT}



if [[ ${HEAD_RT} -ne 0 || ${TAIL_RT} -ne 0 ]]
then
	log "Error with getting head or tail of file. Exiting ... " 1
	exit 10
fi


log "Outputted head and tail of BAM - manual inspection required." 4
