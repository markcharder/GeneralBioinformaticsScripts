#!/usr/bin/env bash

# 23.02.16
# Script used for reformatting bam files (converted from sam derived from INRA genomic Illumina seq alignments using bwa. This is for the genome publication, calling SNPs etc.

MY_PICARD_BIN=/home/mark/bioinformatics/src/picard/dist/picard.jar
MY_WORKING_DIRECTORY=/media/rdrive/CCDM_Program_6-DENTOM-SE01498/samtools/

cd ${MY_WORKING_DIRECTORY}

for file in *.bam; do
    
    if [[ $file =~ U[1-2] ]]; then
        stub=($(echo "$file" | cut -d "." -f 1,2))
    else
        stub=($(echo "$file" | cut -d "." -f 1))
    fi

    java -jar ${MY_PICARD_BIN} AddOrReplaceReadGroups \
        INPUT=$file \
        OUTPUT=${stub}.add_readgroup.bam \
        RGID=INRA_RUN \
        RGLB=INRA_LIBRARY \
        RGPL=ILLUMINA \
        RGPU=140926_SND104_B_L006_HDX-5 \
        RGSM=SS_GENOMIC \
        SORT_ORDER=coordinate

done

for file in *.add_readgroup.bam; do

    if [[ $file =~ U[1-2] ]]; then
        stub=($(echo "$file" | cut -d "." -f 1,2))
    else
        stub=($(echo "$file" | cut -d "." -f 1))
    fi

    java -jar ${MY_PICARD_BIN} MarkDuplicates \
        INPUT=$file \
        OUTPUT=${stub}.mark_duplicates.bam \
        METRICS_FILE=${stub}.metrics.txt \
        CREATE_INDEX=true

done
