#!/usr/bin/env bash

# This script is used to create fasta dictionaries and indexes,
# create bam indexes, and call variants in vcf format from
# the bams and the fastas.

## Set dircetories to use at runtime.
MY_GATK_BINARY=/home/mark/bioinformatics/src/GenomeAnalysisTK.jar
MY_GATK_WORKING_DIRECTORY=/media/rdrive/CCDM_Program_6-DENTOM-SE01498/samtools/
MY_GENOMES_WORKING_DIRECTORY=/home/mark/Research/2016/02.16/data/genomes/
MY_PICARD_BIN=/home/mark/bioinformatics/src/picard/dist/picard.jar

echo "Script started at $(date)..."

## Change to first working directory to process reference files.
cd ${MY_GENOMES_WORKING_DIRECTORY}

echo "Creating fasta dictionaries..."

## Create sequence dictionaries for all fasta files being processed.
for file in *.fasta; do

    stub=($(echo "$file" | cut -d "." -f 1));

    java -jar ${MY_PICARD_BIN} \
        CreateSequenceDictionary \
        R=$file \
        O=${stub}.dict

done

echo "Done!"

echo "Indexing fasta files..."

## Create fasta indexes for all fasta files being processed.
for file in *.fasta; do

    samtools faidx $file

done

echo "Done!"

## Change to second working directory to process bam files.
cd ${MY_GATK_WORKING_DIRECTORY}

echo "Building bam indexes with BuildBamIndex and detecting variants with HaplotypeCaller..."

## Build bam indexes with Picard Tools and call haplotypes with GATK for all bam files.
for file in *mark_duplicates.bam; do

    echo "Processing file ${file}..."

    java -jar ${MY_PICARD_BIN} \
        BuildBamIndex \
        I=$file \
        O=${file}.bai

    test_stub=($(echo "$file" | cut -d "." -f 1))

    if [[ $file =~ U[12] ]]; then
        bam_stub=($(echo "$file" | cut -d "." -f 1,2))
    else
        bam_stub=($(echo "$file" | cut -d "." -f 1))
    fi

    for reference in ${MY_GENOMES_WORKING_DIRECTORY}*.fasta; do 

    genome_stub=($(echo "$reference" | cut -d "." -f 2 | rev | cut -d "/" -f 1 | rev))

        if [ "$test_stub" = "$genome_stub" ]; then

            java -jar ${MY_GATK_BINARY} \
                -R $reference \
                -T HaplotypeCaller \
                -I $file \
                --out ${bam_stub}.vcf \
                --maxNumHaplotypesInPopulation 1
        fi

    done

    echo "Finished processing file ${file}!"

done

## Let user know script has finished.
echo "Script finished at $(date)"
