#!/usr/bin/env bash

set -e

# This 6 step pipeline will extract regions where small RNA reads map exactly and
# will determine and validate hairpins within these regions using RNAFold and HMMIR,
# respectively. It will then map the reads back to the valid hairpins and print
# a coverage bedgraph for user inspection.

genome=$1
reads=$2
annotations=$3

genome_prefix=$(echo "$genome" | rev | cut -d "." -f2 | rev)
reads_prefix=$(echo "$reads" | rev | cut -d "." -f2 | rev)

mkdir $reads_prefix

# Step 1
#
# Collapse fastq file containing small RNA reads and build bowtie index for
# the genome they are to be mapped to. Then map the reads to the genome, only
# outputting exact matches.

fastx_collapser -i $reads -o ${reads_prefix}.fasta
bowtie-build $genome $genome_prefix

bowtie \
-f \
-v 0 \
-a \
--sam \
$genome_prefix \
${reads_prefix}.fasta \
${reads_prefix}.sam

# Step 2
#
# Extract regions flanking sites where small RNA reads map to a depth of 50 x or more.
# Scan these sequences for hairpins using RNAfold and extract only the hairpins, then
# validate them using HMMIR. All procedures wrapped in extractFlanks.pl.

perl extractFlanks.pl \
--fasta $genome \
--sam ${reads_prefix}.sam \
--direction both \
--flank 60 \
--stack 50 \
--output ${reads_prefix}.csv

# Step 3
#
# Filter out those that are predicted to be true hairpins and change Us to Ts.
# Then collapse this fasta file and map the validated hairpins back to the genome.

awk '/TRUE/{split ($0,a,",");print a[1]"\n"a[2]}' ${reads_prefix}.csv | \
perl -pe "s/U/T/g" > ${reads_prefix}.valid.fasta
fastx_collapser -i ${reads_prefix}.valid.fasta -o ${reads_prefix}.uniqueHPs.fasta

bowtie \
-f \
-v 0 \
-a \
--sam \
$genome_prefix \
${reads_prefix}.uniqueHPs.fasta \
${reads_prefix}.uniqueHPs_mapped.sam

# Step 4
#
# Create a bed file from the sam file containing the unique hairpins mapped back
# to the genome.
# Merge the bed file containing the sequences mapped back to the genome to create a
# single hairpin locus for each read stack.
# Then extract these intergenic hairpin sequences.

awk '!/@/{seqlen = length($10);  seqend = $4 + seqlen; print $3"\t"$4"\t"seqend"\n"}' ${reads_prefix}.uniqueHPs_mapped.sam | sort -V -k1,1 -k2,2 > ${reads_prefix}.bed
bedtools merge -d -1 -i ${reads_prefix}.bed > ${reads_prefix}.merged.bed
bedtools intersect -v -a ${reads_prefix}.merged.bed -b $annotations > ${reads_prefix}.intergenic.bed
bedtools getfasta -fi $genome -bed ${reads_prefix}.intergenic.bed -fo ${reads_prefix}.intergenic.fasta

# Step 5
#
# Map the reads back to the hairpin sequences to determine the depth of coverage across the
# hairpins.

bowtie-build ${reads_prefix}.intergenic.fasta ${reads_prefix}.intergenic

bowtie \
-q \
-v 0 \
-a \
--sam \
${reads_prefix}.intergenic \
$reads \
${reads_prefix}.mapped_to_hairpins.sam

# Step 6
#
# Determine coverage across the hairpins and create a bedgraph.

samtools view -bS ${reads_prefix}.mapped_to_hairpins.sam > ${reads_prefix}.mapped_to_hairpins.bam
samtools sort ${reads_prefix}.mapped_to_hairpins.bam ${reads_prefix}.mapped_to_hairpins.sort
bedtools genomecov -bg -ibam ${reads_prefix}.mapped_to_hairpins.sort.bam -g ${reads_prefix}.intergenic.fasta > ${reads_prefix}.bedgraph

# Clean up
mv ${reads_prefix}.fasta ${reads_prefix}.sam ${reads_prefix}.csv ${reads_prefix}.uniqueHPs.fasta ${reads_prefix}.uniqueHPs_mapped.sam ${reads_prefix}.intergenic.bed ${reads_prefix}.intergenic.fasta ${reads_prefix}.merged.bed ${reads_prefix}.bed ${reads_prefix}.valid.fasta ${reads_prefix}.mapped_to_hairpins.sam ${reads_prefix}.mapped_to_hairpins.bam ${reads_prefix}.mapped_to_hairpins.sort.bam ${reads_prefix}.bedgraph $reads_prefix
rm *ebwt *fai *extracts folded.fa hairpins.fa ip.txt rna.ps removed.txt 
