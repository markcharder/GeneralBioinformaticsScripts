#!/usr/bin/env bash

# Pipeline used to define 'effector-like' sequences for each genome.
# Basically, sequences with signal peptides but no transmembrane helices.

softwares=/home/mark/bioinformatics/src
signalp=/home/mark/bioinformatics/src/signalp-4.1/signalp
tmhmm=/home/mark/bioinformatics/src/tmhmm-2.0c/bin/tmhmm
interpro=/media/rdrive/CCDM_Program_6-DENTOM-SE01498/mark_bioinf_work/interproscan-5.17-56.0/interproscan.sh
pfam=/home/mark/bioinformatics/src/PfamScan/pfam_scan.pl
export PERL5LIB=/home/mark/bioinformatics/src/PfamScan/:$PERL5LIB

ln -s ${softwares}/repbase20.05_aaSeq_cleaned_TE.fa .

ARRAY=`cat $1`

for file in ${ARRAY[@]}; do

    out=`echo $file | cut -d"." -f1`

    # Get all protein sequences from gff3.
    perl gff2fasta.pl ${out}.fa ${out}.gff3 ${out}

    #Remove stop codons from protein sequences for signalP.
    perl -pe "s/\*//g" ${out}.pep.fasta > temp
    mv temp ${out}.pep.fasta

    rm *.cds.fasta *.gene.fasta

done

# Make blast database of repbase.
makeblastdb -in repbase20.05_aaSeq_cleaned_TE.fa -out repbase20.05_aaSeq_cleaned_TE -dbtype prot

for file in *.pep.fasta; do

    out=`echo $file | cut -d"." -f1`

    # Blast all proteins against repbase.
    blastp \
    -query ${file} \
    -db repbase20.05_aaSeq_cleaned_TE \
    -max_hsps 1 \
    -max_target_seqs 1 \
    -outfmt "6 qseqid pident qcovs evalue" \
    -out ${out}.blast

    # Filter blast output to remove significant hits to repbase.
    awk '{if ($2 > 30 && $3 > 30 && $4 > e-10){next} else {print $1}}' \
    ${out}.blast > ${out}.nontransposases
    perl /home/mark/Research/2016/04.16/scripts/get_fasta_seqs_from_list.pl \
    ${out}.nontransposases ${file} > ${out}.filtered
    rm *.blast *.nontransposases

    # Split fastas to run SingalP and print out a status on how many sequences in each file.
    gt splitfasta -numfiles 3 ${out}.filtered

    test1=`grep ">" ${out}.filtered.1 | wc`
    test2=`grep ">" ${out}.filtered.2 | wc`
    test3=`grep ">" ${out}.filtered.3 | wc`

    echo "$test1 for ${out}.filtered.1"
    echo "$test2 for ${out}.filtered.1"
    echo "$test3 for ${out}.filtered.1"

    for temp_file in ${out}.filtered.*; do

        $signalp \
        $temp_file >> ${out}.signalP

    done

    rm ${out}.filtered.*

done

for file in *.signalP; do

    out=`echo $file | cut -d"." -f1`

    # Get a list of proteins with signal peptides.
    awk '{if ($10 == "Y") {print $1}}' \
    $file > ${out}.signalPList

    # Use signal peptide-containing protein lest to get entries from fasta.
    perl /home/mark/Research/2016/04.16/scripts/get_fasta_seqs_from_list.pl \
    ${out}.signalPList ${out}.pep.fasta > ${out}.signalp.fasta

    rm *signalP

done

for file in *.signalp.fasta; do

    out=`echo $file | cut -d"." -f1`

    # Use tmhmm to detect transmembrane domains in proteins with signal peptides.
    $tmhmm \
    $file > ${out}.tmhmm

    rm *signalp.fasta

done

for file in *.tmhmm; do

    out=`echo $file | cut -d"." -f1`

    # Get a list of proteins with no transmembrane helices.
    awk '{if ($0 ~ /Number of predicted TMHs:  0/){print $2}}' \
    $file > ${out}.tmhmmList

    # Use this list to collect fasta entries from the signal peptide fasta.
    perl /home/mark/Research/2016/04.16/scripts/get_fasta_seqs_from_list.pl \
    ${out}.tmhmmList ${out}.pep.fasta > ${out}.signalp.tmhmm.fa

    rm *tmhmmList

done

for file in *.signalp.tmhmm.fa; do

    out=`echo $file | cut -d"." -f1`

    # Use pfam to look for protein domains, including GPI anchors.
    perl $pfam \
    --outfile ${out}.pfam \
    -cpu 6 \
    -fasta $file \
    -dir /home/mark/bioinformatics/src/

    # Get a list of proteins with GPI anchored proteins.
    awk '{if ($6 == "GPI-anchored"){print $1}}' \
    ${out}.pfam > ${out}.gpi

    # Use the list of GPI anchored proteins to get a list of proteins without GPI anchors.
    grep -v \
    ${out}.gpi ${out}.pfam | \
    awk '{if ($0 !~ /#/){print $1}}' > ${out}.signalp.tmhmm.gpi
    
    # Use the list of proteins without GPI anchors to get entries from the fasta.
    perl /home/mark/Research/2016/04.16/scripts/get_fasta_seqs_from_list.pl \
    ${out}.signalp.tmhmm.gpi $file > ${out}.effectors.fa

    # Get a list of effectors from the effectors fasta file generated in the previous step.
    grep ">" ${out}.effectors.fa | perl -pe "s/>//g" > ${out}.effectors.list

done

# Clean up remaining guff.
rm -r *.gpi *.tmhmm.fa *.tmhmm *.signalPList *pep* *index *filtered repbase20.05_aaSeq_cleaned_TE.* TMHMM*
