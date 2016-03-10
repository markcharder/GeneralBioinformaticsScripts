#!/usr/bin/env bash

## This script is used to get TAXIDs from NCBI (through the curl command) based on sequence accessions.
ACCESSIONS=$@

for f in ${ACCESSIONS[@]}; do
       echo -n -e "$f\t"
       curl -s "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$f&rettype=fasta&retmode=xml" |\
       grep TSeq_taxid |\
       cut -d '>' -f 2 |\
       cut -d '<' -f 1 |\
       tr -d "\n"
       echo
    done
