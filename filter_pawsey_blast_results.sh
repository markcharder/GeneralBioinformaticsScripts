#!/usr/bin/env bash

## This script was used to filter the re-formatted blast outputs from Pawsey for hits outside Sclerotinia
## With an e-value below e1-10, a % identity more than 30 and an alignment coverage of both sequences more
## than 30 %. The user can specify thresholds for these values.

if [ "$#" -lt 5 ]; then
    echo ""
    echo "Usage:"
    echo "$0 <blast results> <output_file> <% identity threshold> <evalue threshold> <% coverage threshold>"
    echo ""
    exit
fi

in1=$1
in2=$2
in3=$3
in4=$4
in5=$5

calculate_coverage=$(bc -l <<< "$in5/100")

echo "Script begun $(date)"
echo "Parsing blast output..."

awk -v identity_threshold=$in3 \
    -v evalue_threshold="1e-"$in4 \
    -v coverage_threshold=$calculate_coverage \
    '{

    ## Get the information you need from the various fields of the blast table.

    coverage_position = $(NF - 2);

    subject_length_position = $(NF - 3);

    query_length_position = $2;

    evalue_position = $(NF - 1);

    percent_identity_position = $(NF);

    query_coverage = (coverage_position / query_length_position);

    subject_coverage = (coverage_position / subject_length_position);

    ## Check if the description contains "Sclerotinia sclerotiorum",
    ## the e-value is below the threshold,
    ## the % identity is above the threshold,
    ## and the coverage of both query and subject is above the threshold.

    if ($0 !~ /Sclerotinia sclerotiorum/ &&

        evalue_position < evalue_threshold &&

        percent_identity_position > identity_threshold &&

        query_coverage > coverage_threshold &&

        subject_coverage > coverage_threshold){

            print $0}
}' $in1 | sort -V -k1,1 -k5,5 | sort -V -u -k1,1 >> ${in2}.txt ## Sort the output and print only unique queries to remove multiple HSPs.

## Let the user know the script has finished.
echo "Script finished $(date)"
