#!/usr/bin/env perl

## Script used to reformat default output from Pawsey script for blasting with a database split.
## Although Pawsey has provided some extremely useful scripts for running blast in parallel,
## when splitting a blast database, it seems only possible to produce one output (not sure which number in BLAST specifications...?).
## In order to create a tabular blast output, same as "-outfmt '6 qseqid qseqlen sseqid sseqlen length evalue pident'", this script
## was used.
use strict;
use warnings;

## Define global variables
my $query_id;
my $subject_id;
my $evalue;
my $percent_identity;
my $length_count=0;
my $query_length;
my $subject_length;
my $alignment_coverage;

warn "Reformatting blast...\n";

## File is specified as first argument to script.
open my $fh, "<", $ARGV[0] or die "Please specify input file...";
while (my $line = <$fh>){

    chomp $line;

    ## Get the query ID and reset $length_count.
    if ($line =~ /Query=/){
        $length_count=0;
        my @fields = split /\s+/, $line;
        $query_id = $fields[1];}

    ## If line contains length information, increment the length_count.
    ## If length count is one, i.e. this is the query length, define $query_length.
    ## Or else, define $subject_length.
    if ($line =~ /^Length=/){
        $length_count++;
        my @fields = split "=", $line;
        if ($length_count == 1){
            $query_length = $fields[1]}
        else{
            $subject_length = $fields[1];}}

    ## Get the subject ID.
    if ($line =~ />/){
        $line =~ s/>//g;
        $subject_id = $line;}

    ## Get the evalue of the last subject ID.
    if ($line =~ /\s+Score\s+=/){
        my @fields = split /\s+/, $line;
        $fields[8] =~ s/\,//g;
        $evalue = $fields[8];}   

    ## Get the percent_identity of the alignment.
    if ($line =~ /\s+Identities\s+=/){
        my @fields = split /\s+/, $line;
        $fields[4] =~ s/\(|\)|\,//g;
        $percent_identity = $fields[4];}

    ## Get the alignment coverage by subtracting the query start from the query end and adding 1.
    ## Print all information gained so far.
    if ($line =~ /Query\s+/){
        my @fields = split /\s+/, $line;
        $alignment_coverage = $fields[3] - $fields[1] + 1;
        print "$query_id\t$query_length\t$subject_id\t$subject_length\t$alignment_coverage\t$evalue\t$percent_identity\n";}}

close $fh;

warn "Script finished!\n";
