#!/usr/bin/env perl

use strict;
use warnings;

## This script will apply a consistent naming scheme to genes in a gff file. 
## Each gene will be given a user-specified prefix and numbered in multiples of 10 
## (to allow for inclusion of extra tracks if necessary).
## This script assumes that chromosomes or scaffolds are named as "prefix_n" 
## where n is the number of the chromosome or scaffold (the first field of the gff). 
## The script will not work if this field is named otherwise.
## For example, this script was originally written to modify a gff file that contained chromosome info 
## in the first column that read like this "Chr_1, Chr_2, Chr_3.. Chr_n".
## This script also assumes that the gff file is formatted such that each gene track 
## comes before its own exon and cds tracks, so that its exons and cdss 
## can be counted starting at the gene track (setting the count to 0 at this point).

my $usage = "$0 <prefix> <gff file>\n";

if (@ARGV < 2){die ($usage)};

my $prefix = $ARGV[0];
my $input = $ARGV[1];
my $genecount = 0;
my $chrno = 0;
my $exoncount = 0;
my $cdscount = 0;

open my $fh, "<", $input or die ($usage);

while (my $line = <$fh>){

    if ($line =~ /^#/){
        print $line;
        next;}

    chomp $line;
    my @fields = split "\t", $line;
    my @chrattrs = split "_", $fields[0];
    $chrno = $chrattrs[1];

    if ($fields[2] eq "gene"){

        $genecount++;
        $exoncount = 0;
        $cdscount = 0;}

    if ($fields[2] eq "gene"){

        if ($chrno < 10){         

            if ($genecount < 10){
                $line =~ s/ID=.+$/ID=${prefix}0${chrno}g0000${genecount}0;Name=${prefix}0${chrno}g0000${genecount}0/g;}
            elsif (($genecount < 100) && ($genecount >=10)){
                $line =~ s/ID=.+$/ID=${prefix}0${chrno}g000${genecount}0;Name=${prefix}0${chrno}g000${genecount}0/g;}
            elsif (($genecount < 1000) && ($genecount >=100)){
                $line =~ s/ID=.+$/ID=${prefix}0${chrno}g00${genecount}0;Name=${prefix}0${chrno}g00${genecount}0/g;}
            elsif (($genecount < 10000) && ($genecount >=1000)){
                $line =~ s/ID=.+$/ID=${prefix}0${chrno}g0${genecount}0;Name=${prefix}0${chrno}g0${genecount}0/g;}
            elsif ($genecount >= 10000){
                $line =~ s/ID=.+$/ID=${prefix}0${chrno}g${genecount}0;Name=${prefix}0${chrno}g${genecount}0/g;}}

        else {

            if ($genecount < 10){
                $line =~ s/ID=.+$/ID=${prefix}${chrno}g0000${genecount}0;Name=${prefix}${chrno}g0000${genecount}0/g;}
            elsif (($genecount < 100) && ($genecount >=10)){
                $line =~ s/ID=.+$/ID=${prefix}${chrno}g000${genecount}0;Name=${prefix}${chrno}g000${genecount}0/g;}
            elsif (($genecount < 1000) && ($genecount >=100)){
                $line =~ s/ID=.+$/ID=${prefix}${chrno}g00${genecount}0;Name=${prefix}${chrno}g00${genecount}0/g;}
            elsif (($genecount < 10000) && ($genecount >=1000)){
                $line =~ s/ID=.+$/ID=${prefix}${chrno}g0${genecount}0;Name=${prefix}${chrno}g0${genecount}0/g;}
            elsif ($genecount >=10000){
                $line =~ s/ID=.+$/ID=${prefix}${chrno}g${genecount}0;Name=${prefix}${chrno}g${genecount}0/g;}}}

    if ($fields[2] eq "mRNA"){

        if ($chrno < 10){         

            if ($genecount < 10){
                $line =~ s/ID=.+$/ID=mRNA.${prefix}0${chrno}g0000${genecount}0;Name=mRNA.${prefix}0${chrno}g0000${genecount}0;Parent=${prefix}0${chrno}g0000${genecount}0/g;}
            elsif (($genecount < 100) && ($genecount >=10)){
                $line =~ s/ID=.+$/ID=mRNA.${prefix}0${chrno}g000${genecount}0;Name=mRNA.${prefix}0${chrno}g000${genecount}0;Parent=${prefix}0${chrno}g000${genecount}0/g;}
            elsif (($genecount < 1000) && ($genecount >=100)){
                $line =~ s/ID=.+$/ID=mRNA.${prefix}0${chrno}g00${genecount}0;Name=mRNA.${prefix}0${chrno}g00${genecount}0;Parent=${prefix}0${chrno}g00${genecount}0/g;}
            elsif (($genecount < 10000) && ($genecount >=1000)){
                $line =~ s/ID=.+$/ID=mRNA.${prefix}0${chrno}g0${genecount}0;Name=mRNA.${prefix}0${chrno}g0${genecount}0;Parent=${prefix}0${chrno}g0${genecount}0/g;}
            elsif ($genecount >= 10000){
                $line =~ s/ID=.+$/ID=mRNA.${prefix}0${chrno}g${genecount}0;Name=mRNA.${prefix}0${chrno}g${genecount}0;Parent=${prefix}0${chrno}g${genecount}0/g;}}

        else {

            if ($genecount < 10){
                $line =~ s/ID=.+$/ID=mRNA.${prefix}${chrno}g0000${genecount}0;Name=mRNA.${prefix}${chrno}g0000${genecount}0;Parent=${prefix}${chrno}g0000${genecount}0/g;}
            elsif (($genecount < 100) && ($genecount >=10)){
                $line =~ s/ID=.+$/ID=mRNA.${prefix}${chrno}g000${genecount}0;Name=mRNA.${prefix}${chrno}g000${genecount}0;Parent=${prefix}${chrno}g000${genecount}0/g;}
            elsif (($genecount < 1000) && ($genecount >=100)){
                $line =~ s/ID=.+$/ID=mRNA.${prefix}${chrno}g00${genecount}0;Name=mRNA.${prefix}${chrno}g00${genecount}0;Parent=${prefix}${chrno}g00${genecount}0/g;}
            elsif (($genecount < 10000) && ($genecount >=1000)){
                $line =~ s/ID=.+$/ID=mRNA.${prefix}${chrno}g0${genecount}0;Name=mRNA.${prefix}${chrno}g0${genecount}0;Parent=${prefix}${chrno}g0${genecount}0/g;}
            elsif ($genecount >=10000){
                $line =~ s/ID=.+$/ID=mRNA.${prefix}${chrno}g${genecount}0;Name=mRNA.${prefix}${chrno}g${genecount}0;Parent=${prefix}${chrno}g${genecount}0/g;}}}

    if ($fields[2] eq "exon"){

        $exoncount++;

        if ($chrno < 10){         

            if ($genecount < 10){
                $line =~ s/ID=.+$/ID=exon.$exoncount.${prefix}0${chrno}g0000${genecount}0;Name=exon.$exoncount.${prefix}0${chrno}g0000${genecount}0;Parent=mRNA.${prefix}0${chrno}g0000${genecount}0/g;}
            elsif (($genecount < 100) && ($genecount >=10)){
                $line =~ s/ID=.+$/ID=exon.$exoncount.${prefix}0${chrno}g000${genecount}0;Name=exon.$exoncount.${prefix}0${chrno}g000${genecount}0;Parent=mRNA.${prefix}0${chrno}g000${genecount}0/g;}
            elsif (($genecount < 1000) && ($genecount >=100)){
                $line =~ s/ID=.+$/ID=exon.$exoncount.${prefix}0${chrno}g00${genecount}0;Name=exon.$exoncount.${prefix}0${chrno}g00${genecount}0;Parent=mRNA.${prefix}0${chrno}g00${genecount}0/g;}
            elsif (($genecount < 10000) && ($genecount >=1000)){
                $line =~ s/ID=.+$/ID=exon.$exoncount.${prefix}0${chrno}g0${genecount}0;Name=exon.$exoncount.${prefix}0${chrno}g0${genecount}0;Parent=mRNA.${prefix}0${chrno}g0${genecount}0/g;}
            elsif ($genecount >= 10000){
                $line =~ s/ID=.+$/ID=exon.$exoncount.${prefix}0${chrno}g${genecount}0;Name=exon.$exoncount.${prefix}0${chrno}g${genecount}0;Parent=mRNA.${prefix}0${chrno}g${genecount}0/g;}}

        else {

            if ($genecount < 10){
                $line =~ s/ID=.+$/ID=exon.$exoncount.${prefix}${chrno}g0000${genecount}0;Name=exon.$exoncount.${prefix}${chrno}g0000${genecount}0;Parent=mRNA.${prefix}${chrno}g0000${genecount}0/g;}
            elsif (($genecount < 100) && ($genecount >=10)){
                $line =~ s/ID=.+$/ID=exon.$exoncount.${prefix}${chrno}g000${genecount}0;Name=exon.$exoncount.${prefix}${chrno}g000${genecount}0;Parent=mRNA.${prefix}${chrno}g000${genecount}0/g;}
            elsif (($genecount < 1000) && ($genecount >=100)){
                $line =~ s/ID=.+$/ID=exon.$exoncount.${prefix}${chrno}g00${genecount}0;Name=exon.$exoncount.${prefix}${chrno}g00${genecount}0;Parent=mRNA.${prefix}${chrno}g00${genecount}0/g;}
            elsif (($genecount < 10000) && ($genecount >=1000)){
                $line =~ s/ID=.+$/ID=exon.$exoncount.${prefix}${chrno}g0${genecount}0;Name=exon.$exoncount.${prefix}${chrno}g0${genecount}0;Parent=mRNA.${prefix}${chrno}g0${genecount}0/g;}
            elsif ($genecount >=10000){
                $line =~ s/ID=.+$/ID=exon.$exoncount.${prefix}${chrno}g${genecount}0;Name=exon.$exoncount.${prefix}${chrno}g${genecount}0;Parent=mRNA.${prefix}${chrno}g${genecount}0/g;}}}

    if ($fields[2] eq "CDS"){

        $cdscount++;

        if ($chrno < 10){         

            if ($genecount < 10){
                $line =~ s/ID=.+$/ID=cds.$cdscount.${prefix}0${chrno}g0000${genecount}0;Name=cds.$cdscount.${prefix}0${chrno}g0000${genecount}0;Parent=mRNA.${prefix}0${chrno}g0000${genecount}0/g;}
            elsif (($genecount < 100) && ($genecount >=10)){
                $line =~ s/ID=.+$/ID=cds.$cdscount.${prefix}0${chrno}g000${genecount}0;Name=cds.$cdscount.${prefix}0${chrno}g000${genecount}0;Parent=mRNA.${prefix}0${chrno}g000${genecount}0/g;}
            elsif (($genecount < 1000) && ($genecount >=100)){
                $line =~ s/ID=.+$/ID=cds.$cdscount.${prefix}0${chrno}g00${genecount}0;Name=cds.$cdscount.${prefix}0${chrno}g00${genecount}0;Parent=mRNA.${prefix}0${chrno}g00${genecount}0/g;}
            elsif (($genecount < 10000) && ($genecount >=1000)){
                $line =~ s/ID=.+$/ID=cds.$cdscount.${prefix}0${chrno}g0${genecount}0;Name=cds.$cdscount.${prefix}0${chrno}g0${genecount}0;Parent=mRNA.${prefix}0${chrno}g0${genecount}0/g;}
            elsif ($genecount >= 10000){
                $line =~ s/ID=.+$/ID=cds.$cdscount.${prefix}0${chrno}g${genecount}0;Name=cds.$cdscount.${prefix}0${chrno}g${genecount}0;Parent=mRNA.${prefix}0${chrno}g${genecount}0/g;}}

        else {

            if ($genecount < 10){
                $line =~ s/ID=.+$/ID=cds.$cdscount.${prefix}${chrno}g0000${genecount}0;Name=cds.$cdscount.${prefix}${chrno}g0000${genecount}0;Parent=mRNA.${prefix}${chrno}g0000${genecount}0/g;}
            elsif (($genecount < 100) && ($genecount >=10)){
                $line =~ s/ID=.+$/ID=cds.$cdscount.${prefix}${chrno}g000${genecount}0;Name=cds.$cdscount.${prefix}${chrno}g000${genecount}0;Parent=mRNA.${prefix}${chrno}g000${genecount}0/g;}
            elsif (($genecount < 1000) && ($genecount >=100)){
                $line =~ s/ID=.+$/ID=cds.$cdscount.${prefix}${chrno}g00${genecount}0;Name=cds.$cdscount.${prefix}${chrno}g00${genecount}0;Parent=mRNA.${prefix}${chrno}g00${genecount}0/g;}
            elsif (($genecount < 10000) && ($genecount >=1000)){
                $line =~ s/ID=.+$/ID=cds.$cdscount.${prefix}${chrno}g0${genecount}0;Name=cds.$cdscount.${prefix}${chrno}g0${genecount}0;Parent=mRNA.${prefix}${chrno}g0${genecount}0/g;}
            elsif ($genecount >=10000){
                $line =~ s/ID=.+$/ID=cds.$cdscount.${prefix}${chrno}g${genecount}0;Name=cds.$cdscount.${prefix}${chrno}g${genecount}0;Parent=mRNA.${prefix}${chrno}g${genecount}0/g;}}}

    print "$line\n";}
close $fh;
