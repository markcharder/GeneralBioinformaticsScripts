#!/usr/bin/env perl

use strict;
use warnings;

## This script will take a fasta file and determine the average percentage of cysteins per sequence.
## It will then determine whether each individual sequence has a percentage of cysteins higher than the average by a user-defined ratio.
## It will produce two output files, the first containing all sequences, average cystein content and ratio to global average.
## The second contains a subset with a ratio of cystein content to the global average above the user-defined threshold.
## Probably most useful for whole genomic proteomes or at least a random subset.

my %headers;
my $header = "";
my $count;
my $sequence = "";

my $input = $ARGV[0];
my $output = $ARGV[1];
my $output1 = $ARGV[2];
my $threshold = $ARGV[3];
my $usage= "$0 <fasta> <output file 1> <output file 2> <threshold>\n" .
"output file 1 = all sequences\n" .
"output file 2 = a subset with cystein content above a certain threshold\n" .
"threshold = ratio of percent cystein content to average cystein content for input fasta sequences\n";

if (@ARGV != 4){die "\n$usage\n"};

open my $infile, "<", $input or die "\n$usage\n";

while (my $line = <$infile>){

    chomp $line;

    if (($line =~ />/)||(eof($infile))){
        if ($. > 1){
            $count = ($sequence =~ tr/C//);
            my $seqlen = length($sequence);
            push @{$headers{$header}}, $header, $count, $seqlen;
            $sequence = "";}
        
        $header = $line;
        $header =~ s/>//g;
        $header =~ s/mRNA\.//g;}

    if ($line !~ />/){
        $sequence .= $line;}}

close $infile;

my $cyscount = 0;

for my $i (keys %headers){
    my @array = @{$headers{$i}};
    my $cyscontent = $array[1] / $array[2];
    $cyscount = $cyscount + $cyscontent;}

my $numberofsequences = keys %headers;
my $averagecyscontent = $cyscount / $numberofsequences * 100;

open my $outfile, ">", $output or die "\n$usage\n";
print $outfile "Number of sequences analysed: $numberofsequences\n\n" .
      "Average cystein content (%): $averagecyscontent\n\n" .
      "Sequence\tCysteins\t% content\tRatio\n";

my %highcyscontent;

for my $i (sort keys %headers){
    my @array = @{$headers{$i}};
    my $cyscontent = $array[1] / $array[2] * 100;
    my $ratio = $cyscontent / $averagecyscontent;
    if ($ratio > $threshold) {
        push @{$highcyscontent{$array[0]}}, $array[0], $array[1], $cyscontent, $ratio;}
    print $outfile "$array[0]\t$array[1]\t$cyscontent\t$ratio\n";}
close $outfile;

open $outfile, ">", $output1 or die "\n$usage\n";
print $outfile "Sequences with ratio of cystein content to average above $threshold\n\n" .
      "Sequence\tCysteins\t% content\tRatio\n";
for my $i (sort keys %highcyscontent){
    my @array = @{$highcyscontent{$i}};
    print $outfile join ("\t", @array) . "\n";}
close $outfile;
