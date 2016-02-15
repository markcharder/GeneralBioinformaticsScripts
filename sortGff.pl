#!/usr/bin/env perl

use strict;
use warnings;

## This script will sort a gff3 file, taking into account the lines beginning with "#". 
## Each scaffold must be named as "prefix_n" where "prefix" is the name of the scaffold 
## and "n" is the scaffold number. In this case as scaffolds were complete chromosomes,
## they were numbered from largest to smallest.

my $firsthead = 0;
my %chrids;
my $chrnumber = "";

open my $fh, "<", $ARGV[0] or die;

while (my $line = <$fh>){

    chomp $line;

    if ($line =~ /^##g/){
        $firsthead = $line;}

    if ($line =~ /^##s/){

        my @chrattrs = split / |_/, $line;
        $chrnumber = $chrattrs[2];

        push @{$chrids{$chrnumber}}, $firsthead, $line;}

    if ($line =~ /^Chr/){
        push @{$chrids{$chrnumber}}, $line;}

    if ($line =~ /^###/){
        push @{$chrids{$chrnumber}}, $line;}}

close $fh;

for my $i (sort {$a<=>$b} keys %chrids){

    my @array = @{$chrids{$i}};

    my @headerarray = @array[0..1];

    my @featsarray = @array[2..$#array - 1];

    my $end = $array[-1];

    my @newfeatsarray;

    foreach my $j (@featsarray){
        my @fieldsarray = split "\t", $j;
        push @newfeatsarray, \@fieldsarray;}

    my @sortedarray = sort {$a->[3] <=> $b->[3]} @newfeatsarray;

    print join ("\n", @headerarray) . "\n";

    foreach my $j (@sortedarray){
        print join ("\t", @$j) . "\n";}

    print "$end\n";}
