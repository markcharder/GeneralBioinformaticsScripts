#!/usr/bin/env perl

#Usage: reciprocalBlastRetriever.pl <blast> <reciprocal blast>
#This script takes a blast table output from a blastp analysis and a blast table output from another, reciprocal blast analysis and places corresponding reciprocal blasts side by side.

use strict;
use warnings;

my %reciprocalhits;

open my $fh1, "<", $ARGV[0] or die "Please specify first input.";
while (my $line = <$fh1>){
    chomp $line;
    $reciprocalhits{$line} = "";}
close $fh1;

open my $fh2, "<", $ARGV[1] or die "Please specify second input.";
while (my $line = <$fh2>){
    my @fields = split "\t", $line;
    for my $i (keys %reciprocalhits){
        if ($i =~ m/\Q$fields[2]\E.+\Q$fields[0]/){
            print "$i\t$line";}}}
close $fh2;
