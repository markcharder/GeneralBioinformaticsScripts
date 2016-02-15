#!/usr/bin/perl

use strict;
use warnings;

# Use this script to check if a peptide fasta file contains premature stops, missing stops or non-mehionine starts.

my $stop = "";
my $start = "";
my @seqarray = "";
my $header = "";
my $detectstops = 0;
my $linecount = 0;
my $switch = 0;

while (my $line = <>){
    
    chomp $line;
    if (($line =~ />/)&&($. > 1)){
        $start = substr($seqarray[1],0,1);
        $stop = chop $seqarray[-1];
        chop $seqarray[-1];

        foreach my $i (@seqarray){
            if ($i =~ /\*/){
                $switch = 1;
                $detectstops = 1;}}

        unless (@seqarray eq ""){

        if ($detectstops == 1){
            $switch = 1;
            print "$header:\tPremature stop(s)\n";}

        if ($start !~ /M/){
            $switch = 1;
            print "$header:\tNon-methionine start\n";}

        if ($stop !~ /\*/){
            $switch = 1;
            print "$header:\tNo terminal stop\n";}}

        $detectstops = 0;
        $header = $line;
        @seqarray = "";}

    if ($line !~ />/){
        push @seqarray, $line;}

    if ((eof) && ($switch == 0)){
        print "No issues detected\n";}
}
