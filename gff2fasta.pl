#!/usr/bin/env perl

use strict;
use warnings;
use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::Fasta;

## Use this script to obtain three fasta sequences from a gff file and it's corresponding fasta.
## This script uses bioperl.

my $cdsoutput = Bio::SeqIO->new(-format => 'fasta', -file => ">$ARGV[2].cds.fasta");
my $pepoutput = Bio::SeqIO->new(-format => 'fasta', -file=> ">$ARGV[2].pep.fasta");
my $geneoutput = Bio::SeqIO->new(-format => 'fasta', -file => ">$ARGV[2].gene.fasta");

my $infasta = $ARGV[0];
my $db = Bio::DB::Fasta->new($infasta);
print "Fasta parsed\n";

my %cds;
my $mrnaname;
my $frame;

open my $gff, "<", $ARGV[1] or die "Specify gff";
while (my $line = <$gff>){

    chomp $line;

    if ($line =~ /^#/){

        if (eof($gff)){

        my $mergedcds;

        foreach my $key (sort {$a <=> $b} keys %cds){
            my $coord = $cds{$key};
            my @cdscoord = split (" ", $coord);
            my $cdsseq = $db->seq($cdscoord[0],$cdscoord[1],$cdscoord[2]);
            $mergedcds .= $cdsseq;}

        my $cdsout = Bio::Seq->new(
            -seq => $mergedcds,
            -id => $mrnaname,
            -display_id => $mrnaname,
            -alphabet => "dna",);

        if ($frame eq "-") {
           $cdsout = $cdsout->revcom();}

        my $pepout = $cdsout->translate();
        $cdsoutput->write_seq($cdsout);
        $pepoutput->write_seq($pepout);
        next;}

    else{
        next;}}

    my @array = split "\t", $line;
    my $type = $array[2];

    if ($type eq "gene"){

        my @attrs = split ";", $array[8];
        $attrs[0] =~ s/ID=//g;
        my $genename = $attrs[0];
        my $genestart = $array[3];
        my $geneend = $array[4];
        my $geneseq = $db->seq($array[0], $genestart, $geneend);
        my $geneout = Bio::Seq->new(
            -seq => $geneseq,
            -id => $genename,
            -displayid => $genename,
            -alphabet => "dna",);

        if ($array[6] eq "-"){
            $geneout = $geneout->revcom();}
        $geneoutput->write_seq($geneout);}

    if (($type eq "mRNA") and ($. > 4)){

        my $mergedcds;

        foreach my $key (sort {$a <=> $b} keys %cds){
            my $coord = $cds{$key};
            my @cdscoord = split (" ", $coord);
            my $cdsseq = $db->seq($cdscoord[0],$cdscoord[1],$cdscoord[2]);
            $mergedcds .= $cdsseq;}

        my $cdsout = Bio::Seq->new(
            -seq => $mergedcds,
            -id => $mrnaname,
            -display_id => $mrnaname,
            -alphabet => "dna",);

        if ($frame eq "-") {
           $cdsout = $cdsout->revcom();}

        my $pepout = $cdsout->translate();
        $cdsoutput->write_seq($cdsout);
        $pepoutput->write_seq($pepout);

        my @attrs = split ";", $array[8];
        $attrs[0] =~ s/ID=//g;
        $mrnaname = $attrs[0];
        $frame = $array[6];
        %cds = ();}

    elsif ($type eq "mRNA"){
        my @attrs = split (";", $array[8]);
        $attrs[0] =~ s/ID=//g;
        $mrnaname = $attrs[0];
        $frame = $array[6];}

    elsif ($type eq "CDS"){
        my $cdscoord = $array[0] . " " . $array[3] . " " . $array[4];
        $cds{$array[3]} = $cdscoord;}}

close $gff;


