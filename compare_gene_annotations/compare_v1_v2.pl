#!/usr/bin/env perl

use strict;
use warnings;

my $identicalscount=0;
my %identicals;
my %more90;
my $more90=0;
my $pless10=0;
my %less90;
my $pdecrease=0;
my $pincrease=0;
my %decreased;
my %increased;
my $more80=0;
my $more90less10=0;

open my $fh, "<", $ARGV[0] or die;
while (my $line = <$fh>){
    chomp $line;
    $line =~ s/mRNA\.|mRNA//g;
    my @fields = split "\t", $line;
    my @scleroid = split /\|/, $fields[2];
    if ($fields[3] == $fields[1]){
        if ($fields[6] == 100){
            $identicalscount++;
            push @{$identicals{$fields[0]}}, $fields[0], $scleroid[3];}}
        if ($fields[6] != 100){
            if ($fields[6] > 90){
                $more90++;
                push @{$more90{$fields[0]}},$fields[0],$scleroid[3],$fields[6];}
            if (($fields[6] < 90)&&($fields[6] > 80)){
                $more80++;}}
    if ($fields[3] != $fields[1]){
        if ($fields[3] > $fields[1]){
            $pdecrease = 1 - ($fields[1]/$fields[3]);
            push @{$decreased{$fields[0]}},$scleroid[3],$fields[6],$pdecrease;
            if ($fields[6] > 90){
                $more90less10++;}
            if ($pdecrease < 0.1){
                $pless10++;}}
        if ($fields[3] < $fields[1]){
            $pincrease = 1 - ($fields[3]/$fields[1]);
            push @{$increased{$fields[0]}},$scleroid[3],$fields[6],$pincrease;
            if ($fields[6] > 90){
                $more90less10++;}
            if ($pincrease < 0.1){
                $pless10++;}}}}
close $fh;

open my $fh2, '>', $ARGV[1] or die;
print $fh2 "Number of identical sequences: $identicalscount\t" . "(" . ($identicalscount/11132*100) . " %)\n";
print $fh2 "Number of sequences 90-100 % identical: $more90\t" .  "(" . ($more90/11132*100) . " %)\n";
print $fh2 "Number of sequences 80-90 % identical: $more80\t" .  "(" . ($more80/11132*100) . " %)\n";
my $combined = $more90+$identicalscount;
print $fh2 "Number of sequences > 90 % identical: $combined\t" .  "(" . ($combined/11132*100) . " %)\n";
$combined = $more80+$more90+$identicalscount;
print $fh2 "Number of sequences > 80 % identical: $combined\t" .  "(" . ($combined/11132*100) . " %)\n";
print $fh2 "Number of sequences < 10 % different in length: $pless10\t" . "(" . ($pless10/11132*100) . " %)\n";
my $finaltally = $more90less10+$identicalscount;
print $fh2 "Number of sequences > 90 % identical < 10 % difference in length: $finaltally\t" .  "(" . ($finaltally/11132*100) . " %)\n";
print $fh2 "\nIdentical gene models:\n";
for my $i (keys %identicals){
    my @array = @{$identicals{$i}};
    print $fh2 join ("\t", @array) . "\n";}
print $fh2 "\nGene models increased in length:\n";
for my $i (keys %increased){
    my @array = @{$increased{$i}};
    print $fh2 join ("\t", @array) . "\n";}
print $fh2 "\nGene models decreased in length:\n";
for my $i (keys %decreased){
    my @array = @{$decreased{$i}};
    print $fh2 join ("\t", @array) . "\n";}
close $fh2;
