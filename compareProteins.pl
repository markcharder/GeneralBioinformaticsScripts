#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "$0 <blast hits table> <outfile>";

if (@ARGV != 2){die (print "\n$usage\n\n")};

my %identicals;
my %morethanten;
my %lessthanten;
my %lessthan90;
my %morethan90;

open my $blasthits, "<", $ARGV[0] or die (print "\n$usage\n\n");
while (my $line = <$blasthits>){

    chomp $line;
    my @fields = split "\t", $line;
    $fields[0] =~ s/mRNA\.//g;
    my $pdiff;

    if ($fields[3] > $fields[1]){
        $pdiff = $fields[1] / $fields[3];}

    if ($fields[1] > $fields[3]){
        $pdiff = $fields[3] / $fields[1];}

    if ($fields[1] == $fields[3]){
        $pdiff = 0;}

    if ($pdiff < 0.9){
        push @{$morethanten{$fields[0]}}, $fields[0], $fields[2], $fields[5], $fields[6], $pdiff;}

    if ($pdiff > 0.9){
        push @{$lessthanten{$fields[0]}}, $fields[0], $fields[2], $fields[5], $fields[6], $pdiff;}

    if ($fields[6] < 90){
        push @{$lessthan90{$fields[0]}}, $fields[0], $fields[2], $fields[5], $fields[6], $pdiff;}

    if ($fields[6] > 90){
        push @{$morethan90{$fields[0]}}, $fields[0], $fields[2], $fields[5], $fields[6], $pdiff;}

    if (($fields[6] == 100) and ($pdiff == 0)){
        push @{$identicals{$fields[0]}}, $fields[0], $fields[2], $fields[5], $fields[6], $pdiff;}}

close $blasthits;

open my $outfile, ">", $ARGV[1] or die (print "\n$usage\n\n");

my $morethanten = keys %morethanten;
my $lessthanten = keys %lessthanten;
my $lessthan90 = keys %lessthan90;
my $morethan90 = keys %morethan90;
my $identicals = keys %identicals;

my %lessthan90morethan10;
for my $i (sort keys %morethanten){
    my @array = @{$morethanten{$i}};
    if (defined $lessthan90{$array[0]}){
        push @{$lessthan90morethan10{$array[0]}}, $array[0], $array[1], $array[2], $array[3], $array[4];}}

my %morethan90lessthan10;
for my $i (sort keys %lessthanten){
    my @array = @{$lessthanten{$i}};
    if (defined $morethan90{$array[0]}){
        push @{$morethan90lessthan10{$array[0]}}, $array[0], $array[1], $array[2], $array[3], $array[4];}}

my $lessthan90morethan10 = keys %lessthan90morethan10;
my $morethan90lessthan10 = keys %morethan90lessthan10;

print $outfile "Identical sequences: $identicals\n\n" .
               "Length-difference < 10 %: $lessthanten\n\n" .
               "Length-difference < 10 %, identity > 90 %: $morethan90lessthan10\n\n" .
               "Length-difference > 10 %: $morethanten\n\n" .
               "Length-difference > 10 %, identity < 90 %: $lessthan90morethan10\n\n" .
               "\nIdentical sequences:\n" .
               "V2-ID\tV1-ID\te-Value\tIdentity\t%-Difference\n";
for my $i (sort keys %identicals){
    my @array = @{$identicals{$i}};
    print $outfile join ("\t", @array) . "\n";}

print $outfile "\n\nLength-difference < 10 %:\n" .
               "V2-ID\tV1-ID\te-Value\tIdentity\t%-Difference\n";
for my $i (sort keys %lessthanten){
    my @array = @{$lessthanten{$i}};
    print $outfile join ("\t", @array) . "\n";}

print $outfile "\n\nLength-difference < 10 %, identity > 90 %:\n" .
               "V2-ID\tV1-ID\te-Value\tIdentity\t%-Difference\n";
for my $i (sort keys %morethan90lessthan10){
    my @array = @{$morethan90lessthan10{$i}};
    print $outfile join ("\t", @array) . "\n";}

print $outfile "\n\nLength-difference > 10 %:\n" .
               "V2-ID\tV1-ID\te-Value\tIdentity\t%-Difference\n";
for my $i (sort keys %morethanten){
    my @array = @{$morethanten{$i}};
    print $outfile join ("\t", @array) . "\n";}

print $outfile "\n\nLength-difference > 10 %, identity < 90 %:\n" .
               "V2-ID\tV1-ID\te-Value\tIdentity\t%-Difference\n";
for my $i (sort keys %lessthan90morethan10){
    my @array = @{$lessthan90morethan10{$i}};
    print $outfile join ("\t", @array) . "\n";}
close $outfile;
