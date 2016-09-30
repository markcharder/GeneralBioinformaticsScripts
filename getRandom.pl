#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Switch;

my $chr1 = 3951982;
my $chr2 = 3683506;
my $chr3 = 3351453;
my $chr4 = 2873318;
my $chr5 = 2822964;
my $chr6 = 2483831;
my $chr7 = 2434682;
my $chr8 = 2299506;
my $chr9 = 2122865;
my $chr10 = 2105496;
my $chr11 = 2052242;
my $chr12 = 1878461;
my $chr13 = 1845946;
my $chr14 = 1815632;
my $chr15 = 1765292;
my $chr16 = 1419421;

my $bed = "";
my $fasta = "";

GetOptions('bed=s' => \$bed, 'fasta=s' => \$fasta);
my $usage = <<USAGE;


	$0 --bed <bed file> --fasta <fasta file>


USAGE

if (($bed eq "") or ($fasta eq "")){
  die ("\nIncorrect input\n$usage\n");
  }

my $prefix = $bed; 
$prefix =~ s/\.\S+//g;

open my $fh, "<", $bed or die ("\nFile not readable\n$usage\n");
while (my $line = <$fh>){
  chomp $line;
  my @fields = split /\s+/, $line;
  my $seqSize = $fields[2] - $fields[1];
  my $strand = $fields[5];
  my $start;
  switch ($fields[0]){
    case "Chr_1" {$start = int(rand($chr1))}
    case "Chr_2" {$start = int(rand($chr2))}
    case "Chr_3" {$start = int(rand($chr3))}
    case "Chr_4" {$start = int(rand($chr4))}
    case "Chr_5" {$start = int(rand($chr5))}
    case "Chr_6" {$start = int(rand($chr6))}
    case "Chr_7" {$start = int(rand($chr7))}
    case "Chr_8" {$start = int(rand($chr8))}
    case "Chr_9" {$start = int(rand($chr9))}
    case "Chr_10" {$start = int(rand($chr10))}
    case "Chr_11" {$start = int(rand($chr11))}
    case "Chr_12" {$start = int(rand($chr12))}
    case "Chr_13" {$start = int(rand($chr13))}
    case "Chr_14" {$start = int(rand($chr14))}
    case "Chr_15" {$start = int(rand($chr15))}
    case "Chr_16" {$start = int(rand($chr16))}
    }
  my $end = $start + $seqSize;
  open my $fh1, ">>", "$prefix.randomised.bed";
  print $fh1 "$fields[0]\t$start\t$end\t.\t.\t$strand\n";
  close $fh1;
  }
close $fh;

system("bedtools getfasta -fi $fasta -bed $prefix.bed -fo $prefix.fa");
system("bedtools getfasta -fi $fasta -bed $prefix.randomised.bed -fo $prefix.randomised.fa");
system("sizeseq -descending -sequences $prefix.fa -outseq $prefix.sorted.fa");
system("sizeseq -descending -sequences $prefix.randomised.fa -outseq $prefix.randomised.sorted.fa");

my $count = 0;
open $fh, "<", "$prefix.sorted.fa";
while (my $line = <$fh>){
  if ($count < 51){
    open my $fh1, ">>", "$prefix.top50.fa";
    print $fh1 $line;
    close $fh1;
    }
  if ($line =~ />/){
    $count++;
    }
  }
close $fh;

$count = 0;
open $fh, "<", "$prefix.randomised.sorted.fa";
while (my $line = <$fh>){
  if ($count < 51){
    open my $fh1, ">>", "$prefix.randomised.top50.fa";
    print $fh1 $line;
    close $fh1;
    }
  if ($line =~ />/){
    $count++;
    }
  }
close $fh;


my @seqArray;
my $head;

open $fh, "<", "$prefix.randomised.top50.fa" or die;
while (my $line = <$fh>){
  chomp $line;
  if ($line =~ />/){
    if ($. > 1){
      my $ta = () = join("", @seqArray) =~ /TA/g;
      my $at = () = join("", @seqArray) =~ /AT/g;
      my $ratio = $ta/$at;
      open my $fh1, ">>", "$prefix.top50.index";
      print $fh1 "random_$prefix\t$ratio\n";
      close $fh1;
      }
      $head = $line;
      $head =~ s/>//g;
      @seqArray = ();
    }
  else {
    my @seqs = split "", $line;
    push @seqArray, @seqs;
    if (eof($fh)){
      my $ta = () = join("", @seqArray) =~ /TA/g;
      my $at = () = join("", @seqArray) =~ /AT/g;
      my $ratio = $ta/$at;
      open my $fh1, ">>", "$prefix.top50.index";
      print $fh1 "random_$prefix\t$ratio\n";
      close $fh1;
      }
    }
  }
close $fh;

open $fh, "<", "$prefix.top50.fa" or die ("\nFasta $prefix.fa not readable \n$usage\n");
while (my $line = <$fh>){
  chomp $line;
  if ($line =~ />/){
    if ($. > 1){
      my $ta = () = join("", @seqArray) =~ /TA/g;
      my $at = () = join("", @seqArray) =~ /AT/g;
      my $ratio = $ta/$at;
      open my $fh1, ">>", "$prefix.top50.index";
      print $fh1 "$prefix\t$ratio\n";
      close $fh1;
      }
      $head = $line;
      $head =~ s/>//g;
      @seqArray = ();
    }
  else {
    my @seqs = split "", $line;
    push @seqArray, @seqs;
    if (eof($fh)){
      my $ta = () = join("", @seqArray) =~ /TA/g;
      my $at = () = join("", @seqArray) =~ /AT/g;
      my $ratio = $ta/$at;
      open my $fh1, ">>", "$prefix.top50.index";
      print $fh1 "$prefix\t$ratio\n";
      close $fh1;
      }
    }
  }
close $fh;
