#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use threads;

my $input1 = "";
my $input2 = "";
my $type = "";
my $parallel = "TRUE";
my $threads = 2;
my $blast1;
my $blast2;
my $usage = <<USAGE;


$0 --input_a <first fasta> --input_b <second fasta> --type <nucl or prot>

Optional:

--parallel <TRUE/FALSE>
--blast_threads <threads per blast>

USAGE

GetOptions('input_a=s' => \$input1, 'input_b=s' => \$input2, 'type=s' => \$type, 'parallel:s' => \$parallel, 'blast_threads:i' => \$threads );

if (($input1 eq "") || ($input2 eq "") || ($type eq "")){
  die("\n\nInvalid options$usage\n");
  }

sub makedb1 {
  system("makeblastdb -dbtype $type -in $input1 -out $input1 >> logfile.txt 2>&1");
  }

sub makedb2{
  system("makeblastdb -dbtype $type -in $input2 -out $input2 >> logfile.txt 2>&1");
  }

sub nucl {
  my $db = shift @_;
  my $in = shift @_;
  system("blastn -num_threads $threads -db $db -query $in -max_target_seqs 1 -max_hsps 1 -outfmt '6 qseqid sseqid evalue pident qlen slen' -out $in.blast >> logfile.txt 2>&1");
  }

sub prot {
  my $db = shift @_;
  my $in = shift @_;
  system("blastp -num_threads $threads -db $db -query $in -max_target_seqs 1 -max_hsps 1 -outfmt '6 qseqid sseqid evalue pident qlen slen' -out $in.blast >> logfile.txt 2>&1");
  }

if($parallel eq "TRUE"){
  my @inarray1 = ($input1, $input2);
  my @inarray2 = ($input2, $input1);
  my $thr1 = threads->create(\&makedb1);
  my $thr2 = threads->create(\&makedb2);
  $thr1->join;
  $thr2->join;
  if($type eq "nucl"){
    my $bthr1 = threads->create(\&nucl, @inarray1);
    open my $fh, ">>", "logfile.txt";
    print $fh "\nThread 1 created\nRunning blast 1 in one process on $threads threads\n";
    close $fh;
    my $bthr2 = threads->create(\&nucl, @inarray2);
    open $fh, ">>", "logfile.txt";
    print $fh "\nThread 2 created\nRunning blast 2 in one process on $threads threads\n";
    close $fh;
    $bthr1->join;
    $bthr2->join;
    }
  if($type eq "prot"){
    my $bthr1 = threads->create(\&prot, @inarray1);
    open my $fh, ">>", "logfile.txt";
    print $fh "\nThread 1 created\nRunning blast 1 in one process on $threads threads\n";
    close $fh;
    my $bthr2 = threads->create(\&prot, @inarray2);
    open $fh, ">>", "logfile.txt";
    print $fh "\nThread 2 created\nRunning blast 2 in one process on $threads threads\n";
    close $fh;
    $bthr1->join;
    $bthr2->join;
    }
  }

if($parallel eq "FALSE"){
  if($type eq "nucl"){
    makedb1;
    makedb2;
    nucl;
    }
  if($type eq "prot"){
    makedb1;
    makedb2;
    prot;
    }
  }

my %hash1;
my %hash2;
my %identhash;
my %qlenhash;
my %slenhash;

open my $fh, "<", "$input1.blast";
while (my $line = <$fh>){
  my @fields = split /\s+/, $line;
  $hash1{$fields[0]} = $fields[1];
  $identhash{$fields[0]} = $fields[3];
  $qlenhash{$fields[0]} = $fields[4];
  $slenhash{$fields[0]} = $fields[5];
  }
close $fh;

open $fh, "<", "$input2.blast";
while (my $line = <$fh>){
  my @fields = split /\s+/, $line;
  $hash2{$fields[0]} = $fields[1];
  }
close $fh;

open $fh, ">>", "reciprocal_blast.txt";
for my $key (keys %hash1){
  my $pdiff = $qlenhash{$key}/$slenhash{$key};
  if (defined($hash2{$hash1{$key}})){
    if ($hash2{$hash1{$key}} eq $key){
      print $fh "$key\t$hash1{$key}\treciprocal\t$identhash{$key}\t$qlenhash{$key}\t$slenhash{$key}\t$pdiff\n";
      }
    else{
      print $fh "$key\t$hash1{$key}\tnot_reciprocal\t$identhash{$key}\t$qlenhash{$key}\t$slenhash{$key}\t$pdiff\n";
      }
    }
  }
close $fh;
