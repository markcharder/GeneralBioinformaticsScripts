#!/usr/bin/env perl

# This script takes a comma-separated list of cufflinks files and combines them so
# that they become lists of transcripts with a single FPKM value per condition.
# The order of the cufflinks files specification is important. From cufflinks, these files
# are typically named genes_fpkm.tracking.

BEGIN{
  our $start_time = time();
}
my $time = localtime;
warn("Script started at $time\n");

use strict;
use warnings;
use Getopt::Long;

my %names_hash;
my $usage=<<USAGE;

$0 --cufflinks <cufflinks fpkm count files comma separated>

Optional:
--base_directory <directory where files are stored>
--out <output file>

USAGE

my $cufflinks_files = "";
my $basedir = "./";
my $output = "combined_genes.fpkm_tracking";

GetOptions('cufflinks=s' => \$cufflinks_files,
           'base_directory:s' => \$basedir,
           'out:s' => \$output);

my @cufflinks_files = split ",", $cufflinks_files;

chdir $basedir or die("\nError: Directory $basedir does not exist\n\n$usage\n");

if($cufflinks_files eq ""){
  die("\nError: No cufflinks FPKM files specified\n\n$usage\n");
}

sub make_hash {
  my $input = shift;
  warn("\nFile $input being added\n");
  open my $fh, "<", $input or die ("\nError: File, $input, does not exist\n\n$usage\n");
  while(my $line = <$fh>){
    if ($. == 1){
      next;
    }
    chomp $line;
    my @fields = split "\t", $line;
    push @{$names_hash{$fields[0]}}, $fields[9];
  }
  close $fh;
  warn("\nFile $input successfully added\n");
  return %names_hash;
  continue;
}

for (my $i = 0; $i < @cufflinks_files; $i++){
  make_hash($cufflinks_files[$i]);
}

for my $key (keys %names_hash){
  my @fpkm_array = @{$names_hash{$key}};
  open my $outfile, ">>", $output;
  print $outfile "$key\t" . join("\t", @fpkm_array) . "\n";
  close $outfile;
}
