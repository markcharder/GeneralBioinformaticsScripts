#!/usr/bin/env perl

BEGIN{ 
  our $start_time = time();
}
my $time = localtime;
warn("Script started at $time\n");

# Takes as input an FPKM counts table with an arbitrary even number of columns.
# The columns are split into two, and the first half compared against the second half,
# column by column to get genes that have lower/higher FPKM values by at least 50 % in the 
# first compared to the second.
# A threshold FPKM can be set for the sample with the higher fpkm value.

use Getopt::Long;
use strict;
use warnings;

#Define all global variables.
my $fpkm_table = "";
my $upregulated_output = "upregulated.txt";
my $downregulated_output = "downregulated.txt";
my $min_fpkm = 100;

# Define usage.
my $usage = <<USAGE;

Usage:

$0 --fpkm_table <fpkm table> 

Optional:
--upregulated_output <output file of upregulated genes>
--downregulated_output <output file of downregulated genes>
--minimum_fpkm <threshold for higher fpkm value>

USAGE

GetOptions('fpkm_table=s' => \$fpkm_table,
	   'upregulated_output:s' => \$upregulated_output,
	   'downregulated_output:s' => \$downregulated_output,
	   'minimum_fpkm:i' => \$min_fpkm);

# FPKM table is not provided, script dies with usage message.
if ($fpkm_table eq ""){
  die ("\nError: FPKM table not provided\n\n$usage\n");
}

#Open FPKM table and test the conditions defined above.
open my $fh, "<", $fpkm_table or die ("\nError: FPKM table, $fpkm_table, does not exist\n\n$usage\n");
while(my $line = <$fh>){
  if ($. == 1){
    next;
  }
  chomp $line;
  my @fields = split "\t", $line;
  my $id = $fields[0];
  shift @fields;
  my $midpoint = @fields / 2 - 1;
  my @fpkms_test = @fields[0..$midpoint];
  my @fpkms_control = @fields[$midpoint+1..$#fields-1];
  for (my $i = 0; $i < $midpoint; $i++){
    my $addition_factor = $midpoint+1;
    if (($fields[$i] < 0.5 * $fields[$addition_factor]) and ($fields[$addition_factor] > $min_fpkm)){
      open my $output, ">>", $downregulated_output or die ("$usage\n");
      print $output "$line\n";
      close $output;
    }
    if (($fields[$i] > 2 * $fields[$addition_factor]) and ($fields[$i] > $min_fpkm)){
      open my $output, ">>", $upregulated_output or die ("$usage\n");
      print $output "$line\n";
      close $output;
    }
  }
}
close $fh;

my $end_time = localtime;
my $run_time = time() - our $start_time;

# Print run time.
warn("Script ended at $end_time\nTotal run time: $run_time\n");
