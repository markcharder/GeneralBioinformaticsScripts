#/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;

#############################################################################################################################
#
#MAKE SURE VARIABLE $rnaFoldPath IS SET TO THE FULL PATH TO THE RNAfold EXECUTABLE FROM VIENNA PACKAGE (http://www.tbi.univie.ac.at/RNA/)
#OTHER SYSTEM CALLS INCORPORATE MODULES FROM hhmmir (http://biodev.hgen.pitt.edu/kadriAPBC2009.html). $fold_regions, $extractHairpins AND $hhmmir MUST POINT TO THE EXECUTABLE 'fold_regions.pl' AND THE JAR FILES 'ExtractHairpins.jar' AND 'HHMMIR.jar' FROM hhmmir, RESPECTIVELY.
#
#############################################################################################################################

my $extractFolds="fold_regions.pl";
my $extractHairpins="ExtractHairpins.jar";
my $hhmmir="HHMMiR.jar";
my $rnaFoldPath="/usr/local/bin/RNAfold";
my $usage = <<USAGE;

Usage:

--fasta <fasta file> --sam <sam file> --direction <left/right/both> --flank <size of flank(s)> --stack <read stack height> --output <output file>

Refer to README for detailed descriptions of specifications

USAGE
my $input=@ARGV;
my ($fasta, $sam, $direction, $flank, $stack, $output);
GetOptions ("fasta=s"=>\$fasta,
	    "sam=s"=>\$sam,
	    "direction=s"=>\$direction,
	    "flank=i"=>\$flank,
	    "stack=i"=>\$stack, 
	    "output=s"=>\$output) or die ($usage);
die ($usage) unless $input == 12;

open my $samFile, '<', $sam or die "\nPlease check that specified sam file exists\n $usage";
die "\nPlease check sam input is a valid sam file" unless $sam =~ /(\.sam)$/;
my $begin;
my $finish;
my %samin;
my $end;
while (<$samFile>){
	next if /^@/;
	my ($readid,$flag,$fastaref,$start,undef,undef,undef,undef,undef,$sequence) = split;
	my ($readname,$readstack) = split "-",$readid;
	next if $readstack < $stack;
	next if $flag == 4;
	my $readlength=length($sequence);
	$end = ($start + $readlength);
	if ($direction eq "both"){$begin = ($start - $flank); $finish =($end+$flank)}elsif
	($direction eq "left"){$begin = ($start - $flank); $finish = ($begin+$flank+$readlength)} elsif
	($direction eq "right"){$finish = ($end + $flank); $begin = ($finish-$flank-$readlength)};
	push @{$samin{$fastaref}},[$readname,$begin,$finish,$readstack,$start];

}
warn "Sam file parsed\n";

open my $fastaIn, '<', $fasta or die "\nPlease check that specified fasta file exists\n $usage";  
die "\nPlease check fasta input is a valid fasta file" unless $fasta =~ /(\.fasta)$|(\.fa)$/;

my $wholeseq = Bio::SeqIO->new(-file=>$fasta,
			       -format=>'fasta');
while (my $seqobj = $wholeseq->next_seq){
	my $seqid = $seqobj->id;
    	unless(defined $samin{$seqid}){print "check sam matches fasta\n";
	next;	
}

my $seq=$seqobj->seq;
for (@{$samin{$seqid}}){my ($readname, $begin, $finish,$readstack,$start) = @$_;
	warn join("\t",$readname,$begin,$finish,"$readstack...extracted"),"\n";
	open(my $outfile, ">>", "$fasta.extracts");
    print $outfile ">$seqid:$readname:$readstack:".
      "$begin-$finish:$start\n";
    print $outfile substr($seq, $begin, $finish-$begin+1), "\n";
	close $outfile;
  }

}
warn "Genomic regions of flank size $flank with more than $stack identical read mappings have been extracted\n";
my $fileInput = "$fasta.extracts";
system("perl $extractFolds $rnaFoldPath $fileInput");
warn "Folds in genomic extracts determined\n";
system("java -jar $extractHairpins folded.fa hairpins.fa >> removed.txt");
warn "Single-loop hairpin folds extracted\n";
system("java -jar $hhmmir models/posModel.txt models/negModel.txt hairpins.fa $output 0.71");
warn "Hairpins verified using hhmmir for extracted sequences!\n";

