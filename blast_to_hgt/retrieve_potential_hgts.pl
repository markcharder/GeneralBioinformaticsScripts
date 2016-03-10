#!/usr/bin/env perl

## This script was used to parse results from a blast analysis against the entire NCBI nr AA database using
## all genomic proteins.
## It takes accessions of hits and downloads their taxids from NCBI. It then compares them to lists of taxids in files
## specified in a list of files given by the user.
## If a hit has an unrelated taxid, an e-value above the user-defined threshold (deafault 1e-30), and a rank below the 
## user-defined threshold (default 5), it is printed to file.
## Given a list of taxids from closest to most distant (defined in this explanation as 1 - 10), the script will give hits not
## in 1 but in 2, not in 1 or 2 but in 3, not in 1, 2 or 3 but in 4; and so on. The order of the list is important,
## going from 1 - 10, it considers taxids to be related as such 1<2<3<4<5<6<7<8<9<10, where "<" indicates child to parent
## taxa going from left to right.

use POSIX;
use strict;
use warnings;

warn "Script started " . localtime() . "\n";

## Define different inputs based on command line arguments.
my $blast_file = $ARGV[0];
my $taxid_file = $ARGV[1];
my $threshold_evalue = $ARGV[2];
my $threshold_rank = $ARGV[3];

unless (!defined ($ARGV[2])){
    if ($ARGV[2] !~ /\d+e-\d+/){

        $threshold_evalue = 1e-30;
        warn "Please specify threshold e-value with full notation \(e.g. '1e-30'\)\n";
        warn "e-value set to default threshold of 1e-30\n";
    }

}

if (!defined ($ARGV[2])){
    warn "e-value set to default threshold of 1e-30\n";
    $threshold_evalue = 1e-30;
}

unless (!defined ($ARGV[3])){
    if (($ARGV[3] !~ /\d+/) or ($ARGV[3] >= 10)){
        $threshold_rank = 5;
        warn "Threshold rank should be a number below 10\n";
        warn "Threshold rank set to default of 5\n";
    }

}

if (!defined ($ARGV[3])){
        $threshold_rank = 5;
        warn "Threshold rank set to default of 5\n";
}

my $bash_path = "./";
my $usage = "$0 <blast results file> <taxid list> <threshold e-value> <threshold rank>";

if (@ARGV < 2){
    die ("\n\n$usage\n\n")
}

warn "Using threshold e-value of $threshold_evalue and threshold rank of $threshold_rank\n";

## Define global variables for collecting during blast table parsing.
my %accessions;
my $query_id;
my $accession;
my $evalue;
my $accession_count=0;
my $score;

warn "Parsing blast table\n";

## Parse blast table and collect a hash of subject ids (uniquefied by adding an incrementing
## number to the beginning of each.
open my $blast_file_input, "<", $blast_file or die ("\n\n$usage\n\n");
while (my $line = <$blast_file_input>){

    if ($. < 23){
        next;

    }

    chomp $line;

    if ($line =~ /Query=/){

        my @query_info = split /\s+/, $line;
        $query_id = $query_info[1];
        $query_id =~ s/mRNA\.//g;
    }

    if ($line =~ /\|/){

        my @hit_accession = split /\|/, $line;
        $accession_count++;
        $accession = $accession_count . $hit_accession[1];
        my @hit_attributes = split /\s+/, $line;
        $evalue = $hit_attributes[-1];
        $score = $hit_attributes[-2];
    }

    unless ($line =~ /^[BLGEQS]|^\s*$|^\s+/){
        push @{$accessions{$accession}}, $query_id, $accession, $evalue, $score;

    }

}

close $blast_file_input;

## Define global variables for use in collecting taxids from the blast table hash.
my $hash_count = 0;
my $previous_query = "";
my @subjects;
my $hash_size = keys %accessions;
my %taxids;
my %temporary_hash;

warn "Collecting taxids from http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=\$f&rettype=fasta&retmode=xml\n";
my $previous_message = 0;
## Using accessions from the blast table subjects, collect the taxids from above link using
## system call to get_taxids.sh. Make sure this script is in the same directory as get_taxids.sh.
for my $accession (sort {my @aa = $accessions{$a}[0] =~ /^([A-za-z]+)(\d+)([A-za-z]+)(\d+)/;
                     my @bb = $accessions{$b}[0] =~ /^([A-za-z]+)(\d+)([A-Za-z]+)(\d+)/; 
                     $aa[1] <=> $bb[1] or $aa[3] <=> $bb[3]} keys %accessions){

    my @accessions_array = @{$accessions{$accession}};
    my $query = $accessions_array[0];
    my $subject = $accessions_array[1];
    my $evalue_hash = $accessions_array[2];
    my $score_hash = $accessions_array[3];
    my $subject_remove_count = $subject;
    $subject_remove_count =~ s/^\d+//g;
    $hash_count++;
    my $percent_complete = ($hash_count / $hash_size) * 100;
    my $message = ceil($percent_complete);

    if (($query eq $previous_query) or ($hash_count == 1)){
        push @subjects, $subject_remove_count;
        push @{$temporary_hash{$subject}}, $subject, $query, $evalue_hash, $score_hash;
    }

    if ($query ne $previous_query) {

        unless ($hash_count == 1){

            my @taxids = `${bash_path}get_taxids.sh @subjects`;

            foreach my $taxid (@taxids){

                my ($original, $taxid_single) = split /\s+/, $taxid;

                for my $key (keys %temporary_hash){

                    my @temporary_array = @{$temporary_hash{$key}};
                    my $test_name = $temporary_array[0];
                    $test_name =~ s/^\d+//g;

                    if ($original eq $test_name){

                        push @{$taxids{$temporary_array[0]}}, $temporary_array[1], $original, 
                        $taxid_single, $temporary_array[2], $temporary_array[3];

                    }

                }

            }

            @subjects = "";
            undef %temporary_hash;

            push @subjects, $subject_remove_count;
            push @{$temporary_hash{$subject}}, $subject, $query, $evalue_hash, $score_hash;
        }

        if ($message != $previous_message){

            warn "$message % complete\n";
        }

        $previous_message = $message;

    }

    if ($hash_count == $hash_size){

        my @taxids = `${bash_path}get_taxids.sh @subjects`;
        my $test = (keys %temporary_hash);

        foreach my $taxid (@taxids){
            my ($original, $taxid_single) = split /\s+/, $taxid;

            for my $key (keys %temporary_hash){

                my @temporary_array = @{$temporary_hash{$key}};
                my $test_name = $temporary_array[0];
                $test_name =~ s/^\d+//g;

                if ($original eq $test_name){

                    push @{$taxids{$temporary_array[0]}},$temporary_array[1], $original, 
                    $taxid_single, $temporary_array[2], $temporary_array[3];

                }

            }

        }

    undef %temporary_hash;

    }
 
    $previous_query = $query;

}

my @taxon_array;
my @names_array;

## Collect taxids from the files specified in the list (2nd argument to this script).
## must be in the same directory, or full path should be specified in the list.

open my $taxon_list, "<", $taxid_file or die ("\n\n$usage\n\n");
while (my $line = <$taxon_list>){
    chomp $line;
    push @names_array, $line;
}
close $taxon_list;

print "Processing in this order:\n" . join ("\n", @names_array) . "\n";

my $i = 0;

foreach my $taxon (@names_array){
    my $file_name = "$taxon";
    open my $file_handle, "<", $file_name or die ("\n\n$usage\n\n");
    while (my $line = <$file_handle>){
        chomp $line;
        $taxon_array[$i]{$line} = "";
    }
    close $file_handle;
    $i++;
}

warn "Finished processing\n";

for my $i (keys $taxon_array[0]) {
    print "$i\n";
}

## Define global variables for use in testing for potential hgt.
$hash_count = 0;
my $hash_length = keys %taxids;
my $test_evalue;
$previous_query = "";
undef %temporary_hash;

## Test whether taxids from the blast hits are defined in the Leoteomycetes.
## If not, test whether e value is above set threshold and the hit is set threshold or above in rank.
## Defaults are e-30 and 5, respectively.
## Print to output file for each taxon.

warn "Testing for potential HGT events\n";

for my $taxid (sort {my @aa = $taxids{$a}[0] =~ /^([A-za-z]+)(\d+)([A-za-z]+)(\d+)/;
                     my @bb = $taxids{$b}[0] =~ /^([A-za-z]+)(\d+)([A-Za-z]+)(\d+)/; 
                     $aa[1] <=> $bb[1] or $aa[3] <=> $bb[3]} keys %taxids){

    my @taxid_array = @{$taxids{$taxid}};
    my $unique_accession = $taxid;
    my $query = $taxid_array[0];
    my $original = $taxid_array[1];
    my $taxid_single = $taxid_array[2];
    my $test_evalue = $taxid_array[3];
    my $test_score = $taxid_array[4];
    $hash_count++;

    if (($query eq $previous_query) or ($hash_count == 1)){

        push @{$temporary_hash{$unique_accession}}, $query, $original, $taxid_single, $test_evalue, $test_score;

    }

    if ($query ne $previous_query){

        unless ($hash_count == 1){

            my $temporary_hash_count = 0;

            for my $temporary_key (sort {$temporary_hash{$a}[3] <=> $temporary_hash{$b}[3] 
                                   or $temporary_hash{$b}[4] <=> $temporary_hash{$a}[4]}
                                   keys %temporary_hash){

            $temporary_hash_count++;
            my @temporary_array = @{$temporary_hash{$temporary_key}};

                if ($temporary_array[3] < $threshold_evalue){

                    if ($temporary_hash_count < $threshold_rank + 1){

                        if ((!defined ($taxon_array[0]{$temporary_array[2]}))
                        and (defined ($taxon_array[1]{$temporary_array[2]}))){
                            open my $out_file, ">>", "$blast_file.$names_array[0].out";
                            print $out_file join ("\t", @temporary_array) . "\n";
                            close $out_file;

                        }

                        if ((!defined ($taxon_array[1]{$temporary_array[2]}))
                        and (defined ($taxon_array[2]{$temporary_array[2]}))){
                            open my $out_file, ">>", "$blast_file.$names_array[1].out";
                            print $out_file join ("\t", @temporary_array) . "\n";
                            close $out_file;

                        }

                        if ((!defined ($taxon_array[2]{$temporary_array[2]}))
                        and (defined ($taxon_array[3]{$temporary_array[2]}))){
                            open my $out_file, ">>", "$blast_file.$names_array[2].out";
                            print $out_file join ("\t", "@temporary_array.out") . "\n";
                            close $out_file;

                        }

                        if ((!defined ($taxon_array[3]{$temporary_array[2]}))
                        and (defined ($taxon_array[4]{$temporary_array[2]}))){
                            open my $out_file, ">>", "$blast_file.$names_array[3].out";
                            print $out_file join ("\t", @temporary_array) . "\n";
                            close $out_file;

                        }

                        if ((!defined ($taxon_array[4]{$temporary_array[2]}))
                        and (defined ($taxon_array[5]{$temporary_array[2]}))){
                            open my $out_file, ">>", "$blast_file.$names_array[4].out";
                            print $out_file join ("\t", @temporary_array) . "\n";
                            close $out_file;

                        }
 
                        if (!defined ($taxon_array[5]{$temporary_array[2]})){
                            open my $out_file, ">>", "$blast_file.$names_array[5].out";
                            print $out_file join ("\t", @temporary_array) . "\n";
                            close $out_file;

                        }

                    }

                }
    
            }

    $temporary_hash_count = 0;
    undef %temporary_hash;
    push @{$temporary_hash{$unique_accession}}, $query, $original, $taxid_single, $test_evalue, $test_score;

        }

    }

    if ($hash_count == $hash_size){

        my $temporary_hash_count = 0;

        for my $temporary_key (sort {$temporary_hash{$a}[3] <=> $temporary_hash{$b}[3]
                               or $temporary_hash{$b}[4] <=> $temporary_hash{$a}[4]}
                               keys %temporary_hash){

        $temporary_hash_count++;
        my @temporary_array = @{$temporary_hash{$temporary_key}};

            if ($temporary_array[3] < $threshold_evalue){

                if ($temporary_hash_count < $threshold_rank + 1){

                    if ((!defined ($taxon_array[0]{$temporary_array[2]}))
                    and (defined ($taxon_array[1]{$temporary_array[2]}))){
                        open my $out_file, ">>", "$blast_file.$names_array[0].out";
                        print $out_file join ("\t", @temporary_array) . "\n";
                        close $out_file;

                    }

                    if ((!defined ($taxon_array[1]{$temporary_array[2]}))
                    and (defined ($taxon_array[2]{$temporary_array[2]}))){
                        open my $out_file, ">>", "$blast_file.$names_array[1].out";
                        print $out_file join ("\t", @temporary_array) . "\n";
                        close $out_file;

                    }

                    if ((!defined ($taxon_array[2]{$temporary_array[2]}))
                    and (defined ($taxon_array[3]{$temporary_array[2]}))){
                        open my $out_file, ">>", "$blast_file.$names_array[2]";
                        print $out_file join ("\t", "@temporary_array.out") . "\n";
                        close $out_file;

                    }

                    if ((!defined ($taxon_array[3]{$temporary_array[2]}))
                    and (defined ($taxon_array[4]{$temporary_array[2]}))){
                        open my $out_file, ">>", "$blast_file.$names_array[3].out";
                        print $out_file join ("\t", @temporary_array) . "\n";
                        close $out_file;

                    }

                    if ((!defined ($taxon_array[4]{$temporary_array[2]}))
                    and (defined ($taxon_array[5]{$temporary_array[2]}))){
                        open my $out_file, ">>", "$blast_file.$names_array[4].out";
                        print $out_file join ("\t", @temporary_array) . "\n";
                        close $out_file;

                    }
 
                    if (!defined ($taxon_array[5]{$temporary_array[2]})){
                        open my $out_file, ">>", "$blast_file.$names_array[5].out";
                        print $out_file join ("\t", @temporary_array) . "\n";
                        close $out_file;

                    }

                }
            
            }

        }

    }
    

    $previous_query = $query;

}

warn "Job finished " . localtime() . "\n";

