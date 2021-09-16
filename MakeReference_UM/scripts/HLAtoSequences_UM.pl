#!/usr/bin/perl
 
#########################################################################################
#
# Author: Sherman Jia, 2012
#
# HLAtoSequences.pl
#
# This script Converts HLA alleles (in .ped file format) to amino acid or DNA sequences 
#
# Input file should contain: FID, IID, pID, mID, sex, pheno, HLA-A (2), B (2), C (2), 
# DPA1 (2), DPB1 (2), DQA1 (2), DQB1 (2), DRB1 (2), DRB3 (2), DRB4 (2), DRB5 (2)
#
#
# HLAtoSequences_UM.pl
#
# As noted in several places below below, script modified in November 2017 by P. Stuart 
# to fully hanlde modern field-based HLA nomenclature and prcoess 3 additional HLA
# genes (HLA-DRB3, -DRB4, -DRB5)
#########################################################################################

if (scalar(@ARGV) !=3){ die "usage: ./HLAtoSequences.pl HLApedFile HLAdictionary TYPE> OutputPedFile\n";}

$pedname=shift;
$dictionary=shift;
$TYPE=shift;

my %ALLELES = ();
my %INSERTIONS = ();
my %LENGTH  = ();
my %INS_LENGTH = ();

# Loads HLA dictionary of amino acid or DNA sequences.
open(FILE,"$dictionary") || die "can't open $dictionary\n";
while (<FILE>) {
    chomp;
    s/^\s+//;
    @fields = split(/\s+/);

    $ALLELES{$fields[0]}=$fields[1];
    $INSERTIONS{$fields[0]}=$fields[2];
    @tmp = split(/:/,$fields[0]);

    $LENGTH{$tmp[0]}=length($fields[1]);
    $INS_LENGTH{$tmp[0]}=length($fields[2]);
}
close(FILE);

open(FILE,"$pedname") || die "can't open $pedname\n";

# Commented out scalar initiation code because these scalars not used (P. Stuart, 11-28-2017)
# my $A1 = ""; my $A2 = "";
# my $B1 = ""; my $B2 = "";
# my $C1 = ""; my $C2 = "";
# my $DPA1 = ""; my $DPA2 = "";
# my $DPB1 = ""; my $DPB2 = "";
# my $DQA1 = ""; my $DQA2 = "";
# my $DQB1 = ""; my $DQB2 = "";
# my $DRB1 = ""; my $DRB2 = "";

# Read input HLA alleles to be converted
while (<FILE>) {
    chomp;
    s/^\s+//;
    @fields = split(/\s+/); 
    # Prints the first 6 columns (FID, IID, pID, mID, sex, pheno)
    print "$fields[0]\t$fields[1]\t$fields[2]\t$fields[3]\t$fields[4]\t$fields[5]";

    # Prints the amino acid / DNA sequence for the pair of alleles at each locus
    # Added lines for HLA-DRB3,4,5 (P. Stuart, 11-28-2017)
    &PrintSequences("A",$fields[6],$fields[7],$TYPE);
    &PrintSequences("C",$fields[10],$fields[11],$TYPE);
    &PrintSequences("B",$fields[8],$fields[9],$TYPE);
    &PrintSequences("DRB3",$fields[22],$fields[23],$TYPE);
    &PrintSequences("DRB4",$fields[24],$fields[25],$TYPE);
    &PrintSequences("DRB5",$fields[26],$fields[27],$TYPE);
    &PrintSequences("DRB1",$fields[20],$fields[21],$TYPE);
    &PrintSequences("DQA1",$fields[16],$fields[17],$TYPE);
    &PrintSequences("DQB1",$fields[18],$fields[19],$TYPE);
    &PrintSequences("DPA1",$fields[12],$fields[13],$TYPE);
    &PrintSequences("DPB1",$fields[14],$fields[15],$TYPE);
    print "\n";
}
close(FILE);

use strict;
# Function to print amino acid / DNA sequences for pairs of HLA alleles
sub PrintSequences {
    my($locus, $allele1, $allele2, $type) = @_;

    # Get allele name in modern nomenclature (original code commented out because it still assummes all field values < 100; P. Stuart)
    # my $al1 = $locus.":".substr($allele1,0,2).":".substr($allele1,2);
    # my $al2 = $locus.":".substr($allele2,0,2).":".substr($allele2,2);
    
    # Get allele name in modern nomenclature (revised code, P. Stuart, accommodates fields with values > 99)
    
    my @tmp = split(/:/,$allele1);
    my $al1 = $locus.":".$tmp[0].":".$tmp[1];
    my @tmp = split(/:/,$allele2);
    my $al2 = $locus.":".$tmp[0].":".$tmp[1];
    
    # Original code assumes that N is only possible 1 letter extension for allele names in modern nomenclature (P. Stuart)
    # if ((exists($ALLELES{$al1}) || exists($ALLELES{$al1."N"})) && (exists($ALLELES{$al2}) || exists($ALLELES{$al2."N"}))){
	# if (exists($ALLELES{$al1."N"})){$al1=$al1."N";}
	# if (exists($ALLELES{$al2."N"})){$al2=$al2."N";}
	
	# Revised code accomodates six possible 1 letter extensions (N,Q,L,S,C,A) (P. Stuart)
	if ((exists($ALLELES{$al1}) || exists($ALLELES{$al1."N"}) || exists($ALLELES{$al1."Q"}) || exists($ALLELES{$al1."L"}) || exists($ALLELES{$al1."S"}) || exists($ALLELES{$al1."C"}) || exists($ALLELES{$al1."A"})) && 
	    (exists($ALLELES{$al2}) || exists($ALLELES{$al2."N"}) || exists($ALLELES{$al2."Q"}) || exists($ALLELES{$al2."L"}) || exists($ALLELES{$al2."S"}) || exists($ALLELES{$al2."C"}) || exists($ALLELES{$al2."A"}))){
	if (exists($ALLELES{$al1."N"})){$al1=$al1."N";}
	if (exists($ALLELES{$al2."N"})){$al2=$al2."N";}
	if (exists($ALLELES{$al1."Q"})){$al1=$al1."Q";}
	if (exists($ALLELES{$al2."Q"})){$al2=$al2."Q";}
	if (exists($ALLELES{$al1."L"})){$al1=$al1."L";}
	if (exists($ALLELES{$al2."L"})){$al2=$al2."L";}
	if (exists($ALLELES{$al1."S"})){$al1=$al1."S";}
	if (exists($ALLELES{$al2."S"})){$al2=$al2."S";}
	if (exists($ALLELES{$al1."C"})){$al1=$al1."C";}
	if (exists($ALLELES{$al2."C"})){$al2=$al2."C";}
	if (exists($ALLELES{$al1."A"})){$al1=$al1."A";}
	if (exists($ALLELES{$al2."A"})){$al2=$al2."A";}
	
	# Print amino acids / DNA sequence
	$allele1 = $ALLELES{$al1}; 
	$allele2 = $ALLELES{$al2};

	if ($type eq "AA"){
	    # If locus is HLA-B, C, DRB1, DPA1, or DQB1, print the amino acids / DNA in reverse order
	    # Include HLA-DRB3,4,5 sequences to be printed in reverse order (P. Stuart, 11-28-2017)
	    if ($locus eq "B" || $locus eq "C" || $locus eq "DRB1" || $locus eq "DRB3" || $locus eq "DRB4" || $locus eq "DRB5" || $locus eq "DPA1" || $locus eq "DQB1"){
		$allele1 = reverse($ALLELES{$al1}); 
		$allele2 = reverse($ALLELES{$al2});
	    }
	}
    }else{
	# If allele is not found in the dictionary, print 0's
	$allele1 = "";
	$allele2 = "";
	for (my $i = 0; $i < $LENGTH{$locus}; $i++){
	    $allele1 = $allele1 . "0";
	    $allele2 = $allele2 . "0";
	}
    }

    # Print amino acid / DNA sequence for pair of alleles
    for (my $i = 0; $i < length($allele1); $i++){
	print "\t" . substr($allele1,$i,1) . " " . substr($allele2,$i,1);
    }

    # Print amino acid insertions for pair of alleles
    my $ins1 = $INSERTIONS{$al1};
    my $ins2 = $INSERTIONS{$al2};
    if ($ins1 ne "" && $ins2 ne ""){
	print "\t" . $ins1 . " " . $ins2;
    }else{
	for (my $i = 0; $i < $INS_LENGTH{$locus}; $i++){
	    print "\t0 0";
	}
    }
    
}
