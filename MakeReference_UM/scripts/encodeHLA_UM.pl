#!/usr/bin/perl

##################################################################################################
#
# Sherman Jia, 2012
# encodeHLA.pl
#
# This script encodes HLA alleles into bi-allelic markers (for imputation and PLINK analysis)
# The input ped file should contain: FID,IID,pID,mID,SEX,PHENO, 
#                                    2 each of: HLA-A,B,C,DPA1,DPB1,DQA1,DQB1,DRB1,DRB3,DRB4,DRB5
#
# encodeHLA_UM.pl
#
# As noted in several places below, script modified in November 2017 by Phil Stuart to use build 
# 37 instead of build 36 coordinates, handle modern field-based HLA nomenclature, and process 3 
# additional HLA genes (DRB3, DRB4, DRB5).  Script further modifed in April 2021 by Phil Stuart
# to allow use of either hg19 or hg38 coordinates depending upon value of third argument.
##################################################################################################

if (scalar(@ARGV) !=3){ die "usage: ./encodeHLA.pl HLApedFile outputmapname hg_version_humanref > newped\n";}

$pedname=shift;
$mapname=shift;
# third argument to allow user to select either hg19 or hg39 coordinates
$hg=shift;

my $nsamples = 0;
my %alleles = ();
my %genepos = ();

# all gene positions lifted over from hg18 to either user-selected hg19 or hg38 (P. Stuart)
# positions added for HLA-DRB3,4,5 genes (11-28-2017 v3, P. Stuart)'

if ($hg eq "hg19") {
    $genepos{"HLA_A"} =    29911991;
    $genepos{"HLA_C"} =    31238192;
    $genepos{"HLA_B"} =    31323293;
    $genepos{"HLA_DRB3"} = 32461808;
    $genepos{"HLA_DRB4"} = 32476733;
    $genepos{"HLA_DRB5"} = 32491592;
    $genepos{"HLA_DRB1"} = 32552064;
    $genepos{"HLA_DQA1"} = 32608306;
    $genepos{"HLA_DQB1"} = 32631061;
    $genepos{"HLA_DPA1"} = 33037086;
    $genepos{"HLA_DPB1"} = 33049368;
} elsif ($hg eq "hg38") {
    $genepos{"HLA_A"} =    29944214;
    $genepos{"HLA_C"} =    31270415;
    $genepos{"HLA_B"} =    31355516;
    $genepos{"HLA_DRB3"} = 32494031;
    $genepos{"HLA_DRB4"} = 32508956;
    $genepos{"HLA_DRB5"} = 32523815;
    $genepos{"HLA_DRB1"} = 32584287;
    $genepos{"HLA_DQA1"} = 32640529;
    $genepos{"HLA_DQB1"} = 32663284;
    $genepos{"HLA_DPA1"} = 33069309;
    $genepos{"HLA_DPB1"} = 33081591;
} else {
    die "invalid hg version of human reference specified";
}

open(FILE,"$pedname") || die "can't open $pedname\n";
while (<FILE>) {
    #Adds each HLA allele to the hash table
    chomp;
    s/^\s+//;
    @fields = split(/\s+/); 

    # add lines for HLA-DRB3,4,5 (P. Stuart, v3 11-28-2017)
    
    $alleles{ "HLA_A_".$fields[6] } = $genepos{"HLA_A"};     $alleles{ "HLA_A_".$fields[7] } = $genepos{"HLA_A"};
    $alleles{ "HLA_B_".$fields[8] } = $genepos{"HLA_B"};     $alleles{ "HLA_B_".$fields[9] } = $genepos{"HLA_B"};
    $alleles{ "HLA_C_".$fields[10] } = $genepos{"HLA_C"};    $alleles{ "HLA_C_".$fields[11] } = $genepos{"HLA_C"};
    $alleles{ "HLA_DPA1_".$fields[12] } = $genepos{"HLA_DPA1"}; $alleles{ "HLA_DPA1_".$fields[13] } = $genepos{"HLA_DPA1"};
    $alleles{ "HLA_DPB1_".$fields[14] } = $genepos{"HLA_DPB1"}; $alleles{ "HLA_DPB1_".$fields[15] } = $genepos{"HLA_DPB1"};
    $alleles{ "HLA_DQA1_".$fields[16] } = $genepos{"HLA_DQA1"}; $alleles{ "HLA_DQA1_".$fields[17] } = $genepos{"HLA_DQA1"};
    $alleles{ "HLA_DQB1_".$fields[18] } = $genepos{"HLA_DQB1"}; $alleles{ "HLA_DQB1_".$fields[19] } = $genepos{"HLA_DQB1"};
    $alleles{ "HLA_DRB1_".$fields[20] } = $genepos{"HLA_DRB1"}; $alleles{ "HLA_DRB1_".$fields[21] } = $genepos{"HLA_DRB1"};
    $alleles{ "HLA_DRB3_".$fields[22] } = $genepos{"HLA_DRB3"}; $alleles{ "HLA_DRB3_".$fields[23] } = $genepos{"HLA_DRB3"};
    $alleles{ "HLA_DRB4_".$fields[24] } = $genepos{"HLA_DRB4"}; $alleles{ "HLA_DRB4_".$fields[25] } = $genepos{"HLA_DRB4"};
    $alleles{ "HLA_DRB5_".$fields[26] } = $genepos{"HLA_DRB5"}; $alleles{ "HLA_DRB5_".$fields[27] } = $genepos{"HLA_DRB5"};
    
    # all lines using dfields1 and dfields2 are substituted for original substr($fields) lines that follow;
    # also added lines for HLA-DRB3,4,5 (P. Stuart, v3 11-28-2017)
   
    @dfields1 = ();  @dfields2 = ();
    @dfields1 = split(/:/,$fields[6]);  @dfields2 = split(/:/,$fields[7]);
    $alleles{ "HLA_A_".$dfields1[0] } = $genepos{"HLA_A"}; $alleles{ "HLA_A_".$dfields2[0] } = $genepos{"HLA_A"};

    @dfields1 = ();  @dfields2 = ();
    @dfields1 = split(/:/,$fields[8]);  @dfields2 = split(/:/,$fields[9]);
    $alleles{ "HLA_B_".$dfields1[0] } = $genepos{"HLA_B"}; $alleles{ "HLA_B_".$dfields2[0] } = $genepos{"HLA_B"};

    @dfields1 = ();  @dfields2 = ();
    @dfields1 = split(/:/,$fields[10]);  @dfields2 = split(/:/,$fields[11]);
    $alleles{ "HLA_C_".$dfields1[0] } = $genepos{"HLA_C"}; $alleles{ "HLA_C_".$dfields2[0] } = $genepos{"HLA_C"};

    @dfields1 = ();  @dfields2 = ();
    @dfields1 = split(/:/,$fields[12]);  @dfields2 = split(/:/,$fields[13]);
    $alleles{ "HLA_DPA1_".$dfields1[0] } = $genepos{"HLA_DPA1"}; $alleles{ "HLA_DPA1_".$dfields2[0] } = $genepos{"HLA_DPA1"};

    @dfields1 = ();  @dfields2 = ();
    @dfields1 = split(/:/,$fields[14]);  @dfields2 = split(/:/,$fields[15]);
    $alleles{ "HLA_DPB1_".$dfields1[0] } = $genepos{"HLA_DPB1"}; $alleles{ "HLA_DPB1_".$dfields2[0] } = $genepos{"HLA_DPB1"};

    @dfields1 = ();  @dfields2 = ();
    @dfields1 = split(/:/,$fields[16]);  @dfields2 = split(/:/,$fields[17]);
    $alleles{ "HLA_DQA1_".$dfields1[0] } = $genepos{"HLA_DQA1"}; $alleles{ "HLA_DQA1_".$dfields2[0] } = $genepos{"HLA_DQA1"};

    @dfields1 = ();  @dfields2 = ();
    @dfields1 = split(/:/,$fields[18]);  @dfields2 = split(/:/,$fields[19]);
    $alleles{ "HLA_DQB1_".$dfields1[0] } = $genepos{"HLA_DQB1"}; $alleles{ "HLA_DQB1_".$dfields2[0] } = $genepos{"HLA_DQB1"};

    @dfields1 = ();  @dfields2 = ();
    @dfields1 = split(/:/,$fields[20]);  @dfields2 = split(/:/,$fields[21]);
    $alleles{ "HLA_DRB1_".$dfields1[0] } = $genepos{"HLA_DRB1"}; $alleles{ "HLA_DRB1_".$dfields2[0] } = $genepos{"HLA_DRB1"};

    @dfields1 = ();  @dfields2 = ();
    @dfields1 = split(/:/,$fields[22]);  @dfields2 = split(/:/,$fields[23]);
    $alleles{ "HLA_DRB3_".$dfields1[0] } = $genepos{"HLA_DRB3"}; $alleles{ "HLA_DRB3_".$dfields2[0] } = $genepos{"HLA_DRB3"};

    @dfields1 = ();  @dfields2 = ();
    @dfields1 = split(/:/,$fields[24]);  @dfields2 = split(/:/,$fields[25]);
    $alleles{ "HLA_DRB4_".$dfields1[0] } = $genepos{"HLA_DRB4"}; $alleles{ "HLA_DRB4_".$dfields2[0] } = $genepos{"HLA_DRB4"};

    @dfields1 = ();  @dfields2 = ();
    @dfields1 = split(/:/,$fields[26]);  @dfields2 = split(/:/,$fields[27]);
    $alleles{ "HLA_DRB5_".$dfields1[0] } = $genepos{"HLA_DRB5"}; $alleles{ "HLA_DRB5_".$dfields2[0] } = $genepos{"HLA_DRB5"};

    # $alleles{ "HLA_A_".substr($fields[6],0,2) } = $genepos{"HLA_A"};     $alleles{ "HLA_A_".substr($fields[7],0,2) } = $genepos{"HLA_A"};
    # $alleles{ "HLA_B_".substr($fields[8],0,2) } = $genepos{"HLA_B"};     $alleles{ "HLA_B_".substr($fields[9],0,2) } = $genepos{"HLA_B"};
    # $alleles{ "HLA_C_".substr($fields[10],0,2) } = $genepos{"HLA_C"};    $alleles{ "HLA_C_".substr($fields[11],0,2) } = $genepos{"HLA_C"};
    # $alleles{ "HLA_DPA1_".substr($fields[12],0,2) } = $genepos{"HLA_DPA1"}; $alleles{ "HLA_DPA1_".substr($fields[13],0,2) } = $genepos{"HLA_DPA1"};
    # $alleles{ "HLA_DPB1_".substr($fields[14],0,2) } = $genepos{"HLA_DPB1"}; $alleles{ "HLA_DPB1_".substr($fields[15],0,2) } = $genepos{"HLA_DPB1"};
    # $alleles{ "HLA_DQA1_".substr($fields[16],0,2) } = $genepos{"HLA_DQA1"}; $alleles{ "HLA_DQA1_".substr($fields[17],0,2) } = $genepos{"HLA_DQA1"};
    # $alleles{ "HLA_DQB1_".substr($fields[18],0,2) } = $genepos{"HLA_DQB1"}; $alleles{ "HLA_DQB1_".substr($fields[19],0,2) } = $genepos{"HLA_DQB1"};
    # $alleles{ "HLA_DRB1_".substr($fields[20],0,2) } = $genepos{"HLA_DRB1"}; $alleles{ "HLA_DRB1_".substr($fields[21],0,2) } = $genepos{"HLA_DRB1"};

    $nsamples++;
}
close(FILE);
#print STDERR "$nsamples individuals read\n";


my @all_alleles = sort {
    if (int($alleles{$a}) == int($alleles{$b})){
	return $a cmp $b;
    }else{
	return int($alleles{$a}) <=> int($alleles{$b});
    }
} keys %alleles;
my $al = "";

#Writes unique alleles to map file
open(F2,">$mapname") || die "can't open $mapname\n";
foreach my $allele ( @all_alleles ){
    @tmp = split('\_',$allele);
    $al=$tmp[2];
    if ($al ne "NA" && $al ne "" && $al ne "0" && $al ne "0 0") {
	print F2 "6\t" . $allele . "\t0\t" . $genepos{$tmp[0]."_".$tmp[1]} . "\n";
    }
}
close(F2);

open(FILE,"$pedname") || die "can't open $pedname\n";
while (<FILE>) {
    chomp;
    s/^\s+//;
    @fields = split(/\s+/);
    print "$fields[0]\t$fields[1]\t$fields[2]\t$fields[3]\t$fields[4]\t$fields[5]";

    # add lines for HLA-DRB3,4,5 (P. Stuart, v3 11-28-2017)

    &PrintGenotypes("A",$fields[6],$fields[7]);
    &PrintGenotypes("C",$fields[10],$fields[11]);
    &PrintGenotypes("B",$fields[8],$fields[9]);
    &PrintGenotypes("DRB3",$fields[22],$fields[23]);
    &PrintGenotypes("DRB4",$fields[24],$fields[25]);
    &PrintGenotypes("DRB5",$fields[26],$fields[27]);
    &PrintGenotypes("DRB1",$fields[20],$fields[21]);
    &PrintGenotypes("DQA1",$fields[16],$fields[17]);
    &PrintGenotypes("DQB1",$fields[18],$fields[19]);
    &PrintGenotypes("DPA1",$fields[12],$fields[13]);
    &PrintGenotypes("DPB1",$fields[14],$fields[15]);

    print "\n";
}
close(FILE);
# This subroutine is never called so commented out (P. Stuart)
#sub hashValueAscendingNum {
#    my @all_alleles = sort {$alleles{$a}.$a cmp $alleles{$b}.$b} keys %alleles;
#}

use strict;
sub PrintGenotypes {
    my($locus, $allele1, $allele2) = @_;
    my $al1 = $allele1;
    my $al2 = $allele2;
    $allele1="HLA_".$locus."_".$allele1;
    $allele2="HLA_".$locus."_".$allele2;
    my $G1 = "0";
    my $G2 = "0";
    my @tmp = ();
    my $al = "";
    my $tmplocus = "";
    my @dfields1 = ();  # (P. Stuart addition)
    my @dfields2 = ();  # (P. Stuart addition)
    my @dfields = ();  # (P. Stuart addition)

    foreach my $allele ( @all_alleles ){
	my@tmp = split('\_',$allele);
	$al=$tmp[2];
	$tmplocus=$tmp[0]."_".$tmp[1];
	my@dfields = split(/:/,$al);  # (P. Stuart addition)
	my@dfields1 = split(/:/,$al1);  # (P. Stuart addition)
	my@dfields2 = split(/:/,$al2);  # (P. Stuart addition)
	
        #if allele being checked is valid and is within specific locus
	if ($al ne "NA" && $al ne "" && $al ne "0" && $al ne "0 0" && $tmplocus eq "HLA_".$locus){
	    #if there is no missing types at a locus
	    if ($al1 ne "NA" && $al1 ne "" && $al1 ne "0" && $al2 ne "NA" && $al2 ne "" && $al2 ne "0"){
		#assign genotype for allele 1
		if ($al1 eq "NA" || $al1 eq "" || $al1 eq "0"){
		    $G1 = "0";
		}else{
		    #if alleles matches at 4 digit, or if first 2 digit of data matches allele
		    # if ($al1 eq $al || substr($al1,0,2) eq $al){   # (old code)
		    if ($al1 eq $al || $dfields1[0] eq $al){    # (new code by P. Stuart)
			$G1 = "P";
		    #if data is 2 digit and allele is 4 and they match at 2, set as "unknown"
		    # }elsif (length($al1) == 2 && length($al) == 4 && substr($al,0,2) eq $al1){   # (old code)
		    }elsif ($dfields1[1] eq "" && $dfields[1] ne "" && $dfields[0] eq $al1){     # (new code by P. Stuart)
			$G1 = "0";
		    }else{
			$G1 = "A";
		    }
		}

		#assign genotype for allele 2
		if ($al2 eq "NA" || $al2 eq "" || $al2 eq "0"){
		    $G2 = "0";
		}else{
		    #if alleles matches at 4 digit, or if first 2 digit of data matches allele
		    # if ($al2 eq $al || substr($al2,0,2) eq $al){    # (old code)
		    if ($al2 eq $al || $dfields2[0] eq $al){    # (new code by P. Stuart)
			$G2 = "P";
		    #if data is 2 digit and allele is 4 and they match at 2, set as "unknown"
		    # }elsif (length($al2) == 2 && length($al) == 4 && substr($al,0,2) eq $al2){    # (old code)
		    }elsif ($dfields2[1] eq "" && $dfields[1] ne "" && $dfields[0] eq $al2){      # (new code by P. Stuart)
			$G2 = "0";
		    }else{
			$G2 = "A";
		    }
		}
		if ($G1 eq "0" || $G2 eq "0"){
		    print "\t0 0";
		}else{
		    print "\t" . $G1 . " " . $G2;
		}
	    }else{
		#print 0's if alleles are not typed at locus
		print "\t0 0";
	    }
	}	
    }
}
