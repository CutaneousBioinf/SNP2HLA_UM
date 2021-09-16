#!/bin/csh -f

#############################################################################################################################################################################################################################################################
#   SNP2HLA_UM.csh
#
#   DESCRIPTION: This script runs imputation of HLA amino acids and classical alleles using SNP data.
#
#   Author: Sherman Jia (xiaomingjia@gmail.com)
#           + Small modifications by Buhm Han (buhmhan@broadinstitute.org): 8/7/12
#           + Further modifications by Phil Stuart (pstuart@umich.edu):
#             (a) Adapted script to use Beagle 5.2 instead of Beagle 3.0.4 as its imputation
#                 engine, thereby providing substantially more accurate and faster imputation
#                 and phasing.  Verified to work with release 05Apr21.9b7 of Beagle 5.2 and
#                 also with the last version(18May20.d20) of Beagle 5.1.
#             (b) Requires a reference panel created by MakeReference_UM script of this
#                 updated package.
#             (c) Script now uses Plink 1.9 instead of v1.07 because of its improved speed.
#             (d) Values for heap memory, stack memory, number of threads and the number of
#                 burn-in and phasing iterations used by Beagle can optionally be specified
#                 on the command line when launching this script.
#             (e) The command line can now optionally specify the name of a PLINK format
#                 genetic map file with cM units, which can be used by Beagle to enhance the
#                 accuracy of its phasing and imputation.
#             (f) Created parameter MAF_MIN to facilitate changing the minimum allowed MAF
#                 for reference panel SNPs from the current default of 0.01.
#             (g) Script now outputs a single .vcf file holding all imputation results
#                 rather than the suite of eight files output by the older version of this
#                 script.
#
#   IMPORTANT NOTES:
# 
#     (1) First input is Plink dataset (*.bed/bim/fam) containing SNP data for the target sample
#
#     (2) Second input is reference panel (*.bgl.vcf and *.markers in Beagle 5.2 format; 
#         *.fam/.bim/.FRQ.frq in PLINK format)
#
#     (3) SNP2HLA uses ID-based variant matching when comparing target dataset with reference panel.
#         It is the user's responsibility to ensure that the naming scheme for variant IDs in the
#         target dataset is consistent with that used in the reference panel.  If variant IDs are
#         not consistent, the Plink.bim file for the target sample must be modified to achieve
#         variant ID consistency.
#
#     (4) When the target sample and reference panel differ substantially in their ethnic or
#         population makeup, it may be beneficial to increase the TOLERATED_DIFF parameter from
#         its default value of 0.15.  The tolerated difference(after flipping if necessary)
#         between allele frequencies of a variant present in both the target and reference
#         datasets is a QC measure that presumably weeds out variants that are poorly typed
#         in either dataset, but the default value may needlessly eliminate variants well-typed
#         in both datasets when there are bona fide discrepancies in variant allele frequecies
#         between the two source populations.
# 
#   DEPENDENCIES: (download and place in the same folder as this script)
#     1. PLINK (1.9); available at https://www.cog-genomics.org/plink/1.9/
#     2. Beagle (5.2); available at https://faculty.washington.edu/browning/beagle/beagle.html
#     3. merge_tables.pl (Perl script to merge files indexed by a specific column)
#     4. hackvcf (utility by Phil Stuart to hack alleles in .vcf reference panel before input into Beagle 5.2)
#     5. revert_alleles (utility by Phil Stuart to revert hacked non-ACTG alleles in *.vcf imputed output file, necessitated by Beagle 5.2, back to allele identities in original reference panel)
#     6. If genetic_map_file argument specified, PLINK format genetic map on cM scale with build 37 or 38 coordinates (plink.chr6.GRCh37.map or plink.chr6.GRCh38.map), downloaded from http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/
#
#    USAGE: ./SNP2HLA_UM.csh DATA (.bed/.bim/.fam) REFERENCE (.bgl.phased/.markers/.fam/.bim/.FRQ.frq) OUTPUT plink {optional: java_max_memory[gb] java_max_stacksize[mb] nthreads burnin niterations genetic_map_file
#############################################################################################################################################################################################################################################################

if ($#argv < 4) then
    echo "USAGE: ./SNP2HLA_UM.csh DATA (.bed/.bim/.fam) REFERENCE (.bgl.vcf/.markers/.fam/.bim/.FRQ.frq) OUTPUT plink {optional: java_max_memory[gb] java_max_stacksize[mb] number_of_threads number_of_burnin_iterations number_of_phasing_iterations plink_format_genetic_map_file}"; exit 1
endif

set SCRIPTPATH=`dirname $0`
set MERGE=$SCRIPTPATH/merge_tables.pl

# CHECK FOR DEPENDENCIES
if (! -e $SCRIPTPATH/$4) then
    echo "Please copy the executable of version 1.9 of PLINK (https://www.cog-genomics.org/plink/1.9/) into the current directory."; exit 1 
else if (! -e $SCRIPTPATH/beagle.jar) then
    echo "Please install Beagle 5.2 (https://faculty.washington.edu/browning/beagle/beagle.html), rename it as beagle.jar and copy it into the current directory."; exit 1
else if (! -e $SCRIPTPATH/hack_vcf) then
  echo "Please copy hack_vcf (included with this package) into the current directory."; exit 1
else if (! -e $MERGE) then
    echo "Please copy merge_tables.pl (included with this package) into the current directory."; exit 1
else if (! -e $SCRIPTPATH/revert_alleles) then
    echo "Please copy revert_alleles (included with this package) into the current directory."; exit 1
else if ($#argv >= 10 && ! -e $SCRIPTPATH/$10) then
  echo "Please download appropriate Plink-format genetic recombination map from http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/ and copy it into into this directory."; exit 1
endif

# INPUTS
set INPUT=$1
set REFERENCE=$2
set OUTPUT=$3
set PLINK=$4

if ($#argv >= 5) then
    set MEM=$5
	set PLINK_MEM="$5"000""
else
    set MEM=2 # Default java memory 2Gb
	set PLINK_MEM=2000
endif

if ($#argv >= 6) then
    set SS=$6
else
    set SS=5 # Default java stacksize 5 Mb
endif

if ($#argv >= 7) then
    set THREAD=$7
else
    set THREAD=1
endif

if ($#argv >= 8) then
    set BURNIN=$8
else
    set BURNIN=6
endif

if ($#argv >= 9) then
    set ITER=$9
else
    set ITER=12
endif

if ($#argv >= 10) then
    set MAP=$10
endif

#CHECK FOR PRESENCE OF INPUT FILES IN DIRECTORY

if (! -f $INPUT.bed) then
  echo "Input file $INPUT.bed not found"; exit 1
else if (! -f $INPUT.bim) then
  echo "Input file $INPUT.bim not found"; exit 1
else if (! -f $INPUT.fam) then
  echo "Input file $INPUT.fam not found"; exit 1
else if (! -f $REFERENCE.bim) then
  echo "Input file $REFERENCE.bim not found"; exit 1
else if (! -f $REFERENCE.markers) then
  echo "Input file $REFERENCE.markers not found"; exit 1
else if (! -f $REFERENCE.fam) then
  echo "Input file $REFERENCE.fam not found"; exit 1
else if (! -f $REFERENCE.bgl.vcf) then
  echo "Input file $REFERENCE.bgl.vcf not found"; exit 1
else if (! -f $REFERENCE.FRQ.frq) then
  echo "Input file $REFERENCE.FRQ.frq not found"; exit 1
endif

set JAVATMP=$OUTPUT.javatmpdir
mkdir -p $JAVATMP
alias plink './$PLINK --noweb --silent --allow-no-sex --memory $PLINK_MEM --threads $THREAD'
alias beagle 'java -Djava.io.tmpdir=$JAVATMP -Xss$SS\m -Xmx$MEM\g -jar $SCRIPTPATH/beagle.jar'
alias hack_vcf './hack_vcf -chr 6 -revert_opt 1'
alias vcf2phased './vcf2phased -gprobs_opt 1 -beagle_opt 2'
alias revert_alleles './revert_alleles -script_opt 2'

# Functions to run
set CONVERT_IN  = 1
set EXTRACT_MHC = 1
set FLIP        = 1
set IMPUTE      = 1
set CONVERT_OUT = 1
set CLEANUP     = 1

# SET PARAMETERS
set TOLERATED_DIFF = .15
set MAF_MIN = .01
set i = 1

echo "SNP2HLA_UM: Performing HLA imputation for dataset $INPUT";
echo "- Java memory = "$MEM"Gb"
echo "- Java stacksize = "$SS"Mb"
echo "- Beagle number of threads = "$THREAD" threads"
echo "- Beagle number of burnin iterations = "$BURNIN" iterations"
echo "- Beagle number of phasing iterations = "$ITER" iterations"
if ($#argv >= 10) then
    echo "- Beagle genetic map file = "$MAP""
else
    echo "- Beagle genetic map file = none specified"
endif

# Hack non-standard allele identities in .vcf reference panel to a form permitted by Beagle 5.2

if ($CONVERT_IN) then
  echo ""
  echo "[$i] Hacking non-standard allele identities in .vcf reference panel to a form permitted by Beagle 5.2."; @ i++
  hack_vcf -fnmarker $REFERENCE.markers -fnfam $REFERENCE.fam -fnvcf_in $REFERENCE.bgl.vcf -fnvcf_out $REFERENCE.vcf
endif

set MHC=$OUTPUT.MHC

if ($EXTRACT_MHC) then
    echo "[$i] Extracting SNPs from the MHC."; @ i++
    plink --bfile $INPUT --chr 6 --from-mb 29 --to-mb 34 --maf $MAF_MIN --make-bed --out $OUTPUT.MHC
endif
	
if ($FLIP) then
    echo ""
    echo "[$i] Performing SNP quality control."; @ i++

    # Identifying non-A/T non-C/G SNPs to flip
    echo "SNP 	POS	A1	A2" > $OUTPUT.tmp1
    cut -f2,4- $MHC.bim >> $OUTPUT.tmp1
    echo "SNP 	POSR	A1R	A2R" > $OUTPUT.tmp2
    cut -f2,4- $REFERENCE.bim >> $OUTPUT.tmp2
    $MERGE $OUTPUT.tmp2 $OUTPUT.tmp1 SNP |  grep -v -w NA > $OUTPUT.SNPS.alleles

    awk '{if ($3 != $6 && $3 != $7){print $1}}' $OUTPUT.SNPS.alleles > $OUTPUT.SNPS.toflip1
    plink --bfile $MHC --flip $OUTPUT.SNPS.toflip1 --make-bed --out $MHC.FLP

    # Calculating allele frequencies (sed command was altered by P. Stuart to enforce replacement of only exact matches to A1, A2, MAF)
    plink --bfile $MHC.FLP --freq --out $MHC.FLP.FRQ
    sed 's/\<A1\>/A1I/g' $MHC.FLP.FRQ.frq | sed 's/\<A2\>/A2I/g' | sed 's/\<MAF\>/MAF_I/g' > $OUTPUT.tmp

    mv $OUTPUT.tmp $MHC.FLP.FRQ
    $MERGE $REFERENCE.FRQ.frq $MHC.FLP.FRQ.frq SNP | grep -v -w NA > $OUTPUT.SNPS.frq
    sed 's/ /\t/g' $OUTPUT.SNPS.frq | awk '{if ($3 != $8){print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $9 "\t" $8 "\t" 1-$10 "\t*"}else{print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $8 "\t" $9 "\t" $10 "\t."}}' > $OUTPUT.SNPS.frq.parsed
    
    # Finding A/T and C/G SNPs
    awk '{if (($2 == "A" && $3 == "T") || ($2 == "T" && $3 == "A") || ($2 == "C" && $3 == "G") || ($2 == "G" && $3 == "C")){if ($4 > $7){diff=$4 - $7; if ($4 > 1-$7){corrected=$4-(1-$7)}else{corrected=(1-$7)-$4}}else{diff=$7-$4;if($7 > (1-$4)){corrected=$7-(1-$4)}else{corrected=(1-$4)-$7}};print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" diff "\t" corrected}}' $OUTPUT.SNPS.frq.parsed > $OUTPUT.SNPS.ATCG.frq

    # Identifying A/T and C/G SNPs to flip or remove
    awk '{if ($10 < $9 && $10 < .15){print $1}}' $OUTPUT.SNPS.ATCG.frq > $OUTPUT.SNPS.toflip2
    awk '{if ($4 > 0.4){print $1}}' $OUTPUT.SNPS.ATCG.frq > $OUTPUT.SNPS.toremove

    # Identifying non A/T and non C/G SNPs to remove
    awk '{if (!(($2 == "A" && $3 == "T") || ($2 == "T" && $3 == "A") || ($2 == "C" && $3 == "G") || ($2 == "G" && $3 == "C"))){if ($4 > $7){diff=$4 - $7;}else{diff=$7-$4}; if (diff > '$TOLERATED_DIFF'){print $1}}}' $OUTPUT.SNPS.frq.parsed >> $OUTPUT.SNPS.toremove
    awk '{if (($2 != "A" && $2 != "C" && $2 != "G" && $2 != "T") || ($3 != "A" && $3 != "C" && $3 != "G" && $3 != "T")){print $1}}' $OUTPUT.SNPS.frq.parsed >> $OUTPUT.SNPS.toremove
    awk '{if (($2 == $5 && $3 != $6) || ($3 == $6 && $2 != $5)){print $1}}' $OUTPUT.SNPS.frq.parsed >> $OUTPUT.SNPS.toremove

    # Making QCd SNP file (sed command was altered by P. Stuart to enforce replacement of only exact matches to A1, A2, MAF)
    plink --bfile $MHC.FLP --geno 0.2 --exclude $OUTPUT.SNPS.toremove --flip $OUTPUT.SNPS.toflip2 --make-bed --out $MHC.QC
    plink --bfile $MHC.QC --freq --out $MHC.QC.FRQ
    sed 's/\<A1\>/A1I/g' $MHC.QC.FRQ.frq | sed 's/\<A2\>/A2I/g' | sed 's/\<MAF\>/MAF_I/g' > $OUTPUT.tmp
    mv $OUTPUT.tmp $MHC.QC.FRQ.frq
    $MERGE $REFERENCE.FRQ.frq $MHC.QC.FRQ.frq SNP | grep -v -w NA > $OUTPUT.SNPS.QC.frq

    cut -f2 $OUTPUT.SNPS.QC.frq | awk '{if (NR > 1){print $1}}' > $OUTPUT.SNPS.toinclude

    echo "SNP 	POS	A1	A2" > $OUTPUT.tmp1
    cut -f2,4- $MHC.QC.bim >> $OUTPUT.tmp1

    $MERGE $OUTPUT.tmp2 $OUTPUT.tmp1 SNP | awk '{if (NR > 1){if ($5 != "NA"){pos=$5}else{pos=$2}; print "6\t" $1 "\t0\t" pos "\t" $3 "\t" $4}}' > $MHC.QC.bim

    # Extracting SNPs and recoding QC'd file as vcf
    plink --bfile $MHC.QC --extract $OUTPUT.SNPS.toinclude --make-bed --out $MHC.QC.reorder
	plink --bfile $MHC.QC.reorder --recode vcf-iid --a1-allele $REFERENCE.markers 4 1 --out $MHC.QC

    # Remove temporary files
    rm $OUTPUT.tmp1 $OUTPUT.tmp2
    rm $MHC.FLP*
    rm $OUTPUT.SNPS.*
endif

if ($IMPUTE) then
    echo ""
    echo "[$i] Performing HLA imputation."; @ i++

	if ($#argv >= 10) then
        beagle ref=$REFERENCE.vcf gt=$MHC.QC.vcf impute=true gp=true nthreads=$THREAD chrom=6 burnin=$BURNIN iterations=$ITER out=$OUTPUT.bgl map=$MAP
	else
        beagle ref=$REFERENCE.vcf gt=$MHC.QC.vcf impute=true gp=true nthreads=$THREAD chrom=6 burnin=$BURNIN iterations=$ITER out=$OUTPUT.bgl
	endif
endif

if ($CONVERT_OUT) then
    gunzip -d $OUTPUT.bgl.vcf.gz
    echo ""
	echo "[$i] Reverting hacked allele identities in *.vcf imputation output file back to those used by reference panel."; @ i++
	revert_alleles -fnvcf $OUTPUT.bgl.vcf -fnref $REFERENCE.bim 
endif

if ($CLEANUP) then
    rm $REFERENCE.vcf 
    rm $OUTPUT.MHC.* >& /dev/null
    rm $OUTPUT.tmp* >& /dev/null
    rm $OUTPUT.IMPUTED.*.bgl.phased.phased >& /dev/null
    rm $OUTPUT.bgl.log >& /dev/null
    rm -r $JAVATMP >& /dev/null
    rm -f plink.log >& /dev/null
    echo "DONE!"
    echo ""
endif

