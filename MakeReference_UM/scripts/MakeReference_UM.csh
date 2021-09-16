#!/bin/csh -f

########################################################################################################################################################################################################################
#   MakeReference_UM.csh
#
#   DESCRIPTION: This script helps prepare a reference dataset for HLA imputation
#
#   Author : Sherman Jia(xiaomingjia@gmail.com)
#           + Modifications by Phil Stuart (pstuart@umich.edu):
#             (a) Adapted script to use Beagle 5.2 instead of Beagle 3.0.4 as its imputation
#                 engine, thereby providing substantially more accurate and faster imputation
#                 and phasing. Verified to work with release 05Apr21.9b7 of Beagle 5.2 and
#                 also with the last version (18May20.d20) of Beagle 5.1.
#             (b) Uses newly constructed HLA data dictionaries that were updated to a more
#                 recent release (3.30) of the IPD-IMGT/HLA database.  These dictionaries now
#                 include HLA-DRB3, -DRB4 and -DRB5 in addition to the eight classical HLA 
#                 genes (A,B,C,DPA1,DPB1,DQA1,DQB1,DRB1), are based on modern post-2010 HLA
#                 nomenclature, and use ambiguous allele codes to accurately denote bases in
#                 the HLA SNP dictionary for which not all 3-field subtypes of a 2-field
#                 protein allele have the same nucleotide allele.
#             (c) Script now uses Plink 1.9 instead of v1.07 because of its improved speed.
#             (d) Script now has flexibility to handle either of the two most recent builds
#                 of the human genome reference sequence (hg19 & hg38).
#             (e) Script now outputs the phased reference panel in .vcf format rather than
#                 the less standard .bgl.phased format used by Beagle 3.0.4.  Unphased
#                 reference panel data are still output in Plink .bed/bim/fam format.
#             (f) Created parameter MAF_MIN to facilitate changing the minimum allowed MAF
#                 for reference panel SNPs from the default of 0.01.
#             (g) Values for heap memory, stack memory, number of threads and the number of
#                 burn-in and phasing iterations used by Beagle can optionally be specified 
#                 on the command line when launching this script.
#             (h) The command line can now optionally specify the name of a PLINK format
#                 genetic map file with cM units, which can be used by Beagle to enhance the
#                 accuracy of its phasing and imputation.
#                  
#   IMPORTANT NOTES:
#
#     (1) The input Plink .bed/.bim/.fam files contain the SNP data.
#
#     (2) The input HLA .ped file contains sample information along with the 4-digit HLA alleles; 
#         columns are (FID,IID,pID,mID,SEX,PHENO,A,B,C,DPA1,DPB1,DQA1,DQB1,DRB1,DRB3,DRB4,DRB5).  
#         Use zeros as placeholders for missing HLA genotypes; fields must be present for all 11 HLA genes.
#
#     (3) HLA alleles in HLA PED file must use the newer xx:xx nomenclature rather than the older xxxx 
#         nomenclature required by older versions of the MakeReference script (e.g., 06:02 instead 
#         of 0602) in order to accomodate 3-digit fields.  Absent HLA-DRB3, DRB4 or DRB5 alleles 
#         due to the the lack of the gene on one or both chromosomes for an individual should be denoted
#         as allele 99:99. 
#
#     (4) When constructing multiethnic reference panels, the user may want to modify or comment out the
#         following shell script command that removes SNPs violating Hardy-Weinberg equilibrium:
#
#         awk '{if (NR > 1){print}}' $SNP_DATA.FOUNDERS.hardy.hwe | awk ' $9 < 0.000001 { print $2 }' | sort - u > remove.snps.hardy
#
#         This command is meant as a QC step to detect poorly genotyped variants, since these will often
#         have genotype frequencies that deviate strongly from the expected HWE.However, HWE is an
#         appropriate expectation for a single population in genetic equilibrium, but not for a
#         combination of populations where the genotype frequecies for a variant may vary greatly
#         among the populations.
#
#   DEPENDENCIES: (download and place in the same folder as this script)
#     1. PLINK (1.9); available at https://www.cog-genomics.org/plink/1.9/
#     2. Beagle (5.2); available at https://faculty.washington.edu/browning/beagle/beagle.html
#     3. tped2vcf (utility by Phil Stuart for Plink .tped -> .vcf format)
#     4. HLAtoSequences_UM.pl (Perl script to generate amino acid sequences from HLA types; modified by P. Stuart from original HLAtoSequences.pl)
#     5. encodeVariants.pl (Perl script to encode multi-allelic PLINK markers (AAs and SNPs) into bi-allelic markers)
#     6. encodeHLA_UM.pl (Perl script to encode HLA alleles into biallelic markers; modified by P. Stuart from original encodeHLA.pl)
#     7. HLA dictionary files (HLA_DICTIONARY_AA.hg*.imgt3300.map/txt; HLA_DICTIONARY_SNPS.hg*.imgt3300.map/txt; newly constructed by P. Stuart and provided in both hg19 and hg38 coordinates)
#     8. revert_alleles (Utility by Phil Stuart to revert hacked non-ACTG alleles in *.vcf output file, necessitated by Beagle 5.2 allele coding restrictions, back to allele identities in reference panel)
#     9. purge_SNPs (Utility by Phil Stuart to delete SNPs with ambiguous or missing alleles)
#    10. If genetic_map_file argument specified, PLINK format genetic map on cM scale with build 37 or 38 coordinates (plink.chr6.GRCh37.map or plink.chr6.GRCh38.map), downloaded from http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/
#
#   USAGE: ./MakeReference_UM.csh SNPS[.bed/.bim/.fam] HLA.ped OUTPUT plink hg_version max_memory[gb] max_stacksize[mb] nthreads nburnin niterations genetic_map_file
########################################################################################################################################################################################################################

if ($#argv < 5) then
    echo "USAGE: ./MakeReference_UM.csh SNPS (.bed/.bim/.fam) HLAfile (FID,IID,pID,mID,SEX,PHENO,A,B,C,DPA1,DPB1,DQA1,DQB1,DRB1,DRB3,DRB4,DRB5) OUTPUT (.bgl.vcf/.markers/.bed) plink hg_version_human_reference {optional: java_max_memory[gb] java_max_stacksize[mb] number_of_threads number_of_burnin-iterations number_of_phasing-iterations plink_format_genetic_map_file}"; exit 1
endif

set SCRIPTPATH=`dirname $0`

# CHECK FOR DEPENDENCIES
if (! -e $SCRIPTPATH/$4) then
    echo "Please copy the binary for v1.9 of Plink (https://www.cog-genomics.org/plink/1.9/) into this directory."; exit 1
else if (! -e $SCRIPTPATH/beagle.jar) then
    echo "Please download the binary for Beagle 5.2 (https://faculty.washington.edu/browning/beagle/beagle.html), rename it as beagle.jar, and copy it into this directory."; exit 1
else if (! -e $SCRIPTPATH/tped2vcf) then
    echo "Please copy tped2vcf (included with this package) into $SCRIPTPATH/"; exit 1
else if (! -e $SCRIPTPATH/HLAtoSequences_UM.pl) then
    echo "Please copy HLAtoSequences_UM.pl (included with this package) into this directory."; exit 1
else if (! -e $SCRIPTPATH/encodeVariants.pl) then
    echo "Please copy encodeVariants.pl (included with this package) into this directory."; exit 1
else if (! -e $SCRIPTPATH/encodeHLA_UM.pl) then
    echo "Please copy encodeHLA_UM.pl (included with this package) into this directory."; exit 1
else if ($5 == "hg19" && ! -e $SCRIPTPATH/HLA_DICTIONARY_AA.hg19.imgt3300.map) then
    echo "Please copy HLA_DICTIONARY_AA.hg19.imgt3300.map (included with this package) into this directory."; exit 1
else if ($5 == "hg19" && ! -e $SCRIPTPATH/HLA_DICTIONARY_AA.hg19.imgt3300.txt) then
    echo "Please copy HLA_DICTIONARY_AA.hg19.imgt3300.txt (included with this package) into this directory."; exit 1
else if ($5 == "hg19" && ! -e $SCRIPTPATH/HLA_DICTIONARY_SNPS.hg19.imgt3300.map) then
    echo "Please copy HLA_DICTIONARY_SNPS.hg19.imgt3300.map (included with this package) into this directory."; exit 1
else if ($5 == "hg19" && ! -e $SCRIPTPATH/HLA_DICTIONARY_SNPS.hg19.imgt3300.txt) then
    echo "Please copy HLA_DICTIONARY_SNPS.hg19.imgt3300.txt (included with this package) into this directory."; exit 1
else if ($5 == "hg38" && ! -e $SCRIPTPATH/HLA_DICTIONARY_AA.hg38.imgt3300.map) then
    echo "Please copy HLA_DICTIONARY_AA.hg38.imgt3300.map (included with this package) into this directory."; exit 1
else if ($5 == "hg38" && ! -e $SCRIPTPATH/HLA_DICTIONARY_AA.hg38.imgt3300.txt) then
    echo "Please copy HLA_DICTIONARY_AA.hg38.imgt3300.txt (included with this package) into this directory."; exit 1
else if ($5 == "hg38" && ! -e $SCRIPTPATH/HLA_DICTIONARY_SNPS.hg38.imgt3300.map) then
    echo "Please copy HLA_DICTIONARY_SNPS.hg38.imgt3300.map (included with this package) into this directory."; exit 1
else if ($5 == "hg38" && ! -e $SCRIPTPATH/HLA_DICTIONARY_SNPS.hg38.imgt3300.txt) then
    echo "Please copy HLA_DICTIONARY_SNPS.hg38.imgt3300.txt (included with this package) into this directory."; exit 1
else if (! -e $SCRIPTPATH/revert_alleles) then
    echo "Please copy revert_alleles (included with this package) into this directory."; exit 1
else if (! -e $SCRIPTPATH/purge_SNPs) then
    echo "Please copy purge_SNPs (included with this package) into into this directory."; exit 1
else if ($#argv >= 11 && ! -e $SCRIPTPATH/$11) then
    echo "Please download appropriate Plink-format genetic recombination map from http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/ and copy it into into this directory."; exit 1
endif

set SNP_DATA=$1
set HLA_DATA=$2
set OUTPUT=$3
set PLINK=$4
set HG=$5

if ($#argv >= 6) then
    set MEM=$6
else
    set MEM=2 # Default java memory 2 Gb
endif

if ($#argv >= 7) then
    set SS=$7
	set PLINK_MEM="$6"000""
else
    set SS=5 # Default java stacksize 5 Mb
	set PLINK_MEM=5000
endif

if ($#argv >= 8) then
    set THREAD=$8
else
    set THREAD=1
endif

if ($#argv >= 9) then
    set BURNIN=$9
else
    set BURNIN=6
endif

if ($#argv >= 10) then
    set ITER=$10
else
    set ITER=12
endif

if ($#argv >= 11) then
    set MAP=$11
endif

#CHECK FOR PRESENCE OF SNP AND HLA INPUT FILES IN DIRECTORY

if (! -f $SNP_DATA.bed) then
  echo "Input file $SNP_DATA.bed not found"; exit 1
else if (! -f $SNP_DATA.bim) then
  echo "Input file $SNP_DATA.bim not found"; exit 1
else if (! -f $SNP_DATA.fam) then
  echo "Input file $SNP_DATA.fam not found"; exit 1
else if (! -f $HLA_DATA) then
  echo "Input file $HLA_DATA not found"; exit 1
endif

alias plink './$PLINK --silent --allow-no-sex --memory $PLINK_MEM --threads $THREAD'
alias beagle 'java -Xss$SS\m -Xmx$MEM\g -jar $SCRIPTPATH/beagle.jar'
alias tped2vcf './tped2vcf -chr 6 -missing 0 -recode_opt 1 -revert_opt 1'
alias revert_alleles './revert_alleles -script_opt 2'
alias purge_SNPs './purge_SNPs'

set MAF_MIN = 0.01

set ENCODE_AA          = 1
set ENCODE_HLA         = 1
set ENCODE_SNPS        = 1
set EXTRACT_FOUNDERS   = 1
set MERGE              = 1
set QC                 = 1
set PREPARE            = 1
set PHASE              = 1
set CLEANUP            = 1

set i=1

echo "MakeReference_UM: Creating reference panel $OUTPUT"
echo "- Java memory = "$MEM"Gb"
echo "- Java stacksize = "$SS"Mb"
echo "- Beagle number of threads = "$THREAD" threads"
echo "- Beagle number of burnin interations = $BURNIN interations"
echo "- Beagle number of phasing iterations = "$ITER" iterations"
if ($#argv >= 11) then
    echo "- Beagle genetic map file = "$MAP""
else
    echo "- Beagle genetic map file = none specified"
endif
echo ""

# Encode HLA amino acids
if ($ENCODE_AA) then
    echo "[$i] Generating amino acid sequences from HLA types.";  @ i++
    ./HLAtoSequences_UM.pl $HLA_DATA HLA_DICTIONARY_AA.$HG.imgt3300.txt AA > $OUTPUT.AA.ped
    cp HLA_DICTIONARY_AA.$HG.imgt3300.map $OUTPUT.AA.map

    echo "[$i] Encoding amino acids positions." ;  @ i++
    ./encodeVariants.pl $OUTPUT.AA.ped $OUTPUT.AA.map $OUTPUT.AA.CODED

    plink --file $OUTPUT.AA.CODED --missing-genotype 0 --make-bed --out $OUTPUT.AA.TMP
    awk '{if ($5 == "0" || $5 == "x" || $6 == "x"){print $2}}' $OUTPUT.AA.TMP.bim | grep -v INS | cut -f2 > to_remove
    plink --bfile $OUTPUT.AA.TMP --exclude to_remove --make-bed --out $OUTPUT.AA.CODED
    rm $OUTPUT.AA.TMP*; rm to_remove
    rm $OUTPUT.AA.???
endif

# Encode classical HLA alleles into binary format
if ($ENCODE_HLA) then
    echo "[$i] Encoding HLA alleles.";  @ i++
    ./encodeHLA_UM.pl $HLA_DATA $OUTPUT.HLA.map $HG > $OUTPUT.HLA.ped
    plink --file $OUTPUT.HLA --make-bed --out $OUTPUT.HLA
endif

# Encode HLA SNPs
if ($ENCODE_SNPS) then
    echo "[$i] Generating DNA sequences from HLA types.";  @ i++
    ./HLAtoSequences_UM.pl $HLA_DATA HLA_DICTIONARY_SNPS.$HG.imgt3300.txt SNPS > $OUTPUT.SNPS.ped
    cp HLA_DICTIONARY_SNPS.$HG.imgt3300.map $OUTPUT.SNPS.map

    echo "[$i] Encoding SNP positions." ;  @ i++
    ./encodeVariants.pl $OUTPUT.SNPS.ped $OUTPUT.SNPS.map $OUTPUT.SNPS.CODED
    plink --file $OUTPUT.SNPS.CODED --missing-genotype 0 --make-bed --out $OUTPUT.SNPS.TMP

    awk '{if ($5 == "0" || $5 == "x" || $6 == "x"){print $2}}' $OUTPUT.SNPS.TMP.bim | grep -v INS | cut -f2 > to_remove
    plink --bfile $OUTPUT.SNPS.TMP --exclude to_remove --make-bed --out $OUTPUT.SNPS.CODED
	purge_SNPs -fnbim $OUTPUT.SNPS.CODED.bim -fnout to_remove2
	plink --bfile $OUTPUT.SNPS.CODED --exclude to_remove2 --make-bed --out $OUTPUT.SNPS.CODED >& /dev/null
    rm $OUTPUT.SNPS.TMP*; rm to_remove; rm to_remove2
    rm $OUTPUT.SNPS.???
endif

if ($EXTRACT_FOUNDERS) then
    echo "[$i] Extracting founders."; @ i++
    plink --bfile $SNP_DATA --filter-founders --mind 0.3 --alleleACGT --make-bed --out $SNP_DATA.FOUNDERS

    # Initial QC on Reference SNP panel
    plink --bfile $SNP_DATA.FOUNDERS --hardy --out $SNP_DATA.FOUNDERS.hardy
    plink --bfile $SNP_DATA.FOUNDERS --freq --out $SNP_DATA.FOUNDERS.freq
    plink --bfile $SNP_DATA.FOUNDERS --missing --out $SNP_DATA.FOUNDERS.missing
    awk '{if (NR > 1){print}}' $SNP_DATA.FOUNDERS.hardy.hwe | awk ' $9 < 0.000001 { print $2 }' | sort -u > remove.snps.hardy 
    awk '{if (NR > 1){print}}' $SNP_DATA.FOUNDERS.freq.frq | awk ' $5 < '$MAF_MIN' { print $2 } ' > remove.snps.freq
    awk '{if (NR > 1){print}}' $SNP_DATA.FOUNDERS.missing.lmiss | awk ' $5 > 0.05 { print $2 } ' > remove.snps.missing
    cat remove.snps.* | sort -u > all.remove.snps

    plink --bfile $SNP_DATA.FOUNDERS --exclude all.remove.snps --make-bed --out $SNP_DATA.FOUNDERS.QC

    # Founders are identified here as individuals with "0"s in mother and father IDs in .fam file

    plink --bfile $OUTPUT.HLA --filter-founders --maf 0.00005 --make-bed --out $OUTPUT.HLA.FOUNDERS
    plink --bfile $OUTPUT.SNPS.CODED --filter-founders --maf 0.00005 --make-bed --out $OUTPUT.SNPS.FOUNDERS
    plink --bfile $OUTPUT.AA.CODED --filter-founders --maf 0.00005 --make-bed --out $OUTPUT.AA.FOUNDERS

	rm remove.snps.* >& /dev/null
endif

# Merging SNP, HLA, and amino acid datasets 
if ($MERGE) then
    echo "[$i] Merging SNP, HLA, and amino acid datasets.";  @ i++
    echo "$OUTPUT.HLA.FOUNDERS.bed $OUTPUT.HLA.FOUNDERS.bim $OUTPUT.HLA.FOUNDERS.fam" > merge_list
    echo "$OUTPUT.AA.FOUNDERS.bed $OUTPUT.AA.FOUNDERS.bim $OUTPUT.AA.FOUNDERS.fam" >> merge_list
    echo "$OUTPUT.SNPS.FOUNDERS.bed $OUTPUT.SNPS.FOUNDERS.bim $OUTPUT.SNPS.FOUNDERS.fam" >> merge_list
    plink --bfile $SNP_DATA.FOUNDERS.QC --merge-list merge_list --make-bed --out $OUTPUT.MERGED.FOUNDERS >& /dev/null
    # plink --bfile $SNP_DATA.FOUNDERS.QC --merge-list merge_list --make-bed --out $OUTPUT.MERGED.FOUNDERS
    rm $OUTPUT.HLA.??? >& /dev/null
    rm $OUTPUT.HLA.????? >& /dev/null
    rm $OUTPUT.AA.CODED.??? >& /dev/null
	rm $OUTPUT.AA.CODED.????? >& /dev/null
    rm $OUTPUT.SNPS.CODED.??? >& /dev/null
    rm $OUTPUT.SNPS.CODED.???? >& /dev/null
	rm $OUTPUT.SNPS.CODED.????? >& /dev/null
    rm merge_list >& /dev/null
endif

if ($QC) then
    echo "[$i] Performing quality control.";  @ i++
    plink --bfile $OUTPUT.MERGED.FOUNDERS --freq --out $OUTPUT.MERGED.FOUNDERS.FRQ
    awk '{if (NR > 1 && ($5 < 0.0001 || $5 > 0.9999)){print $2}}' $OUTPUT.MERGED.FOUNDERS.FRQ.frq > all.remove.snps
    awk '{if (NR > 1){if (($3 == "A" && $4 == "P") || ($4 == "A" && $3 == "P")){print $2 "\tP"}}}' $OUTPUT.MERGED.FOUNDERS.FRQ.frq > allele.order

    # QC: Maximum per-SNP missing > 0.5, MAF > 0.1%
    plink --bfile $OUTPUT.MERGED.FOUNDERS --reference-allele allele.order --exclude all.remove.snps --geno 0.5 --make-bed --out $OUTPUT

    # Calculate allele frequencies
    plink --bfile $OUTPUT --keep-allele-order --freq --out $OUTPUT.FRQ
    rm $SNP_DATA.FOUNDERS.* >& /dev/null
    rm $OUTPUT.MERGED.FOUNDERS.* >& /dev/null
    rm $OUTPUT.*.FOUNDERS.??? >& /dev/null
	rm $OUTPUT.*.FOUNDERS.????? >& /dev/null
    rm allele.order >& /dev/null
    rm all.remove.snps >& /dev/null
endif

if ($PREPARE) then
    # Prepare files for Beagle phasing
    echo "[$i] Preparing files for Beagle. ";  @ i++
    awk '{print $2 " " $4 " " $5 " " $6}' $OUTPUT.bim > $OUTPUT.markers
    plink --bfile $OUTPUT --keep-allele-order --recode --alleleACGT --out $OUTPUT
	plink --bfile $OUTPUT --recode transpose --out $OUTPUT

    echo "[$i] Converting to vcf format.";  @ i++
    tped2vcf -fnmarker $OUTPUT.markers -fnfam $OUTPUT.fam -fngtype $OUTPUT.tped -fnvcf $OUTPUT.vcf
endif 

if ($PHASE) then
    echo "[$i] Phasing reference using Beagle.";  @ i++
	if ($#argv >= 11) then
        beagle gt=$OUTPUT.vcf nthreads=$THREAD chrom=6 burnin=$BURNIN iterations=$ITER out=$OUTPUT.bgl map=$MAP
	else
        beagle gt=$OUTPUT.vcf nthreads=$THREAD chrom=6 burnin=$BURNIN iterations=$ITER out=$OUTPUT.bgl
	endif

    gunzip -d $OUTPUT.bgl.vcf.gz
    echo ""
	echo "[$i] Reverting hacked allele identities in *.bgl.vcf file back to those traditionally used by SNP2HLA."; @ i++
    revert_alleles -fnvcf $OUTPUT.bgl.vcf -fnref $OUTPUT.bim
endif

if ($CLEANUP) then
    rm $OUTPUT.ped >& /dev/null
    rm $OUTPUT.map >& /dev/null
    rm $OUTPUT.tped >& /dev/null
    rm $OUTPUT.tfam >& /dev/null
	rm $OUTPUT.vcf >& /dev/null
    rm $OUTPUT.FRQ.log >& /dev/null
    rm $OUTPUT.log >& /dev/null
    rm $OUTPUT.bgl.log >& /dev/null
endif

echo "[$i] Done."
