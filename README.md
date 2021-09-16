# SNP2HLA_UM
**SNP2HLA_UM** can construct HLA reference panels and impute HLA genotypes.


## Introduction

SNP2HLA_UM is an updated version of the original [v1.0.3 SNP2HLA package](http://software.broadinstitute.org/mpg/snp2hla/).

The most important changes include:

1. using Beagle v5.2 instead of v3.0.4 for faster and more accurate phasing and imputation;

1. using PLINK 1.90 instead of PLINK v1.07 for increased speed;

1. updated HLA SNP and amino acid dictionaries based on a newer release (3.30) of the [IPD-IMGT/HLA
    database](https://www.ebi.ac.uk/ipd/imgt/hla/) and now including HLA-DRB3/4/5 in addition to the eight classical HLA
    genes of the original SNP2HLA package.

## Contents of package

1. **MakeReference_UM** is a collection of scripts and data files to build a reference panel for HLA imputation;

1. **SNP2HLA_UM** is a collection of scripts to impute HLA SNP, amino acid and classical allele genotypes;

1. **Reference panels:** three pre-built HLA reference panels in hg19 coordinates based on 233 individuals of South Asian ancestry
(UM-SAS_REF), 2504 individuals from phase 3 of the 1000 Genomes Project (1KGP_REF), and a panel that merges these
two panels (UM-SAS+1KGP_REF). The HLA genes in these panels were typed at true 2-field (4-digit) resolution, 
unlike the majority of HLA panels where typing was performed at only G-group resolution.

## Requirements

SNP2HLA_UM runs in a Linux environment; it can run on a Windows machine using the Windows Subsystem for Linux.

Linux packages
- csh
- openjdk (version 8.0 or later)
- libgfortran4 (runtime library for GNU fortran apps)

Other software
- [PLINK v1.9](https://www.cog-genomics.org/plink/1.9/)
- [Beagle v5.2](https://faculty.washington.edu/browning/beagle/beagle.html)

Genetic maps (optional)
- [plink.chr6.GRCh37.map](http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/)
- [plink.chr6.GRCh38.map](http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/)

## Installation

1. Install any of the required Linux packages that aren't already on your system: on Debian-based distributions `sudo apt update` followed by
`sudo apt install csh`, `sudo apt install openjdk`, and/or `sudo apt install libgfortran4`.

1. If you encounter difficulties running the provided Fortran executables (hack_vcf, purge_SNPs, revert_alleles, tped2vcf) with the
GNU Fortran runtime libraries, an alternative is to install the GNU Fortran compiler <`sudo apt install gfortran`> and then compile the
Fortran source files in folder /MakeReference_UM/fortran_code or /SNP2HLA_UM/fortran_code.

1. Download PLINK v1.9 and Beagle v5.2 and place both executable files into the directory where you will run SNP2HLA_UM. 
Rename the beagle file as beagle.jar.

1. Optionally, download the genetic map file for the genome build you are working with and place it into the same directory. 
Beagle recommends using these genetic map files.

1. **MakeReference_UM;** copy from /MakeReference_UM/dictionaries and /MakeReference_UM/scripts into the working directory:

    - HLA_DICTIONARY_hg19.tar.gz or HLA_DICTIONARY_hg38.tar.gz
    - MakeReference_UM.csh
    - encodeHLA_UM.pl
    - encodeVariants.pl
    - HLAtoSequences_UM.pl
    - purge_SNPs
    - revert_alleles
    - tped2vcf

1. **SNP2HLA_UM;** copy from /SNP2HLA_UM/scripts into the working directory:

    - SNP2HLA_UM.csh
    - hack_vcf
    - merge_tables.pl
    - revert_alleles

1. Uncompress the tar.gz file for the data dictionaries if running MakeReference_UM

1. If necessary, change to executable the permissions of the .csh files, .pl files, and programs purge_SNPs, revert_alleles,
tped2vcf and hack_vcf.

## Input and output files

Examples of input and output files can be found in directories /MakeReference_UM/example_input, /MakeReference_UM/example_output,
/SNP2HLA_UM/example_input and /SNP2HLA_UM/example_output.

#### MakeReference_UM input

- Plink .bed/bim/fam files with genotype data for SNPs in the MHC region (chr6:29-34 Mb)
- HLA pedigree data file holding sample information and 2-field HLA genotypes for 11 HLA genes; columns are in the following order:
FID, IID, pID, mID, SEX, PHENO, HLA-A, HLA-B, HLA-C, HLA-DPA1, HLA-DPB1, HLA-DQA1, HLA-DRB1, HLA-DRB3, HLA-DRB4, HLA-DRB5.  Untyped
alleles are denoted with a 0; missing copies of HLA-DRB3/4/5 alleles are denoted with a 99:99.  All columns must be present, even
for HLA genes with no typing for the entire dataset.

#### MakeReference_UM output

- .bgl.vcf file holding phased data for the newly created reference panel
- Plink .bed/bim/fam files holding unphased genotypes for the newly created reference panel
- Plink .FRQ.frq file with minor allele frequencies for each variant in the reference panel
- .markers file holding variant ID, bp position and alleles for each variant in the reference panel

#### SNP2HLA_UM input

- six files for the HLA reference panel created by MakeReference_UM
- Plink ./bed/bim/fam files holding genotype data in the MHC region for the target dataset to be
imputed; SNP2HLA_UM matches variant IDs between the target dataset and reference panel, so it is
crucial that variants are named consistently in both. 

#### SNP2HLA_UM output

- .vcf file holding phased and imputed HLA and MHC genotypes for the target dataset.

## Test examples

### MakeReference_UM

```
$ ./MakeReference_UM.csh \
      HAPMAP_CEU_hg19 \
      HAPMAP_CEU_HLA_new.ped \
      HAPMAP_CEU_hg19_REF \
      plink \
      hg19 \
      5 \
      15 \
      4 \
      10 \
      20 \
`#    plink.chr6.GRCh37.map` \
   > MakeReference_hg19.log
```
Notes:
- To run this example copy the input files from directory /MakeReference_UM/example_input/hg19 into your working directory
that has the data dictionaries and programs and scripts required to run MakeReference_UM.
- The first 5 command line arguments are mandatory and specify, respectively, the basename of the Plink .bed/bim/fam files holding
MHC SNP gentoypes, the name of the HLA genotype pedigree file, the basename for the six output reference panel files, the name of
the PLINK v1.9 executable, and the human reference build (hg19 or hg38) being used.
- The last 6 command line arguments are optional; they specify, respectively, the amount of Java heap memory in Gb, the maximum
Java stacksize in Mb, the number of CPU threads to be used by Beagle, the number of Beagle burnin interations, the number
of Beagle phasing iterations, and the name of the genetic map file.  Defaults for these arguments are 2 Gb, 5 Mb, 
1 thread, 6 burnin iterations, 12 phasing iterations, and no genetic map file.
- The command line argument for the genetic map file is commented out to avoid requiring downloading it for this test example.
- If everything is working properly, running this example should produce seven output files with content closely matching those in directory 
/MakeReference_UM/example_output/hg19.

### SNP2HLA_UM

```
$ ./SNP2HLA_UM.csh \
      1958BC_hg19 \
      HAPMAP_CEU_hg19_REF \
      1958BC_hg19_imputed \
      plink \
      5 \
      15 \
      4 \
      10 \
      20 \
`#    plink.chr6.GRCh37.map` \
  > SNP2HLA_hg19.log
```
Notes:
- To run this example copy the input files from directory /SNP2HLA_UM/example_input/hg19 into your working directory
that contains the programs and scripts required to run SNP2HLA_UM.
- The first 4 command line arguments are mandatory and specify, respectively, the basename of the Plink .bed/bim/fam files holding
MHC SNP gentoypes for the target dataset, the basename of the six reference panel files (./bed/bim/fam/bgl.vcf/FRQ.frq/markers)
produced by MakeReference_UM, and the name of the PLINK v1.9 executable.
- The last 6 command line arguments are optional; they specify, respectively, the amount of Java heap memory in Gb, the maximum
Java stacksize in Mb, the number of CPU threads to be used by Beagle, the number of Beagle burnin interations, the number
of Beagle phasing iterations, and the name of the genetic map file.  Defaults for these arguments are 2 Gb, 5 Mb, 
1 thread, 6 burnin iterations, 12 phasing iterations, and no genetic map file.
- The command line argument for the genetic map file is commented out to avoid requiring downloading it for this test example.
- If everything is working properly, running this example should produce two output files with content closely matching those in directory 
/SNP2HLA_UM/example_output/hg19.


## Support
Before contacting us, please read all of the instructions carefully and try out the test usage examples.  If you still have
issues or questions, contact:

Philip Stuart \
pstuart@umich.edu \
Department of Dermatology \
University of Michigan \
Ann Arbor, Michigan 48109, USA

## References

If you use SNP2HLA_UM in your work please cite *both* of the following studies:

- Jia X, Han B, Onengut-Gumuscu S, Chen WM, Concannon PJ, Rich SS, Raychaudhuri S, de Bakker PI (2013) Imputing amino acid polymorphisms in 
human leukocyte antigens. PLoS One. Jun 6;8(6):e64683.

- Stuart PE, Tsoi LC, Nair RP, Ghosh M, Kabra M, Shaiq PA, Raja GK, Qamar R, Thelma BK, Patrick MT, Parihar A, Singh S, Khandpur S,
Kumar U, Wittig M, Degenhardt F, Tejasvi T, Voorhees JJ, Weidinger S, Franke A, Abecasis GR, Sharma VK, Elder JT (2021) Transethnic
analysis of psoriasis susceptibility in South Asians and Europeans enhances fine-mapping in the MHC and genomewide. HGG Advances,
under review.

If you use any of the pre-built reference panels, please cite the Stuart et al. study.

If you use either of the pre-built reference panels with phase 3 1000 Genomes data, please also cite these two studies:

- Abi-Rached L, Gouret P, Yeh JH, Di Cristofaro J, Pontarotti P, Picard C, Paganini J (2018) Immune diversity sheds light on missing
variation in worldwide genetic diversity panels. PLoS One. Oct 26;13(10):e0206512.

- 1000 Genomes Project Consortium, Auton A, Brooks LD, Durbin RM, Garrison EP, Kang HM, Korbel JO, Marchini JL, McCarthy S, McVean GA, 
Abecasis GR (2015) Nature Oct 1;526(7571):68-74.

If you find SNP2HLA_UM useful, you may also be interested in [**HLA-TAPAS**](https://github.com/immunogenomics/HLA-TAPAS), 
another recent updating of the original SNP2HLA package.

