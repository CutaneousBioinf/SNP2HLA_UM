      program tped2vcf
!--------------------------------------------------------------------------------------
!     Program tped2vcf is a utility program to be used as part of my modification
!     of the SNP2HLA package of Jia et al (2013) that allows it to use Beagle v5.2
!     rather than Beagle 3.0.4.  This utility converts a Plink pedigree file with 
!     unphased genotypes (*.tped) into v4.2 vcf format, which is needed for my 
!     modification of the MakeReference C shell script.    
!
!     The unphased genotypes input file (Plink *.tped format) is generated 
!     internally by the MakeReference script using Plink with a --transpose flag.
!     The tped format has one row per variant, with the first four values in the 
!     row being chromosome number, marker ID, recombination distance (usually 0) 
!     and basepair position, followed by genotypes for that variant for all
!     samples.  Marker IDs may be up to 100 characters in length.  Genotypes must
!     be listed as two alleles, separated by at least one space and/or tab; each 
!     allele is a character string of up to 100 characters in length that may be 
!     composed of any combination of characters in the standard ASCII set.
!
!     In addition to the unphased genotypes .tped file, this program also requires
!     as input a beagle v3 markers input file.  This file has no header lines. 
!     Each row contains four fields--marker ID, basepair postion, allele 1 
!     identity, and allele 2 indentity.  The marker IDs may be up to 100 characters
!     in length; allele IDs may be up tp 100 characters in length and may use any 
!     character in the standard ASCII set.  All markers must be biallelic; the 
!     user is asked to specify which character designates a missing allele in
!     the unphased genotype input file.
!
!     Other input parameters supplied by the user include the chromosome ID for 
!     the input markers, and whether or not markers with non-ACTG alleles and/or 
!     so-called duplicated markers should be recoded into a format acceptable to
!     beagle 5.2.  If the user requests this last option, the program first looks 
!     for all sets of markers that are identical for both basepair position and 
!     ref/alt alleles (i.e., amino acid and 2-digit or 4-digit HLA alleles with 
!     same position and P/A alleles, or a microarray HLA SNP and its MakeReference-
!     computed SNP-dictionary equivalent, or very rarely an AA and HLA SNP with 
!     same position and happenchance matching of 1 letter AA and nucleotide codes) 
!     and artificially distinguishes them as distinct markers by appropriately 
!     changing the ref/alt alleles.  Duplicate sets with up to 4 markers are 
!     converted to all possible unique combinations of ref=A and alt=AX alleles 
!     where X is all 4 permutations of A/C/T/G, sets with 5-16 markers are 
!     converted to ref=A and alt=AXX alleles where XX is all 16 permutations of 
!     A/C/T/G, sets with 17-64 markers are converted to ref=A and alt=AXXX alleles 
!     where XXX is all 64 permutations of A/C/T/G, and sets with 65-256 markers 
!     are converted to ref=A and alt=AXXXX where XXXXX is all 256 permutations 
!     of A/C/T/G.  After duplicated markers are made distinguishable, the program 
!     then checks if any markers have alleles that aren't a combination of A/C/T/G 
!     or the missing character code 0 (usually amino acid positions with only 
!     two possible residues).  These, if any, are converted to A/AA alleles.  
!     Because this conversion can create new sets of apparently identical markers, 
!     particularly for HLA-DRB3/4/5 AA and SNP variants with identical positions 
!     and one Z allele, another round of deduplication is necessary. 
!
!     The deduplication scheme outlined above can cause an inadvertent mismatch
!     between variants that are in both the target and the reference panel.  If the 
!     reference panel contains a microarray-genotyped SNP in an HLA gene that is 
!     also present as a computed SNP variant based on IMGT-SNP dictionary sequences, 
!     the hacked alleles assigned to the microarray-genotyped version of the SNP 
!     in the reference panel will not match the true alleles for this SNP in the 
!     target dataset.  To remedy this situation, after deduplication all markers 
!     that are not microarray-genotyped variants (i.e., all markers whose IDs don't 
!     begin with 'AA_', 'SNP_', 'HLA_', 'INS_') are reassigned their original allele 
!     identities before deduplication.  Because this procedure could possibly 
!     recreate apparently duplicate markers (i.e., same bp position and alleles;
!     although unlikely, this could happen when an indel has a single base insertion
!     with a padded A preceding the insertion point), the option to impose the 
!     procedure is user-selectable.
!     
!     Output is a file in v4.2 vcf format.  The header lines at the top of the
!     file list basic information, including the vcf format version, the date of
!     file creation, the source of file creation, the contig ID (chrom no.), the
!     length in bp spanned by the input markers, descriptions of the INFO and FORMAT
!     fields, and a header line that includes all sample IDs.  All genotypes in
!     the output file are converted to 0/1 format, where 0 is A1 and 1 is A2.  
!     A "GT" is listed in the format field as the output is genotypes.  Alleles for 
!     each genotype are separated a "\" to denote they are unphased.  Missing 
!     genotypes are denoted as "./." per vcf format specifications.
!
!     Written by Philip Stuart in standard Fortran 90.                 
!                                                                      
!     Revision History:												 
!                                                                      
!     v 1.0  04/08/21   First version.  Represents a reduction of the old
!                       beagle2vcf utility (v1.7), which was originally used to
!                       convert either an old-style .bgl.phased Beagle v3.0.4
!                       output file or a Plink .tped froamt file into a vcf file.
!
!     References:
!
!     Jia X, Han B, Onengut-Gumuscu S, Chen W-M, Concannon PJ, Rich SS, Raychaudhuri S, 
!        de Bakker PIW (2013) Imputing amino acid polymorphisms in human leukocyte
!        antigens. PLoS One 8(6) e64683.
!
!     SNP2HLA software, version 1.0.3, downloaded in May 2016 from 
!        https://www.broadinstitute.org/impg/snp2hla
!--------------------------------------------------------------------------------------
      implicit none

      integer i,item,ioerr,imk,istart,iend,j,k, &
              nskip,nmarker,nsample,ndup, &
              chr,recode_opt,revert_opt,gtype_opt,i_position,pos_first,pos_last
      integer,allocatable :: position(:)
      character(1) accept,missing,delimiter
      character(2) aset1(4,2)
      character(3) aset2(16,2)
      character(4) aset3(64,2)
      character(5) aset4(256,2)
      character(6) aset5(1024,2)
      character(6) string
      character(8) date
      character(35) answer(2)
      character(45) ctmp
      character(200) mkid,mkid1,mkid2,a1,a2
      character(200), allocatable :: samples(:), ralleles(:,:), ralleles_original(:,:),alleles(:,:)
      character(1000) fnmarker,fnfam,fngtype,fnvcf
      character(1001) arg
      logical exist
      
      data aset1(:,1) /4*'A '/
      data aset1(:,2) /'AA','AC','AG','AT'/
      
      data aset2(:,1) /16*'A  '/
      data aset2(:,2) /'AAA','AAC','AAG','AAT','ACA','ACC','ACT','ACG','ATA','ATC','ATT','ATG','AGA','AGC','AGT','AGG'/
      
      data aset3(:,1) /64*'A   '/
      data aset3(:,2) /'AAAA','AAAC','AAAG','AAAT','AACA','AACC','AACT','AACG','AATA','AATC','AATT','AATG','AAGA','AAGC','AAGT','AAGG', &
                       'ACAA','ACAC','ACAG','ACAT','ACCA','ACCC','ACCT','ACCG','ACTA','ACTC','ACTT','ACTG','ACGA','ACGC','ACGT','ACGG', &
                       'ATAA','ATAC','ATAG','ATAT','ATCA','ATCC','ATCT','ATCG','ATTA','ATTC','ATTT','ATTG','ATGA','ATGC','ATGT','ATGG', &
                       'AGAA','AGAC','AGAG','AGAT','AGCA','AGCC','AGCT','AGCG','AGTA','AGTC','AGTT','AGTG','AGGA','AGGC','AGGT','AGGG'/

      data aset4(:,1) /256*'A    '/
      data aset4(:,2) /'AAAAA','AAAAC','AAAAG','AAAAT','AAACA','AAACC','AAACT','AAACG','AAATA','AAATC','AAATT','AAATG','AAAGA','AAAGC','AAAGT','AAAGG', &
                       'AACAA','AACAC','AACAG','AACAT','AACCA','AACCC','AACCT','AACCG','AACTA','AACTC','AACTT','AACTG','AACGA','AACGC','AACGT','AACGG', &
                       'AATAA','AATAC','AATAG','AATAT','AATCA','AATCC','AATCT','AATCG','AATTA','AATTC','AATTT','AATTG','AATGA','AATGC','AATGT','AATGG', &
                       'AAGAA','AAGAC','AAGAG','AAGAT','AAGCA','AAGCC','AAGCT','AAGCG','AAGTA','AAGTC','AAGTT','AAGTG','AAGGA','AAGGC','AAGGT','AAGGG', &
                       'ACAAA','ACAAC','ACAAG','ACAAT','ACACA','ACACC','ACACT','ACACG','ACATA','ACATC','ACATT','ACATG','ACAGA','ACAGC','ACAGT','ACAGG', &
                       'ACCAA','ACCAC','ACCAG','ACCAT','ACCCA','ACCCC','ACCCT','ACCCG','ACCTA','ACCTC','ACCTT','ACCTG','ACCGA','ACCGC','ACCGT','ACCGG', &
                       'ACTAA','ACTAC','ACTAG','ACTAT','ACTCA','ACTCC','ACTCT','ACTCG','ACTTA','ACTTC','ACTTT','ACTTG','ACTGA','ACTGC','ACTGT','ACTGG', &
                       'ACGAA','ACGAC','ACGAG','ACGAT','ACGCA','ACGCC','ACGCT','ACGCG','ACGTA','ACGTC','ACGTT','ACGTG','ACGGA','ACGGC','ACGGT','ACGGG', &
                       'ATAAA','ATAAC','ATAAG','ATAAT','ATACA','ATACC','ATACT','ATACG','ATATA','ATATC','ATATT','ATATG','ATAGA','ATAGC','ATAGT','ATAGG', &
                       'ATCAA','ATCAC','ATCAG','ATCAT','ATCCA','ATCCC','ATCCT','ATCCG','ATCTA','ATCTC','ATCTT','ATCTG','ATCGA','ATCGC','ATCGT','ATCGG', &
                       'ATTAA','ATTAC','ATTAG','ATTAT','ATTCA','ATTCC','ATTCT','ATTCG','ATTTA','ATTTC','ATTTT','ATTTG','ATTGA','ATTGC','ATTGT','ATTGG', &
                       'ATGAA','ATGAC','ATGAG','ATGAT','ATGCA','ATGCC','ATGCT','ATGCG','ATGTA','ATGTC','ATGTT','ATGTG','ATGGA','ATGGC','ATGGT','ATGGG', &
                       'AGAAA','AGAAC','AGAAG','AGAAT','AGACA','AGACC','AGACT','AGACG','AGATA','AGATC','AGATT','AGATG','AGAGA','AGAGC','AGAGT','AGAGG', &
                       'AGCAA','AGCAC','AGCAG','AGCAT','AGCCA','AGCCC','AGCCT','AGCCG','AGCTA','AGCTC','AGCTT','AGCTG','AGCGA','AGCGC','AGCGT','AGCGG', &
                       'AGTAA','AGTAC','AGTAG','AGTAT','AGTCA','AGTCC','AGTCT','AGTCG','AGTTA','AGTTC','AGTTT','AGTTG','AGTGA','AGTGC','AGTGT','AGTGG', &
                       'AGGAA','AGGAC','AGGAG','AGGAT','AGGCA','AGGCC','AGGCT','AGGCG','AGGTA','AGGTC','AGGTT','AGGTG','AGGGA','AGGGC','AGGGT','AGGGG'/

      data aset5(:,1) /1024*'A     '/
      data aset5(:,2) /'AAAAAA','AAAAAC','AAAAAG','AAAAAT','AAAACA','AAAACC','AAAACT','AAAACG','AAAATA','AAAATC','AAAATT','AAAATG','AAAAGA','AAAAGC','AAAAGT','AAAAGG', &
                       'AAACAA','AAACAC','AAACAG','AAACAT','AAACCA','AAACCC','AAACCT','AAACCG','AAACTA','AAACTC','AAACTT','AAACTG','AAACGA','AAACGC','AAACGT','AAACGG', &
                       'AAATAA','AAATAC','AAATAG','AAATAT','AAATCA','AAATCC','AAATCT','AAATCG','AAATTA','AAATTC','AAATTT','AAATTG','AAATGA','AAATGC','AAATGT','AAATGG', &
                       'AAAGAA','AAAGAC','AAAGAG','AAAGAT','AAAGCA','AAAGCC','AAAGCT','AAAGCG','AAAGTA','AAAGTC','AAAGTT','AAAGTG','AAAGGA','AAAGGC','AAAGGT','AAAGGG', &
                       'AACAAA','AACAAC','AACAAG','AACAAT','AACACA','AACACC','AACACT','AACACG','AACATA','AACATC','AACATT','AACATG','AACAGA','AACAGC','AACAGT','AACAGG', &
                       'AACCAA','AACCAC','AACCAG','AACCAT','AACCCA','AACCCC','AACCCT','AACCCG','AACCTA','AACCTC','AACCTT','AACCTG','AACCGA','AACCGC','AACCGT','AACCGG', &
                       'AACTAA','AACTAC','AACTAG','AACTAT','AACTCA','AACTCC','AACTCT','AACTCG','AACTTA','AACTTC','AACTTT','AACTTG','AACTGA','AACTGC','AACTGT','AACTGG', &
                       'AACGAA','AACGAC','AACGAG','AACGAT','AACGCA','AACGCC','AACGCT','AACGCG','AACGTA','AACGTC','AACGTT','AACGTG','AACGGA','AACGGC','AACGGT','AACGGG', &
                       'AATAAA','AATAAC','AATAAG','AATAAT','AATACA','AATACC','AATACT','AATACG','AATATA','AATATC','AATATT','AATATG','AATAGA','AATAGC','AATAGT','AATAGG', &
                       'AATCAA','AATCAC','AATCAG','AATCAT','AATCCA','AATCCC','AATCCT','AATCCG','AATCTA','AATCTC','AATCTT','AATCTG','AATCGA','AATCGC','AATCGT','AATCGG', &
                       'AATTAA','AATTAC','AATTAG','AATTAT','AATTCA','AATTCC','AATTCT','AATTCG','AATTTA','AATTTC','AATTTT','AATTTG','AATTGA','AATTGC','AATTGT','AATTGG', &
                       'AATGAA','AATGAC','AATGAG','AATGAT','AATGCA','AATGCC','AATGCT','AATGCG','AATGTA','AATGTC','AATGTT','AATGTG','AATGGA','AATGGC','AATGGT','AATGGG', &
                       'AAGAAA','AAGAAC','AAGAAG','AAGAAT','AAGACA','AAGACC','AAGACT','AAGACG','AAGATA','AAGATC','AAGATT','AAGATG','AAGAGA','AAGAGC','AAGAGT','AAGAGG', &
                       'AAGCAA','AAGCAC','AAGCAG','AAGCAT','AAGCCA','AAGCCC','AAGCCT','AAGCCG','AAGCTA','AAGCTC','AAGCTT','AAGCTG','AAGCGA','AAGCGC','AAGCGT','AAGCGG', &
                       'AAGTAA','AAGTAC','AAGTAG','AAGTAT','AAGTCA','AAGTCC','AAGTCT','AAGTCG','AAGTTA','AAGTTC','AAGTTT','AAGTTG','AAGTGA','AAGTGC','AAGTGT','AAGTGG', &
                       'AAGGAA','AAGGAC','AAGGAG','AAGGAT','AAGGCA','AAGGCC','AAGGCT','AAGGCG','AAGGTA','AAGGTC','AAGGTT','AAGGTG','AAGGGA','AAGGGC','AAGGGT','AAGGGG', &
                       'ACAAAA','ACAAAC','ACAAAG','ACAAAT','ACAACA','ACAACC','ACAACT','ACAACG','ACAATA','ACAATC','ACAATT','ACAATG','ACAAGA','ACAAGC','ACAAGT','ACAAGG', &
                       'ACACAA','ACACAC','ACACAG','ACACAT','ACACCA','ACACCC','ACACCT','ACACCG','ACACTA','ACACTC','ACACTT','ACACTG','ACACGA','ACACGC','ACACGT','ACACGG', &
                       'ACATAA','ACATAC','ACATAG','ACATAT','ACATCA','ACATCC','ACATCT','ACATCG','ACATTA','ACATTC','ACATTT','ACATTG','ACATGA','ACATGC','ACATGT','ACATGG', &
                       'ACAGAA','ACAGAC','ACAGAG','ACAGAT','ACAGCA','ACAGCC','ACAGCT','ACAGCG','ACAGTA','ACAGTC','ACAGTT','ACAGTG','ACAGGA','ACAGGC','ACAGGT','ACAGGG', &
                       'ACCAAA','ACCAAC','ACCAAG','ACCAAT','ACCACA','ACCACC','ACCACT','ACCACG','ACCATA','ACCATC','ACCATT','ACCATG','ACCAGA','ACCAGC','ACCAGT','ACCAGG', &
                       'ACCCAA','ACCCAC','ACCCAG','ACCCAT','ACCCCA','ACCCCC','ACCCCT','ACCCCG','ACCCTA','ACCCTC','ACCCTT','ACCCTG','ACCCGA','ACCCGC','ACCCGT','ACCCGG', &
                       'ACCTAA','ACCTAC','ACCTAG','ACCTAT','ACCTCA','ACCTCC','ACCTCT','ACCTCG','ACCTTA','ACCTTC','ACCTTT','ACCTTG','ACCTGA','ACCTGC','ACCTGT','ACCTGG', &
                       'ACCGAA','ACCGAC','ACCGAG','ACCGAT','ACCGCA','ACCGCC','ACCGCT','ACCGCG','ACCGTA','ACCGTC','ACCGTT','ACCGTG','ACCGGA','ACCGGC','ACCGGT','ACCGGG', &
                       'ACTAAA','ACTAAC','ACTAAG','ACTAAT','ACTACA','ACTACC','ACTACT','ACTACG','ACTATA','ACTATC','ACTATT','ACTATG','ACTAGA','ACTAGC','ACTAGT','ACTAGG', &
                       'ACTCAA','ACTCAC','ACTCAG','ACTCAT','ACTCCA','ACTCCC','ACTCCT','ACTCCG','ACTCTA','ACTCTC','ACTCTT','ACTCTG','ACTCGA','ACTCGC','ACTCGT','ACTCGG', &
                       'ACTTAA','ACTTAC','ACTTAG','ACTTAT','ACTTCA','ACTTCC','ACTTCT','ACTTCG','ACTTTA','ACTTTC','ACTTTT','ACTTTG','ACTTGA','ACTTGC','ACTTGT','ACTTGG', &
                       'ACTGAA','ACTGAC','ACTGAG','ACTGAT','ACTGCA','ACTGCC','ACTGCT','ACTGCG','ACTGTA','ACTGTC','ACTGTT','ACTGTG','ACTGGA','ACTGGC','ACTGGT','ACTGGG', &
                       'ACGAAA','ACGAAC','ACGAAG','ACGAAT','ACGACA','ACGACC','ACGACT','ACGACG','ACGATA','ACGATC','ACGATT','ACGATG','ACGAGA','ACGAGC','ACGAGT','ACGAGG', &
                       'ACGCAA','ACGCAC','ACGCAG','ACGCAT','ACGCCA','ACGCCC','ACGCCT','ACGCCG','ACGCTA','ACGCTC','ACGCTT','ACGCTG','ACGCGA','ACGCGC','ACGCGT','ACGCGG', &
                       'ACGTAA','ACGTAC','ACGTAG','ACGTAT','ACGTCA','ACGTCC','ACGTCT','ACGTCG','ACGTTA','ACGTTC','ACGTTT','ACGTTG','ACGTGA','ACGTGC','ACGTGT','ACGTGG', &
                       'ACGGAA','ACGGAC','ACGGAG','ACGGAT','ACGGCA','ACGGCC','ACGGCT','ACGGCG','ACGGTA','ACGGTC','ACGGTT','ACGGTG','ACGGGA','ACGGGC','ACGGGT','ACGGGG', &
                       'ATAAAA','ATAAAC','ATAAAG','ATAAAT','ATAACA','ATAACC','ATAACT','ATAACG','ATAATA','ATAATC','ATAATT','ATAATG','ATAAGA','ATAAGC','ATAAGT','ATAAGG', &
                       'ATACAA','ATACAC','ATACAG','ATACAT','ATACCA','ATACCC','ATACCT','ATACCG','ATACTA','ATACTC','ATACTT','ATACTG','ATACGA','ATACGC','ATACGT','ATACGG', &
                       'ATATAA','ATATAC','ATATAG','ATATAT','ATATCA','ATATCC','ATATCT','ATATCG','ATATTA','ATATTC','ATATTT','ATATTG','ATATGA','ATATGC','ATATGT','ATATGG', &
                       'ATAGAA','ATAGAC','ATAGAG','ATAGAT','ATAGCA','ATAGCC','ATAGCT','ATAGCG','ATAGTA','ATAGTC','ATAGTT','ATAGTG','ATAGGA','ATAGGC','ATAGGT','ATAGGG', &
                       'ATCAAA','ATCAAC','ATCAAG','ATCAAT','ATCACA','ATCACC','ATCACT','ATCACG','ATCATA','ATCATC','ATCATT','ATCATG','ATCAGA','ATCAGC','ATCAGT','ATCAGG', &
                       'ATCCAA','ATCCAC','ATCCAG','ATCCAT','ATCCCA','ATCCCC','ATCCCT','ATCCCG','ATCCTA','ATCCTC','ATCCTT','ATCCTG','ATCCGA','ATCCGC','ATCCGT','ATCCGG', &
                       'ATCTAA','ATCTAC','ATCTAG','ATCTAT','ATCTCA','ATCTCC','ATCTCT','ATCTCG','ATCTTA','ATCTTC','ATCTTT','ATCTTG','ATCTGA','ATCTGC','ATCTGT','ATCTGG', &
                       'ATCGAA','ATCGAC','ATCGAG','ATCGAT','ATCGCA','ATCGCC','ATCGCT','ATCGCG','ATCGTA','ATCGTC','ATCGTT','ATCGTG','ATCGGA','ATCGGC','ATCGGT','ATCGGG', &
                       'ATTAAA','ATTAAC','ATTAAG','ATTAAT','ATTACA','ATTACC','ATTACT','ATTACG','ATTATA','ATTATC','ATTATT','ATTATG','ATTAGA','ATTAGC','ATTAGT','ATTAGG', &
                       'ATTCAA','ATTCAC','ATTCAG','ATTCAT','ATTCCA','ATTCCC','ATTCCT','ATTCCG','ATTCTA','ATTCTC','ATTCTT','ATTCTG','ATTCGA','ATTCGC','ATTCGT','ATTCGG', &
                       'ATTTAA','ATTTAC','ATTTAG','ATTTAT','ATTTCA','ATTTCC','ATTTCT','ATTTCG','ATTTTA','ATTTTC','ATTTTT','ATTTTG','ATTTGA','ATTTGC','ATTTGT','ATTTGG', &
                       'ATTGAA','ATTGAC','ATTGAG','ATTGAT','ATTGCA','ATTGCC','ATTGCT','ATTGCG','ATTGTA','ATTGTC','ATTGTT','ATTGTG','ATTGGA','ATTGGC','ATTGGT','ATTGGG', &
                       'ATGAAA','ATGAAC','ATGAAG','ATGAAT','ATGACA','ATGACC','ATGACT','ATGACG','ATGATA','ATGATC','ATGATT','ATGATG','ATGAGA','ATGAGC','ATGAGT','ATGAGG', &
                       'ATGCAA','ATGCAC','ATGCAG','ATGCAT','ATGCCA','ATGCCC','ATGCCT','ATGCCG','ATGCTA','ATGCTC','ATGCTT','ATGCTG','ATGCGA','ATGCGC','ATGCGT','ATGCGG', &
                       'ATGTAA','ATGTAC','ATGTAG','ATGTAT','ATGTCA','ATGTCC','ATGTCT','ATGTCG','ATGTTA','ATGTTC','ATGTTT','ATGTTG','ATGTGA','ATGTGC','ATGTGT','ATGTGG', &
                       'ATGGAA','ATGGAC','ATGGAG','ATGGAT','ATGGCA','ATGGCC','ATGGCT','ATGGCG','ATGGTA','ATGGTC','ATGGTT','ATGGTG','ATGGGA','ATGGGC','ATGGGT','ATGGGG', &
                       'AGAAAA','AGAAAC','AGAAAG','AGAAAT','AGAACA','AGAACC','AGAACT','AGAACG','AGAATA','AGAATC','AGAATT','AGAATG','AGAAGA','AGAAGC','AGAAGT','AGAAGG', &
                       'AGACAA','AGACAC','AGACAG','AGACAT','AGACCA','AGACCC','AGACCT','AGACCG','AGACTA','AGACTC','AGACTT','AGACTG','AGACGA','AGACGC','AGACGT','AGACGG', &
                       'AGATAA','AGATAC','AGATAG','AGATAT','AGATCA','AGATCC','AGATCT','AGATCG','AGATTA','AGATTC','AGATTT','AGATTG','AGATGA','AGATGC','AGATGT','AGATGG', &
                       'AGAGAA','AGAGAC','AGAGAG','AGAGAT','AGAGCA','AGAGCC','AGAGCT','AGAGCG','AGAGTA','AGAGTC','AGAGTT','AGAGTG','AGAGGA','AGAGGC','AGAGGT','AGAGGG', &
                       'AGCAAA','AGCAAC','AGCAAG','AGCAAT','AGCACA','AGCACC','AGCACT','AGCACG','AGCATA','AGCATC','AGCATT','AGCATG','AGCAGA','AGCAGC','AGCAGT','AGCAGG', &
                       'AGCCAA','AGCCAC','AGCCAG','AGCCAT','AGCCCA','AGCCCC','AGCCCT','AGCCCG','AGCCTA','AGCCTC','AGCCTT','AGCCTG','AGCCGA','AGCCGC','AGCCGT','AGCCGG', &
                       'AGCTAA','AGCTAC','AGCTAG','AGCTAT','AGCTCA','AGCTCC','AGCTCT','AGCTCG','AGCTTA','AGCTTC','AGCTTT','AGCTTG','AGCTGA','AGCTGC','AGCTGT','AGCTGG', &
                       'AGCGAA','AGCGAC','AGCGAG','AGCGAT','AGCGCA','AGCGCC','AGCGCT','AGCGCG','AGCGTA','AGCGTC','AGCGTT','AGCGTG','AGCGGA','AGCGGC','AGCGGT','AGCGGG', &
                       'AGTAAA','AGTAAC','AGTAAG','AGTAAT','AGTACA','AGTACC','AGTACT','AGTACG','AGTATA','AGTATC','AGTATT','AGTATG','AGTAGA','AGTAGC','AGTAGT','AGTAGG', &
                       'AGTCAA','AGTCAC','AGTCAG','AGTCAT','AGTCCA','AGTCCC','AGTCCT','AGTCCG','AGTCTA','AGTCTC','AGTCTT','AGTCTG','AGTCGA','AGTCGC','AGTCGT','AGTCGG', &
                       'AGTTAA','AGTTAC','AGTTAG','AGTTAT','AGTTCA','AGTTCC','AGTTCT','AGTTCG','AGTTTA','AGTTTC','AGTTTT','AGTTTG','AGTTGA','AGTTGC','AGTTGT','AGTTGG', &
                       'AGTGAA','AGTGAC','AGTGAG','AGTGAT','AGTGCA','AGTGCC','AGTGCT','AGTGCG','AGTGTA','AGTGTC','AGTGTT','AGTGTG','AGTGGA','AGTGGC','AGTGGT','AGTGGG', &
                       'AGGAAA','AGGAAC','AGGAAG','AGGAAT','AGGACA','AGGACC','AGGACT','AGGACG','AGGATA','AGGATC','AGGATT','AGGATG','AGGAGA','AGGAGC','AGGAGT','AGGAGG', &
                       'AGGCAA','AGGCAC','AGGCAG','AGGCAT','AGGCCA','AGGCCC','AGGCCT','AGGCCG','AGGCTA','AGGCTC','AGGCTT','AGGCTG','AGGCGA','AGGCGC','AGGCGT','AGGCGG', &
                       'AGGTAA','AGGTAC','AGGTAG','AGGTAT','AGGTCA','AGGTCC','AGGTCT','AGGTCG','AGGTTA','AGGTTC','AGGTTT','AGGTTG','AGGTGA','AGGTGC','AGGTGT','AGGTGG', &
                       'AGGGAA','AGGGAC','AGGGAG','AGGGAT','AGGGCA','AGGGCC','AGGGCT','AGGGCG','AGGGTA','AGGGTC','AGGGTT','AGGGTG','AGGGGA','AGGGGC','AGGGGT','AGGGGG'/

!     Initialize program parameters to default values

      fnmarker='test.markers'
      fnfam='test.fam'
      fngtype='test.phased'
      fnvcf='test.vcf'
      chr=6
      missing='0'
      recode_opt=1
      revert_opt=2
      answer(1)=' (1) Yes'
      answer(2)=' (2) No'
      accept=''
      
!     Query user for starting parameters.  If no arguments on command line, bring up
!     menu for interactive input.  Otherwise, parse command line in Unix-style to determine
!     starting parameters.

      if (command_argument_count()>0) then
         do i=1,command_argument_count(),2
           call get_command_argument(i,arg)
           select case (arg)
             case('-fnmarker')
               call get_command_argument(i+1,fnmarker)
	           inquire(file=fnmarker,exist=exist)
               if (.not.exist) then
                 write(*,'(/a,a,a)') ' Input file ',trim(fnmarker),' not in directory'
                 stop
               end if
             case('-fnfam')
               call get_command_argument(i+1,fnfam)
	           inquire(file=fnfam,exist=exist)
               if (.not.exist) then
                 write(*,'(/a,a,a)') ' Input file ',trim(fnfam),' not in directory'
                 stop
               end if
             case('-fngtype')
               call get_command_argument(i+1,fngtype)
	           inquire(file=fngtype,exist=exist)
               if (.not.exist) then
                 write(*,'(/a,a,a)') ' Input file ',trim(fngtype),' not in directory'
                 stop
               end if
             case('-fnvcf')
               call get_command_argument(i+1,fnvcf)
             case('-missing')
                call get_command_argument(i+1,string)
             case('-chr')
               call get_command_argument(i+1,string)
               read (string,'(i2)') chr
             case('-recode_opt','-recode')
               call get_command_argument(i+1,string)
               read (string,'(i1)') recode_opt
             case('-revert_opt')
               call get_command_argument(i+1,string)
               read (string,'(i1)') revert_opt
             case default
               write (*,'(a,a,/)') 'Unrecognized command-line option: ', arg
               stop
           end select
         end do
       else

!        Print menu on console allowing user to enter input information 
!        and select program options.  Keep cycling through menu until 
!        user indicates (with a '21') that all selections have been made.

         do while (item/=21)

           nskip=14
           write(*,9110) fnmarker,fnfam,fngtype,missing,fnvcf,chr,answer(recode_opt),answer(revert_opt)
           write(*,9190) accept

9110       format ('  1  MARKER INPUT FILE',T45,A35/ &	
                   '  2  FAMILY INPUT FILE',T45,A35/ & 
                   '  3  GENOTYPE INPUT FILE',T45,A35/ &
                   '  4  MISSING ALLELE CODE (GTYPE FILE)',T45,A1/ &
			       '  5  VCF OUTPUT FILE',T45,A35/ &
                   '  6  CHROMOSOME ID IN VCF OUTPUT',T45,I2/ &
                   '  7  RECODE NON-ACTG & DUPL MARKERS (1-2)?',T45,A35/ &
                   '  8  ARRAY SNPs: KEEP OLD ALLELES (1-2)?',T45,A35/)
9190      format  (' 21  ACCEPT PROGRAM',T45,A1)

           do i=1,nskip
             write(*,*)
           end do
           write(*,9195,advance='NO')
9195       format (' CHOOSE AN ITEM [1.....21]: ')

!          Read in user's entry for selected menu item.

            read(*,*) item
            select case (item)
              case (1)
                write (*,9201,ADVANCE='NO') 
9201            format (' Name of marker input file (Beagle 3 format): ')
                read (*,'(a)') fnmarker
	            inquire(FILE=fnmarker,EXIST=EXIST)
	            do while (.not.exist)
	              write (*,9201,ADVANCE='NO')
	              read (*,'(a)') fnmarker
	              inquire(FILE=fnmarker,EXIST=EXIST)
	            end do
              case (2)
                write (*,9202,ADVANCE='NO') 
9202            format (' Name of family input file (linkage format): ')
                read (*,'(a)') fnfam
	            inquire(FILE=fnfam,EXIST=EXIST)
	            do while (.not.exist)
	              write (*,9202,ADVANCE='NO')
	              read (*,'(a)') fnfam
	              inquire(FILE=fnfam,EXIST=EXIST)
	            end do
              case (3)
                write (*,9203,ADVANCE='NO') 
9203            format (' Name of phased genotype input file (Beagle 3 format): ')
                read (*,'(a)') fngtype
	            inquire(FILE=fngtype,EXIST=EXIST)
	            do while (.not.exist)
	              write (*,9203,ADVANCE='NO')
	              read (*,'(a)') fngtype
	              inquire(FILE=fngtype,EXIST=EXIST)
	            end do
              case (4)
	            write (*,9204,ADVANCE='NO') 
9204            format (' Code for missing allele in Beagle phased genotype file: ' )
	            read (*,*) missing
              case (5)
                write (*,9205,ADVANCE='NO') 
9205            format (' Name of phased vcf output file: ')
                read (*,'(a)') fnvcf
              case (6)
	            write (*,9206,ADVANCE='NO') 
9206            format (' Chromosome identifier for variants in output vcf file (1-21): ' )
	            read (*,*) chr
                do while (chr<1.or.chr>23)
                  write (*,9206,ADVANCE='NO')
	              read (*,*) chr
                end do
              case (7)
	            write (*,9207,ADVANCE='NO') 
9207            format (' Recode markers with non-ACTG alleles? (1-2): ' )
	            read (*,*) recode_opt
                do while (recode_opt/=1.and.recode_opt/=2)
                  write (*,9207,ADVANCE='NO')
	              read (*,*) recode_opt
                end do
              case (8)
	            write (*,9208,ADVANCE='NO') 
9208            format (' Keep original alleles for microarray-genotyped SNPs? (1-2): ' )
	            read (*,*) revert_opt
                do while (revert_opt/=1.and.revert_opt/=2)
                  write (*,9208,ADVANCE='NO')
	              read (*,*) revert_opt
                end do
              case (21)
                exit
            end select
         end do
      endif

      write(*,'(/A//)') ' tped2vcf v 1.0'
      
!     Open family input file, and determine number of samples.
      
      open (unit=15,file=fnfam)
      
      nsample=0
      do
        read(15,*,iostat=ioerr) ctmp
        if (ioerr<0) exit
        nsample=nsample+1
      end do
      
!     Allocate dynamic array, rewind file, and read in all sample IDs
      
      allocate (samples(nsample))
      rewind(15)
      do i=1,nsample
        read(15,*) ctmp,samples(i)
      end do
      
      close (15)
      
!     Open marker input file, and determine number of target markers.
      
      open (unit=15,file=fnmarker)
      
      nmarker=0
      do
        read(15,*,iostat=ioerr) ctmp
        if (ioerr<0) exit
        nmarker=nmarker+1
      end do
      
!     Rewind marker file.  If user requested recoding of alleles, read
!     into arrays position and allele identities for all markers, and then
!     dedupulicate so-called identical markers by giving each within a duplicate
!     set its own unique reference and alternate alleles, following this
!     by recoding any remaining markers with non-ACTG alleles as A/AA SNPs.
!     Otherwise, simply determine postion of first and last marker in file.
      
      rewind(15)
      
      if (recode_opt==1) then
          
!       First read in positions and alleles for all markers
          
        allocate (position(nmarker),ralleles(nmarker,2))
        do i=1,nmarker
          read(15,*) ctmp,position(i),ralleles(i,1),ralleles(i,2)
        end do
        pos_first=position(1)
        pos_last=position(nmarker)
        if (revert_opt==1) then
          allocate(ralleles_original(nmarker,2))
          ralleles_original=ralleles
        end if
        
!       Now deduplicate apparently identical markers (i.e., identical bp positions and allele identities)
        
        i=0
        do
          i=i+1
          if (i>nmarker) exit
          if (count(position==position(i).and.ralleles(:,1)==ralleles(i,1).and.ralleles(:,2)==ralleles(i,2))>1) then
            istart=i
            do j=nmarker,istart,-1
              if (position(j)==position(i).and.ralleles(j,1)==ralleles(i,1).and.ralleles(j,2)==ralleles(i,2)) then
                iend=j
                exit
              end if
            end do
            ndup=iend-istart+1
            select case (ndup)
              case (2:4)
                ralleles(istart:iend,:)=aset1(1:ndup,:)
              case (5:16)
                ralleles(istart:iend,:)=aset2(1:ndup,:)
              case (17:64)
                ralleles(istart:iend,:)=aset3(1:ndup,:)
              case (65:256)
                ralleles(istart:iend,:)=aset4(1:ndup,:)
              case (257:1024)
                ralleles(istart:iend,:)=aset5(1:ndup,:)
              case default
                write(*,*) ' Too many markers with same position'
                stop
            end select
            i=i+ndup-1
          end if
        end do
!        deallocate (position)
        
!       Recode any remaining markers with non-ACTG alleles as A/AA SNPs

        do i=1,nmarker
          if (scan(ralleles(i,1),'ACTG0')==0.or.scan(ralleles(i,2),'ACTG0')==0) then
            ralleles(i,1)='A '
            ralleles(i,2)='AA'
          end if
        end do

!       Repeat deduplication of apparently identical markers (i.e., identical bp positions and allele identities),
!       because previous step to recode remaining markers with non-ACTG alleles as A/AA SNPs
!       in some instances may have created new sets of apparently identical markers (e.g., DRB3/4/5 AA and SNP
!       markers with the same position where one allele is a Z = absence of gene will usually not appear as duplicates 
!       to start with because the non-Z alleles are likely different, but the previous step will convert the alleles
!       to A/AA because of the non-ACTG Z.
        
        i=0
        do
          i=i+1
          if (i>nmarker) exit
          if (count(position==position(i).and.ralleles(:,1)==ralleles(i,1).and.ralleles(:,2)==ralleles(i,2))>1) then
            istart=i
            do j=nmarker,istart,-1
              if (position(j)==position(i).and.ralleles(j,1)==ralleles(i,1).and.ralleles(j,2)==ralleles(i,2)) then
                iend=j
                exit
              end if
            end do
            ndup=iend-istart+1
            select case (ndup)
              case (2:4)
                ralleles(istart:iend,:)=aset1(1:ndup,:)
              case (5:16)
                ralleles(istart:iend,:)=aset2(1:ndup,:)
              case (17:64)
                ralleles(istart:iend,:)=aset3(1:ndup,:)
              case (65:256)
                ralleles(istart:iend,:)=aset4(1:ndup,:)
              case (257:1024)
                ralleles(istart:iend,:)=aset5(1:ndup,:)
              case default
                write(*,*) ' Too many markers with same position'
                stop
            end select
            i=i+ndup-1
          end if
        end do
        deallocate (position)
      
!       If requested by user, revert alleles for microarray-genotyped variants to their original values,
!       unless original alleles have non-ACTG characters.
      
        if (revert_opt==1) then
          rewind(15)
          do i=1,nmarker
            read(15,*) mkid
            if (mkid(1:4)/='SNP_'.and.mkid(1:4)/='HLA_'.and.mkid(1:3)/='AA_'.and.mkid(1:4)/='INS_'.and. &
                scan(ralleles_original(i,1),'ACTG0')/=0.and.scan(ralleles_original(i,2),'ACTG0')/=0) then
              ralleles(i,:)=ralleles_original(i,:)
            end if
          end do
          deallocate (ralleles_original)
        end if
      else
        do i=1,nmarker
          read(15,*) ctmp,item
          if (i==1) pos_first=item
          if (i==nmarker) pos_last=item
        end do
      end if
      
      close (15)
      
!     Write header information at top of vcf output file.  Note the use of tab-delimited output 
!     for the header line, in accordance with vcf 4.2 standards.
      
      open (unit=17, file=fnvcf)
      call DATE_AND_TIME(date)
      write(17,9300) date,chr,pos_last-pos_first+1
9300  format ('##fileformat=VCFv4.2'/ &
              '##fileDate=',a8/ &
              '##source=beagle2vcf.f90_v1.2'/ &
              '##contig=<ID=',i2,',length=',i9,'>'/ &
              '##INFO=<ID=PR,Number=0,Type=Flag,Description="Provisional reference allele, may not be based on real reference genome">'/ &
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
       write(17,'(*(g0,:,"'//achar(9)//'"))') '#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',(trim(samples(i)),i=1,nsample)
       
!     Set genotype allele delimiter to "/" to denote unphased genotype data
       
      delimiter='/'
      
!     Open both input files.  
      
      open (unit=15,file=fnmarker)
      open (unit=16,file=fngtype)
      
!     Read in one record at a time from both files.  Combine information from the
!     two records to create appropriate record to be output to vcf file.  Note the
!     use of tab-delimited output in accordance with vcf 4.2 standards.
      
      allocate (alleles(nsample,2))
      
      do imk=1,nmarker

        read(15,*) mkid1,i_position,a1,a2
        read(16,*) ctmp,mkid2,ctmp,ctmp,((alleles(i,j),j=1,2),i=1,nsample)
        
        if (mkid1/=mkid2) then
          write (*,9310) imk 
9310      format (' Marker number ',i6,' differs in two input files.  Halting program.')
          stop  
        end if

        where(alleles==a1) 
          alleles='0'
        elsewhere(alleles==a2)
          alleles='1'
        elsewhere(alleles=='0')
          alleles='.'
        elsewhere
!         Force invalid allele values to missing.
          alleles='.'
        end where
        
        if (recode_opt==1) then
          write(17,'(*(g0,:,"'//achar(9)//'"))') chr,i_position,trim(mkid1),trim(ralleles(imk,1)),trim(ralleles(imk,2)),'.','.','PR','GT',(trim(alleles(i,1))//delimiter//trim(alleles(i,2)),i=1,nsample)
        else
          write(17,'(*(g0,:,"'//achar(9)//'"))') chr,i_position,trim(mkid1),trim(a1),trim(a2),'.','.','PR','GT',(trim(alleles(i,1))//delimiter//trim(alleles(i,2)),i=1,nsample)
        end if          
        
      end do
      
      close(15); close(16); close(17)

   end program tped2vcf