    program purge_SNPs
!--------------------------------------------------------------------------------------
!     Program purge_SNPs is a utility to be used with my modification of the 
!     MakeReference program of SNP2HLA (Jia 2013; SNP2HLA software).
!
!     As part of my modification of SNP2HLA, I updated the amino acid and SNP data
!     dictionaries to reflect the current (r3.27.0) data release of the IPD-IMGT/HLA 
!     database.  Compared to the original data dictionaries supplied with SNP2HLA, 
!     this added many new alleles and filled in some of the missing portions of 
!     protein and DNA sequence for previously known alleles.  I also corrected two
!     issues with the original SNP dictionary: (1) DNA sequence for a four-digit
!     allele is now listed only if at least one of the six or eight-digit subtypes of
!     that four-digit allele has sequence present in the IMGT/HLA database (instead
!     of arbitrarily filling in unknown sequence with another four-digit allele 
!     subtype of the same two-digit class; (2) ambiguous nucleotide codes are now
!     used to represent positions where two or more different nucleotides occured 
!     among all six and eight-digit subtypes (instead of arbitrarily using the
!     sequence of just one of the subtypes).
!
!     My correction of these two situations can lead to the calling of many
!     intragenic HLA SNP variants where at least one of the alleles for the variant
!     is ambigous or missing.  Since the reference panel is usually smaller than
!     the data being imputed, most SNPs with missing or ambiguous alleles in the 
!     reference will also have missing or ambigous alleles in the imputed dataset,
!     which is not ideal for downstream analysis.  This program creates a list of
!     these problematic SNPs so they can be filtered out during creation of the
!     reference panel.
!
!     The required input file is the Plink file $OUTPUT.SNPS.CODED.bim created
!     by the MakeReference C-shell script.  This file must have the six columns of
!     data expected for Plink *.bim format (i.e., chrom, SNP_ID, recomb_distance, 
!     bp_position, allele1, allele2).  The output file is a text file with a list
!     of the IDs for all SNPs to be excluded.
!     
!     Written by Philip Stuart in standard Fortran 2003.                 
!                                                                      
!     Revision History:												 
!                                                                      
!     v 1.0  04/04/17   First version.
!
!     v 1.1  08/28/19   Increased length of character scalars fnbim and fnout 
!                       from 250 to 1000 and arg from 251 to 1001; also 
!                       increased length of allocatale character array snp from
!                       21 to 30.
!
!     References:
!
!     Jia X, Han B, Onengut-Gumuscu S, Chen W-M, Concannon PJ, Rich SS,
!        Raychaudhuri S, de Bakker PIW (2013)Imputing amino acid polymorphisms in
!       humanleukocyte antigens, PLoS One 8(6):e64683.
!
!     Robinson J, Halliwell JA, Hayhurst JH, Flicek P, Parham P, Marsh SGE (2015)
!        The IPD and IPD-IMGT/HLA Database: allele variant databases. Nucleic
!        Acids Research 43:D423-431.
!
!     SNP2HLA software, version 1.0.3, downloaded in May 2016 from 
!        https://www.broadinstitute.org/impg/snp2hla
!--------------------------------------------------------------------------------------
      implicit none

      integer i,j,ioerr,ic1,ic2,ns,nv,nvsm,dist1,dist2
      integer, allocatable :: nvs(:)
      character(1) atmp
      character(1), allocatable :: a1(:), a2(:)
      character(30), allocatable :: snp(:)
      character(1000) fnbim,fnout
      character(1001) arg
      character(1000) set_nt
      logical exist
      
!     Initialize program parameters to default values

      fnbim='HAPMAP_CEU_REF_test3.SNPS.CODED.bim'
      fnout='to_remove2'
      
!     Parse command line to determine starting parameters.

      do i=1,command_argument_count(),2
        call get_command_argument(i,arg)
        select case (arg)
          case('-fnbim')
            call get_command_argument(i+1,fnbim)
	        inquire(file=fnbim,exist=exist)
            if (.not.exist) then
              write(*,'(/a,a,a)') ' Input *.bim file ',trim(fnbim),' not in directory'
              stop
            end if
          case('-fnout')
            call get_command_argument(i+1,fnout)
          case default
            write (*,'(a,a,/)') 'Unrecognized command-line option: ', arg
            stop
        end select
      end do
      
      write(*,'(/A//)') ' purge_SNPs v 1.1'
      
!     Open input plink *.bim file and determine number of coding variants (nv),
!     SNPs (ns), number of coding variants per SNP (nvs), and maximum number of
!     coding variants for any one SNP (nvsm).
      
      open (unit=15,file=fnbim)
      nv=0
      do
        read(15,*,iostat=ioerr) atmp
        if (ioerr<0) exit
        nv=nv+1
      end do
      
      rewind (15)
      ns=1
      read(15,*) atmp,atmp,atmp,dist1
      do i=2,nv
        read(15,*) atmp,atmp,atmp,dist2
        if (dist2/=dist1) then
          ns=ns+1
          dist1=dist2
        end if
      end do
      
      rewind (15)
      allocate(nvs(ns))
      ic1=1
      ic2=0
      read(15,*) atmp,atmp,atmp,dist1
      do i=2,nv
        read(15,*) atmp,atmp,atmp,dist2
        ic2=ic2+1
        if (dist2/=dist1) then
          nvs(ic1)=ic2
          ic1=ic1+1
          dist1=dist2
          ic2=0
        end if
      end do
      nvs(ns)=ic2+1
      nvsm=maxval(nvs)
      
      open(unit=17,file=fnout)
      
      rewind (15)
      allocate (snp(nvsm),a1(nvsm),a2(nvsm))
      do i=1,ns
        do j=1,nvs(i)
          read(15,*) atmp,snp(j),atmp,atmp,a1(j),a2(j)
        end do
        if (nvs(i)==1) then
          if (scan(a1(1),'MRWSYKVHDBNx')>0.or.scan(a2(1),'MRWSYKVHDBNx')>0) write(17,'(a)') trim(snp(1))
        else
          ic1=index(snp(1),'_',back=.true.)+1
          set_nt=''
          do j=1,nvs(i)
            set_nt=trim(set_nt)//snp(j)(ic1:len_trim(snp(j)))
          end do
          if (scan(set_nt,'MRWSYKVHDBNx')>0) then
            do j=1,nvs(i)
              write(17,'(a)') trim(snp(j))
            end do
          end if
        end if
      end do
      
      close (15)
      close (17)
            
          
      end program purge_SNPs
    