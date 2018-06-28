# Neoepitope prediction (St Jude)

**Background**

Cancers are caused by somatically acquired alterations including single nucleotide variations (SNVs), small insertion/deletions (indels), translocations, and other types of rearrangements. The genes affected by these mutations may produce altered proteins, some of which may lead to the emergence of tumor-specific immunogenic epitopes. We developed an analytical workflow for identification of putative neoepitopes based on somatic missense mutations and gene fusions using whole genome sequencing data. The workflow has been used to characterize neoepitope landscape of 23 subtypes of pediatric cancer in the Pediatric Cancer Genome Project (PCGP)[1].


**Repository**

Here you'll find the source code for the applets implementing individual pipeline stages of the neoepitope prediction workflow. 

  The stjude-hlatype applet is used for predicting the HLA class I alleles. User can select to provide fastq (paired or single end reads) or a BAM file as input. When using BAMs as input, the reads surrounding the HLA loci and unmapped reads will be extracted. The reads will be fed into Optitype for HLA typing. The default setting of Optitype is used. The output of the HLA type can be combined with the stjude-epitope app (see below) to perform affinity prediction of neoepitopes.
  
  The stjude-epitope applet is used to extract peptides covering an array of tiling peptides (size defined by users) overlapping each missense mutation or gene fusion. Fusion junctions can be identified using RNAseq by CICERO (Li et al. unpublished data). NetMHCcons (Karosiene et al.) is subsequently used to predict affinities of the peptide array for each HLA receptor in each sample. The neoepitope with affinity lower than the threshold will be highlighted in the output file (default 500 nM).

  The shell and python scripts for building the applets and instantiating DNAnexus are included the root and bin directory. The applets can be used to construct DNAnexus [workflows](https://wiki.dnanexus.com/UI/Workflows). You'll need the [DNAnexus SDK](https://wiki.dnanexus.com/Command-Line-Client/Quickstart) installed and set up to run these scripts.

**Prerequisite**  
The workflow incoporated the following softwares:  
OptiType 1.0 [2]  
netMHCcons 1.1 [3]

The stjude-epitope applet requires human reference genome (Hg19) which is not included. Users can obtain the genome sequence from [UCSC](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/) and place it into the following folder: stjude-epitope/resources/usr/genome


**References**  
[1] Chang T-C, Carter R, Li YJ, Li YX, Wang H, Edmonson M, Chen X, Arnold P, Geiger T, Wu G, Peng JM, Dyer M, Downing J, Green D, Thomas P, Zhang JH: The Neoepitope Landscape in Pediatric Cancers.  
[2] Szolek A, Schubert B, Mohr C, Sturm M, Feldhahn M, Kohlbacher O: OptiType: precision HLA typing from next-generation sequencing data. Bioinformatics 2014, 30:3310-3316.  
[3] Karosiene E, Lundegaard C, Lund O, Nielsen M: NetMHCcons: a consensus method for the major histocompatibility complex class I predictions. Immunogenetics 2012, 64:177-186.





