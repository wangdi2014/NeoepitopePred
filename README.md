# Neoepitope prediction (St Jude)

**Background**

Cancers are caused by somatically acquired alterations including single nucleotide variations (SNVs), small insertion/deletions (indels), translocations, and other types of rearrangements. The genes affected by these mutations may produce altered proteins, some of which may lead to the emergence of tumor-specific immunogenic epitopes. We developed an analytical workflow for identification of putative neoepitopes based on somatic missense mutations and gene fusions using whole genome sequencing data. The workflow was used to characterize the neoepitope landscape of 23 subtypes of pediatric cancer analyzed in the Pediatric Cancer Genome Project (PCGP).


**This repository contains the source code for the workflow of HLA typeing and neoepitope prediction. The pipelines are available for use on DNAnexus.**

Here you'll find the source code for two applets implementing individual pipeline stages of the neoepitope prediction workflow. 

  The stjude-hlatype applet is used for predict the HLA class I alleles. User can select to provide fastq (paired or single end reads) or BAM file as input. When using BAMs as input, the reads surrounding the HLA loci and unmapped reads will be extracted. The reads will be fed into Optitype for HLA typing. The default setting of Optitype is used. The output of the HLA type can be combined with the stjude-epitope app (see below) to perform affinity prediction of neoepitopes.
  
  The stjude-epitope applet is used to extract a peptide covering an array of tiling peptides (size defined by users) overlapping each missense mutation or gene fusion. Fusion junctions can be identified using RNAseq by CICERO (Li et al. unpublished data). NetMHCcons (Karosiene et al.) is subsequently used to predict affinities of the peptide array for each HLA receptor in each sample. The neoepitope with affinity lower than the threshold will be highlighted in output file (default 500 nM).



The shell and python scripts for building the applets and instantiating DNAnexus are included the root and bin directory. The applets can be used to construct a workflow [workflows](https://wiki.dnanexus.com/UI/Workflows) using them. You'll need the [DNAnexus SDK](https://wiki.dnanexus.com/Command-Line-Client/Quickstart) installed and set up to run these scripts.

**Reference**
T.-C. Chang, R. A. Carter, Y. Li, Y. Li, H. Wang, M. N. Edmonson, X. Chen, P. Arnold, T. L.  Geiger, G. Wu, J. Peng, M. Dyer, J. R. Downing, D. R. Green, P. G. Thomas, J. Zhang. The Neoepitope Landscape in Pediatric Cancers.






