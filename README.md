# Neoepitope prediction (St Jude)

**This repository contains the source code for the workflow of HLA typeing and neoepitope prediction. The pipelines are available for use on DNAnexus.**

Here you'll find the source code for two applets implementing individual pipeline stages, and python scripts in the root directory to build the applets and instantiate DNAnexus [workflows](https://wiki.dnanexus.com/UI/Workflows) using them. You'll need the [DNAnexus SDK](https://wiki.dnanexus.com/Command-Line-Client/Quickstart) installed and set up to run these scripts.

The stjude-hlatype app is used for predict the HLA class I alleles. User can select to provide fastq (paired or single end reads) or BAM file as input. When using BAMs as input, the reads surrounding the HLA loci and unmapped reads will be extracted. The reads will be fed into Optitype for HLA typing. The default setting of Optitype is used. The output of the HLA type can be combined with the stjude-epitope app (see below) to perform affinity prediction of neoepitopes.

The stjude-epitope app is used to extract a peptide covering an array of tiling peptides (size defined by users) overlapping each missense mutation or gene fusion. Fusion junctions can be identified using RNAseq by CICERO (Li et al. unpublished data). NetMHCcons (Karosiene et al.) is subsequently used to predict affinities of the peptide array for each HLA receptor in each sample. The neoepitope with affinity lower than the threshold will be highlighted in output file (default 500 nM).




