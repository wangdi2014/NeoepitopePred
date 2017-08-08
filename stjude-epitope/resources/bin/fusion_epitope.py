#!/usr/bin/python

import pandas as pd
import re
import os
import argparse

from Bio import Seq


def print_sample_number(df):
	print {_ind:len(_rows) for _ind, _rows in df.loc[:, ['sample', 'geneA', 'geneB']].drop_duplicates().groupby(['geneA', 'geneB']).groups.iteritems()}

def print_full(df, row=5, col=5):
	with pd.option_context('display.max_rows', row, 'display.max_columns', col):
		print(df)  

def fusion_peptide_extraction2(cicero_fusion_df, mrna_to_protein_df, psize=9):
	frag_size=psize*2-1
	df_list=[]
	clean_tmp()
	for (sample, geneA, geneB), group in cicero_fusion_df.groupby(['sample', 'geneA', 'geneB']):
                print "Processing {}".format(sample)
                if group.shape[0] > 1:
                        print "more than one entry for {}".format(sample)
			#print_full(group, row=4, col=100)
                ### questionable 
                ### identify the minimum of read counts of both sides for each enrty if more than one
                ### choose the entry with the maximum value of the minimum values of both sides
                col_min_series=group[['readsA', 'readsB']].apply(min,axis=1)
                #print "[info] column minimum"
                #print group[['readsA', 'readsB']]
                #print col_min_series

                top_rows_logical = col_min_series == col_min_series.max()
                new_group = group[top_rows_logical]
                print "[info] new group"
                print top_rows_logical
                #print_full(new_group, row=4, col=100)
		
		if top_rows_logical.sum() > 1:
                        print "more than one entry for {} tie for best row. Use occurence".format(sample)
                        #print new_group
                        new_group = new_group.iloc[[0]]
                #print_full(new_group, row=4, col=100)
		#print group.shape

		### my code to consider all potential fusion			
		for i in range(0,group.shape[0]):
			new_group=group.iloc[[i]]
			print new_group
		#### my code end
                	print new_group.index
                	print new_group.loc[new_group.index[0], 'sv_refseqA']
                	if (new_group.loc[new_group.index[0], 'sv_refseqA'] != new_group.loc[new_group.index[0], 'sv_refseqA']) or (new_group.loc[new_group.index[0], 'sv_refseqB'] != new_group.loc[new_group.index[0], 'sv_refseqB']):
                        	print "###\n###\n###\nNo reference accessions are given\n###\n###\n###"
                        	continue
			#print new_group.loc[new_group.index[0], 'sv_refseqA']
                	#print new_group.loc[new_group.index[0], 'sv_refseqB']
			### get protein seq for the mRNA from mrna2protein_df
			protA_seq = mrna2protein_df.loc[re.sub(",.+", "", new_group.loc[new_group.index[0], 'sv_refseqA']), 'protein']
                	protB_seq = mrna2protein_df.loc[re.sub(",.+", "", new_group.loc[new_group.index[0], 'sv_refseqB']), 'protein']
		
			#translate the fusion contig in 3 frames and acquire peptides
                	contig =  new_group.loc[new_group.index[0], 'contig']
			### peptides of frame 1,2,3 with flanking peptides
			### check frame
			#print contig
			in_frame_contig, pAset, pBset = inframe_check(contig, protA_seq, protB_seq, psize)
			geneA=new_group.loc[new_group.index[0], 'geneA']
			geneB=new_group.loc[new_group.index[0], 'geneB']
			peptide=peptide_src_check(in_frame_contig, pAset, pBset, geneA, geneB, psize)
			#print peptide
			if len(peptide) > 0:
				write_seq(sample, geneA, geneB, peptide)

def extra_samples(psize):
	pax3_fusion_16mer = 'IGNGLSPQSKFIRVQN'
	pax7_fusion_16mer = 'VSNGLSPQSKFIRVQN'
	ewsr_fli_fusion_31mer_7_to_6='PSQYSQQSSSYGQQNPSYDSVRRGAWGNN'
	ewsr_fli_fusion_31mer_7_to_8='PSQYSQQSSSYGQQNPYQILGPTSSRLAN'
	ewsr_fli_fusion_31mer_7_to_5='PSQYSQQSSSYGQQSSLLAYNTTSHTDQS'
	ewsr_fli_fusion_31mer_10_to_7='KPGGPMDEGPDLDLGPPLGGAQTISKNTE'
	#VIItoVI = ["SJEWS001301_D1", "SJEWS001302_D1", "SJEWS001303_D1", "SJEWS001304_D1", "SJEWS001305_D1", 'SJEWS001313_D1', 'SJEWS001316_D1', 'SJEWS001320_D1', 'SJEWS001318_D1', 'SJEWS001308_D1', 'SJEWS001311_D1','SJEWS001314_D1']
	#VIItoVIII = ["SJEWS001318_D1"]
	#VIItoV = ['SJEWS001308_D1', 'SJEWS001311_D1', 'SJEWS001312_D1', 'SJEWS001314_D1', 'SJEWS001317_D1']
	#XtoVII = ['SJEWS001306_D1','SJEWS001307_D1']
	### consider only samples with inframe ORF
	VIItoVI=["SJEWS001305_D1", 'SJEWS001320_D1']
	VIItoVIII = ["SJEWS001318_D1"]
	VIItoV = ['SJEWS001308_D1', 'SJEWS001311_D1', 'SJEWS001314_D1']
	XtoVII = ['SJEWS001306_D1']

	if psize:
		ewsr_fli_fusion_31mer_7_to_6=ewsr_fli_fusion_31mer_7_to_6[13-psize+1:13] + ewsr_fli_fusion_31mer_7_to_6[14] + ewsr_fli_fusion_31mer_7_to_6[15: 15 + psize-1]
		ewsr_fli_fusion_31mer_7_to_8=ewsr_fli_fusion_31mer_7_to_8[13-psize+1:13] + ewsr_fli_fusion_31mer_7_to_8[14] + ewsr_fli_fusion_31mer_7_to_8[15: 15 + psize-1]
		ewsr_fli_fusion_31mer_7_to_5=ewsr_fli_fusion_31mer_7_to_5[13-psize+1:13] + ewsr_fli_fusion_31mer_7_to_5[14] + ewsr_fli_fusion_31mer_7_to_5[15: 15 + psize-1]
		ewsr_fli_fusion_31mer_10_to_7=ewsr_fli_fusion_31mer_10_to_7[13-psize+1:13] + ewsr_fli_fusion_31mer_10_to_7[14] + ewsr_fli_fusion_31mer_10_to_7[15: 15 + psize-1]
	write_seq('SJRHB007_D', 'PAX3', 'FOXO1', pax3_fusion_16mer)
	write_seq('SJRHB008_D', 'PAX3', 'FOXO1', pax3_fusion_16mer)
	write_seq('SJRHB009_D', 'PAX7', 'FOXO1', pax7_fusion_16mer)
	write_seq('SJRHB010_D', 'PAX7', 'FOXO1', pax7_fusion_16mer)
	for sample in VIItoVI:
		write_seq(sample, 'EWSR1', 'FLI1', ewsr_fli_fusion_31mer_7_to_6)
	for sample in VIItoVIII:	
		write_seq(sample, 'EWSR1', 'FLI1', ewsr_fli_fusion_31mer_7_to_8)
	for sample in VIItoV:
		write_seq(sample, 'EWSR1', 'FLI1', ewsr_fli_fusion_31mer_7_to_5)
	for sample in XtoVII:
		write_seq(sample, 'EWSR1', 'FLI1', ewsr_fli_fusion_31mer_10_to_7)

def write_seq(sample,gene1,gene2, sequence):
	fn=sample+'.fusion_flanking.seq'
	f=open(fn,'a')
	f.write('>'+gene1+'_'+gene2+'\n'+sequence+'\n')
	f.close()

def peptide_src_check(inframe_ctg, proteinA_set, proteinB_set, gene_A, gene_B, psize):
	#Determine the source of the peptides in the fusion contig.
        #Then, extract the peptides that are novel and that are found between 
        #the A-matching peptides and the B-matching peptides
        previous_encounter = None
        encountered_a_logical = None
        in_a_logical = None
        encountered_b_logical = None
        in_b_logical = None

        fusion_contig_peptide_source = []
	contig_peptide_seq = []
	junction_pos=[]
	final_junction_pos=[]
        for _ind in range(0, len(inframe_ctg) - psize + 1):
        	loc= ''
                pep_query = inframe_ctg[_ind:(_ind + psize)]
                if pep_query in proteinA_set:
			if encountered_a_logical and (len(junction_pos) > 0) and gene_A == gene_B:
                        	final_junction_pos = junction_pos
			encountered_a_logical=True
			in_a_logical = True
                        in_b_logical = False
                        loc = 'A'
			junction_pos=[]
                elif pep_query in proteinB_set:
			encountered_b_logical=True
			in_a_logical = False
                        in_b_logical = True
                        loc = 'B'
			if (previous_encounter == 'A') and (len(junction_pos) > 0) and gene_A != gene_B:
                         	final_junction_pos = junction_pos
                                junction_pos=[]
                else:
			if in_a_logical:
                        	previous_encounter = 'A'
                        elif in_b_logical:
                        	previous_encounter = 'B'
                        in_a_logical = False
                        in_b_logical = False
			loc = '-'
			junction_pos.append(_ind)
                fusion_contig_peptide_source.extend(loc)
		contig_peptide_seq.extend(inframe_ctg[_ind-1:_ind])
	#print "Sample: {}, variant: {}".format(sample, geneA + "_" + geneB)
	print fusion_contig_peptide_source
	print contig_peptide_seq
	print final_junction_pos
	target_peptide=''
	if final_junction_pos:
		target_peptide=inframe_ctg[min(final_junction_pos):max(final_junction_pos)+psize]
	return target_peptide
        #print "List of junction peptides: {}".format(" ".join(final_junction_nonamer_list))
        #print contig
        #print "#####"

def clean_tmp(dir='.', ext='.fusion_flanking.seq'):
        for File in os.listdir(dir):
                if File.endswith(ext):
                        print("[warning] Remove all temp files from the current folder")
			os.remove(os.path.join(dir, File))

def inframe_check(seq_for_check, proteinA_seq, proteinB_seq, psize):
	protA_peptides_set = set([proteinA_seq[_ind:(_ind + psize)] for _ind in range(0, len(proteinA_seq) - psize + 1)])
	protB_peptides_set = set([proteinB_seq[_ind:(_ind + psize)] for _ind in range(0, len(proteinB_seq) - psize + 1)])	
	all_peptides_set = protA_peptides_set | protB_peptides_set		
	### peptides of frame 1,2,3 with flanking peptides	
	f1_contig = str(Seq.Seq(seq_for_check).translate())
	f1_contig_peptides_set = set([f1_contig[_ind:(_ind + psize)] for _ind in range(0, len(f1_contig) - psize + 1)])
	f2_contig = str(Seq.Seq(seq_for_check[1:]).translate())
        f2_contig_peptides_set = set([f2_contig[_ind:(_ind + psize)] for _ind in range(0, len(f2_contig) - psize + 1)])
	f3_contig = str(Seq.Seq(seq_for_check[2:]).translate())
        f3_contig_peptides_set = set([f3_contig[_ind:(_ind + psize)] for _ind in range(0, len(f3_contig) - psize + 1)])
	
	#determine which frame is correct by removing the non-fusion gene peptides from the gene fusion peptides
        #and ensuring there is a maximum of (peptide_length - 1) fusion-specific peptides
        f1_unique_peptides_set = f1_contig_peptides_set - all_peptides_set
        #print "frame1 : {}".format(len(f1_unique_peptides_set))
        f2_unique_peptides_set = f2_contig_peptides_set - all_peptides_set
        #print "frame2 : {}".format(len(f2_unique_peptides_set))
        f3_unique_peptides_set = f3_contig_peptides_set - all_peptides_set
        #print "frame3 : {}".format(len(f3_unique_peptides_set))
        min_novel_contigs = min(len(f1_unique_peptides_set), len(f2_unique_peptides_set), len(f3_unique_peptides_set))

	if(len(f1_unique_peptides_set) == min_novel_contigs):
		print "[info] detected as frame 1"
        	inframe_contig = f1_contig
        elif(len(f2_unique_peptides_set) == min_novel_contigs):
		print "[info] detected as frame 2"
        	inframe_contig = f2_contig
        else:
		print "[info] detected as frame 3"
                inframe_contig = f3_contig
	return inframe_contig, protA_peptides_set, protB_peptides_set

if __name__=='__main__':

	parser=argparse.ArgumentParser(description='Extract peptides for netMHCcons')
	parser.add_argument('-s', '--size_of_peptide', help='Peptide size for prediction of affinity', default=9, required=True)
	parser.add_argument('-f', '--fusion_file', help='Peptide size for prediction of affinity', required=True)
	args=parser.parse_args()

	### filenames 
	#fusion_f1='/nfs_exports/genomes/1/projects/HLA/PCGP/fusion_proteins/fusion_prots_final_minus_three_fusions_03012016.tsv'
	#fusion_f2='/nfs_exports/genomes/1/projects/HLA/GENE_FUSION/all_three_fusions_catted.txt'
	#fusion_f3='/rgs01/project_space/zhanggrp/PCGP_viral/cmpb/HLA/MyPipeline/fusion/CICERO.Epitope.Nov01.txt'
	fusion_f=str(args.fusion_file)
	mrna2protein_f='/usr/bin/fusion_data/refgene2protein.tab'
	#sample_ls='/nfs_exports/genomes/1/projects/HLA/PCGP/diagnostic_sample_ids.txt'

	### print parameters
	print("Peptide size: %s" % str(args.size_of_peptide))
	print("Used fusion files: %s" % str(fusion_f))
	print("mRNA to protein conversion file: %s" % mrna2protein_f)
	#print("Sample list: %s" % sample_ls)
	

	### load fusion data frame from Cicero
	fusion_df=pd.read_csv(fusion_f, header=0, sep="\t")
	print("[Info] fusion dataset")
	print_sample_number(fusion_df)


	### load mrna to protein translation from Edmanson 
	mrna2protein_df=pd.read_csv(mrna2protein_f, header=0, sep='\t', index_col=0)
	print("\n[info] number of mRNA that can be translated into protein: %s" % len(mrna2protein_df))



	### merge fusion dataframes and filter the samples in the interested list
	#fusion_df_ls=fusion_df[fusion_df['sample'].isin(sample_ids['sample'])]
	fusion_df_ls=fusion_df


	### extract peptide with expected size

	#fusion_peptide=fusion_peptide_extraction(fusion_df_ls, mrna2protein_df, psize=size)
	#print_full(fusion_peptide)
	print("[Run] start extracting peptides with size of %s" % str(args.size_of_peptide))
	fusion_peptide_extraction2(fusion_df_ls, mrna2protein_df, psize=int(args.size_of_peptide))
	
	### made custom peptide sequences for eatra fusions
	#extra_samples(psize=int(args.size_of_peptide))

