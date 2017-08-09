#!/usr/bin/python
from __future__ import print_function ### for use 'end=' in print function
__author__='TCC_2017'
from MutPepExpression import data_io
from MutPepExpression.args_handling import parse_args
from MutPepExpression import netmhccons_formatter
from MutPepExpression.args_handling import Utils
from StringIO import StringIO
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from pyliftover import LiftOver

import sys
import argparse
import numpy as np
import pandas as pd
import re
import os
import commands
import os.path
import pysam
import random
import copy
import string
import logging
import tempfile

#### Functions for netMHCcons

def f2ary(fn, head=None, reannotate=True, version_check=True, mRNA_acc_check=True, missense=True):
	file=pd.read_csv(fn, header=int(head)-1,sep='\t')
	file=file.drop_duplicates()
	#print(file)
	if head is None:
		if file.shape[1] != 16:
			print ("[Error] the mutation data frame not in correct format ... quit")
			quit()
		else:
			print ("[Note] header was added to the mutation dataframe")
			file.columns=['GeneName', 'SJQuality', 'Sample', 'Chr', 'Position_hg19', 'Class', 'AAChange', 'ProteinGI', 'mRNA_acc', '#Mutant_In_Tumor', '#Total_In_Tumor', '#Mutant_In_Normal', '#Total_In_Normal', 'ReferenceAllele', 'MutantAllele', 'Flanking']
			print ("[Note] header was added as: %s" % list(file))
	if version_check:
		print("Checking mutation dataframe version")
		file=PCGP_mut_df_genome_build_check(file,pos_col=4)
		print(file)
	if reannotate:
		print("Reannotation")
		file=mut_annotate(file)
	if missense:
		print("Keep only missense SNV")
		file=file[file['Class'].str.contains('missense', case=False)]
	if mRNA_acc_check:
		print("Drop mutation without mRNA accessions")
		#print(list(file))
		#file=file[ pd.notnull(file['mRNA_acc'])]
		file=file[file['mRNA_acc'].str.contains('_')]
		#file=file[file['annovar_acc'].str.contains('_')]	
		#print("Reannotating the mutations without mRNA accessions")
		#file=file[ pd.isnull(file['mRNA_acc'])]
		#mut_annotate(file)
		#mRNA_acc(file)

	#print(list(file))
	#print(file)
	if(file.shape[0]==0):
		print("No mutations qualified for continue")
		quit()
	else:
		return file

def mRNA_acc(df):
	database='gene'
	Entrez.email='ti-cheng.chang@stjude.org'
	batchSize=100
	retmax=10**9
	for index, row in df.iterrows():
		gn=row['GeneName']
		### Esearch
		esearch_handle=Entrez.esearch(db=database, term=gn+' AND srcdb_refseq[PROP]')
		esearch_result=Entrez.read(esearch_handle)
		### there is a bug in the biopython parser.py, that may raise the error of "local variable 'url' referenced before assignment". Need to modify it based on solution on web:
		### https://github.com/biopython/biopython/commit/0c4bd6aec730c78dd11c4d7d79a68eda45125beb
		esearch_handle.close()
		#print (esearch_result)
		mRNA_gi_number=esearch_result["IdList"]
		print ("Gene: %s, mRNA: %s" % (gn, mRNA_gi_number))
	'''		
		### Elink
			retrieved_mRNA_uids = []
			efetch_handle=Entrez.elink(dbfrom="nucleotide", db='nuccore', id=mRNA_gi_number)
			efetch_result=Entrez.read(efetch_handle)
			efetch_handle.close()
		print( efetch_result)

			for each_record in elink_result:
					#print(acc, each_record)
				 		try:
							protein_id=each_record["LinkSetDb"][0]["Link"][0]["Id"]
							retrieved_mRNA_uids.append(protein_id)
					except IndexError:
							print ("[Error] no protein link: %s" % (acc))
		print (retrieved_mRNA_uids)
	'''
	quit()
		

	
def np_print_ary(dt):
	print(np.matrix(dt))

def PCGP_mut_df_genome_build_check(df,pos_col=4):
	col_check_hg18= [ col for col in df.columns if 'hg18' in col.lower() ]
        col_check_hg38= [ col for col in df.columns if 'hg38' in col.lower() ]
	if len(col_check_hg18) > 0 or len(col_check_hg38) > 0:
		if (len(col_check_hg18) == 1 and len(col_check_hg38) == 0) or (len(col_check_hg18) == 0 and len(col_check_hg38) == 1):
                        if len(col_check_hg18) == 1:
			    fd=col_check_hg18[0]
                            col_check=col_check_hg18
			    print("[Warning] following columns from hg18 genome build: %s" % fd)
			    lo=LiftOver('hg18', 'hg19')
                        elif len(col_check_hg38) == 1:
                            fd=col_check_hg38[0]
                            col_check=col_check_hg38
                            print("[Warning] following columns from hg38 genome build: %s" % fd)
                            lo=LiftOver('hg38', 'hg19')
			pos=[]
			#print(df)
                        print(fd)
			for idx, row in df.iterrows():
				conversion=lo.convert_coordinate(row['Chr'], row[col_check[0]])
				if conversion:
					newpos=lo.convert_coordinate(row['Chr'], row[col_check[0]])[0]
					pos.append(newpos[1])
				else:
					newpos=(row['Chr'],-1)
					pos.append(0)
				#newpos=lo.convert_coordinate(row['Chr'], row[col_check[0]])[0]
				#pos.append(newpos[1])
			
			df['Position_hg19']=pos
			return df	
				
		else:
			print("[Error] only one column allowed for conversion: %s ... quit" % col_check)
			quit()
	else:
		#print("No change")
		cols=df.columns.values
		cols[pos_col]='Position_hg19'
		df.columns=cols
		return df
def get_genes(df, psize=9):
	print("number of rows: %d\nnumberof columns: %d\nHeader: %s" % (df.shape[0], df.shape[1], list(df)))
	### method 1: loop over all rows and get sequences
	cnt=0
	fn_ls=[]
	for index, row in df.iterrows():
		### sample0, variant2, chr3, pos4, type5, refa6, alta7, mrnaacc10
		cnt=cnt+1
		
		#acc=row['annovar_acc']
		#mut=row['annovar_peptide_change']
		acc=row['mRNA_acc']
		mut=row['AAChange']
		sam=row['Sample']
		chr=row['Chr']
		pos=row['Position_hg19']
		#genename=row['annovar_gene']
		genename=row['GeneName']
		#mut_type=row['class']
		mut_type=row['Class']
		#annotation=str(row['annotation']).split(',')	
		### check if mutation format correct	
		pattern=re.compile("[A-Z]\d+[A-Z]")
		if pattern.match(mut):
			print("[Mut format validated] %s" % mut)
		else:
			print("[Mut format invalidated] %s" % mut)
			continue

		if str(mut_type.lower()) == 'missense':
			if '_' not in mut:
				mut = genename+'_'+mut
			if not pd.isnull(acc):
				### use entrez to get sequnce on the fly
				#problematic_acc=['NM_001146344', 'NM_001002252', 'NM_006108']
				#if acc in problematic_acc:
				#	mutseq1=entrez(acc, mut, psize)
				#	mutseq2=entrez_alternative(acc, mut, chr, psize)
				#	print("[Warn mutseq] (web) %s; (local) %s" % ( mutseq1, mutseq2))
				#	mut_seq = entrez(acc, mut, psize)
					
				### parse refflat file to get sequnce locally
				#else:
					#print("%s %s %s %s" % (acc, mut, chr, psize))
				mut_seq = entrez_alternative(acc, mut, chr, pos, psize)
			else:
				print("No accession now")
			if mut_seq and len(mut_seq) > 0:
				print(mut_seq)
				fn=sam+'.snv_flanking.seq'
				if cnt == 1:
					f=open(fn, 'w')
				else:
					f=open(fn, 'a') 		
				#f.write('>'+str(df[i][sam_idx])+'.'+str('.'.join(df[i][2:5])+'\n'+mut_seq+'\n'))
				#f.write('>'+mut+'.'+str(df[i][sam_idx])+'.'+str(df[i][chr_idx])+'.'+str(df[i][pos_idx])+'\n'+str(mut_seq)+'\n')
				f.write('>'+mut+'\n'+str(mut_seq)+'\n')
				f.close()
				if not fn in fn_ls:
					fn_ls.append(fn)

	if 'mRNA_acc' not in df:
		df['mRNA_acc']=df['annovar_acc']
			
	#print(df)	
	return fn_ls, df



def mut_annotate(df, annovar='/nfs_exports/apps/gnu-apps/NextGen/gwbin/annotate_variation.pl', build='hg19', db='/nfs_exports/genomes/1/Homo_sapiens/Hg19/ANNOVAR/'):
	temp_name = next(tempfile._get_candidate_names()) + '.temp'
	f=open(temp_name,'w')
	for index, row in df.iterrows():
		chr=row['Chr']
		pos=row['Position_hg19']
		ref='0'
		alt='0'
		if pd.notnull(row['ReferenceAllele']):
			ref=row['ReferenceAllele']
		if pd.notnull(row['MutantAllele']):
			alt=row['MutantAllele']
		content='%s\t%s\t%s\t%s\t%s\n' % (chr, pos, pos, ref, alt)
		f.write(content)
	f.close()
	cmd='%s -geneanno %s -dbtype refGene -neargene 1000 -buildver %s %s' % (annovar, temp_name, build, db)
	print("[cmd annovar] " + cmd)
	status, output = commands.getstatusoutput(cmd)
	ann=pd.read_csv(temp_name+'.variant_function', sep="\t", header=None)
	ann2=pd.read_csv(temp_name+'.exonic_variant_function', sep="\t", header=None)
	os.remove(temp_name)
	os.remove(temp_name+'.variant_function')
	os.remove(temp_name+'.exonic_variant_function')
	os.remove(temp_name+'.log')
	ann2.columns=['line', 'class', 'annotation', 'chr', 'pos','pos2', 'ref allele', 'alt allele'] 
	ann2 = ann2[['class', 'annotation', 'chr', 'pos']]

	df=pd.merge(df, ann2, left_on=['Chr', 'Position_hg19'], right_on=['chr', 'pos'], how='left')
	for idx, row in df.iterrows():
		#print(row['annotation'])
		annotation=str(row['annotation']).split(',')
		Annovar_gene=''
		Annovar_acc=''
		Annovar_pvar=''
		for item in annotation:
			if len(item) > 0:
				items=item.split(':')
				#print("[Items] %s" % items)
				if Annovar_gene=='' and items[0].lower() !='nan' and items[0].lower() !='unknown':
					Annovar_gene=items[0]
					Annovar_acc=items[1]
					Annovar_pvar=items[4]
				if items[0] == row['GeneName'] and 'mRNA_acc' in row and items[1] == row['mRNA_acc']:
					Annovar_gene=items[0]
					Annovar_acc=items[1]
					Annovar_pvar=items[4]	
		### recheck if any entries cannot be reannotate
		#print("[New annotation] %s %s %s" %( Annovar_gene, Annovar_acc, Annovar_pvar))
		if Annovar_gene=='' and row['Class'].lower()=='missense':
			if len(row['AAChange']) > 0 and 'mRNA_acc' in row:
				print ("[Warn] Reannotation failed for %s, use original annotation instead" % row['GeneName'])
				Annovar_gene=row['GeneName']
				Annovar_acc=row['mRNA_acc']
				Annovar_pvar=row['AAChange']
				df.ix[idx,'class']=row['Class']
			else:
				print ("[Warn] No annotation available for %s" % row['GeneName'])		
		#print("[New annotation] %s %s %s" %( Annovar_gene, Annovar_acc, Annovar_pvar))		
		#print("gene: %s, var: %s" % (Annovar_gene, Annovar_acc))
		df.ix[idx,'annovar_gene']=Annovar_gene
		df.ix[idx,'annovar_acc']=Annovar_acc
		Annovar_pvar=Annovar_pvar.lstrip('p.')
		df.ix[idx,'annovar_peptide_change']=Annovar_pvar
	df['class']=df['class'].replace(['nonsynonymous SNV'], 'missense')
	return df

def entrez_alternative (acc, mut, chr, pos, psize=9, bam_build='hg19', check_f=True):
	bam_build_to_ref_locations_dict = {
		'hg18':'/home/dnanexus/genome/hg18.fa',
		'hg19':'/home/dnanexus/genome/hg19.fa',
                #'hg19':'hg19.fa',
		'GRCh37-lite':'/home/dnanexus/genome/GRCh37-lite.fa',
                'HG19_Broad_variant':'/home/dnanexus/genome/Homo_sapiens_assembly19.fasta',
                'g1k-human-build37':'/home/dnanexus/genome/GRCh37-lite.fa',
                'NCBI36_WUGSC_variant':'/home/dnanexus/genome/NCBI36_WUGSC_variant.fa'
	}
	refflat_filename = data_io.get_refflat_filename(bam_build)
	reference_fa_filename = data_io.get_ref_fa_filename(bam_build, bam_build_to_ref_locations_dict)
	fastaobj = pysam.Fastafile(reference_fa_filename)
	refflat_df = data_io.read_refflat_as_df(refflat_filename, bam_build)

	### make sure the comptability of the dataframe for "get_df_of_peptide_and_genomic_rna_coverage_for_sample 
	#sample_muts_df = data_io.ensure_refflat_muts_df_accession_compatibility(muts_df, refflat_df)
	compatible_refflat_df = data_io.get_build_compatibile_refflat(refflat_df, bam_build)
	#compatible_sample_muts_df = data_io.get_build_compatibile_muts_df(sample_muts_df, bam_build)
	#print(list(compatible_refflat_df))
	#print("Finding positions for %s" % acc)
	
	
	if acc in compatible_refflat_df.index:
		acc_sel=compatible_refflat_df.loc[acc,:]
		if isinstance(acc_sel, pd.DataFrame):
			acc_sel=acc_sel.loc[acc_sel['chrom']==str(chr)]
			acc_sel=acc_sel.loc[acc,:]
		if isinstance(acc_sel, pd.DataFrame):
			print("%s, %s" % (pos, type(acc_sel)))
			print(acc_sel)
			acc_sel=acc_sel.loc[(acc_sel['cdsStart'] < pos) & (pos < acc_sel['cdsEnd'])]
			if acc_sel.shape[0] > 0:
				acc_sel=acc_sel.loc[acc,:]
			else:
				print ("[Error] refflat file return nothing for %s, %s, after filtering" %(acc, mut))
				return None
		if isinstance(acc_sel, pd.DataFrame):
			print ("[Error] multiple mRNA in the refflat files")
			return None
			#quit()
		#print(type(acc_sel))
		#print(acc_sel)
		chr=acc_sel['chrom']
		#print(acc_sel['exonStarts'])
		eS=acc_sel['exonStarts'].split(',')
		eE=acc_sel['exonEnds'].split(',')
		cS=int(acc_sel['cdsStart'])
		cE=int(acc_sel['cdsEnd'])
		exon=zip(eS,eE)
		cds=""
		#print("\tstrand: %s" % acc_sel['strand'])
		#print(cS)
		#print(cE)
		#print(exon)
		for e in exon:
			start=int(e[0])
			end=int(e[1])
			seq=""
			if start <= cS and end > cS:
				if end <= cE:
					#print ("cds start: exon: %d, cds: %d, %d" %(start, cS, end))
					seq=fastaobj.fetch(chr, cS-1, end)
					cds=cds+seq
				else:
					#print ("single exon CDS: cds: %d, cds:%d" %(cS, cE))
					seq=fastaobj.fetch(chr, cS-1, cE)
					cds=cds+seq
			elif start > cS	and end < cE:
				#print ("internal exon: %d, %d" %(start, end))
				seq=fastaobj.fetch(chr, start-1, end)
				cds=cds+seq
			elif start < cE	and end >= cE:
				#print ("terminal exon: %d, cds: %d, exon:%d" %(start, cE, end))
				seq=fastaobj.fetch(chr, start-1, cE)
				cds=cds+seq
			#else:
				#print ("UTR: %d, %d" %(start, end))
			#print(seq)
		if acc_sel['strand'] == '-':
			cds=reversecomplement_dna_from_string(cds)
		#print(cds)
		pt=translate_dna_from_string(cds)
		pt=pt.rstrip('\*.*$')
		#print(pt)
		gene_mut=re.compile("(.*)_(.)(\d+)(.)").split(mut.lstrip())
		gene=gene_mut[1]
		ref=gene_mut[2]
		pos=int(gene_mut[3])
		alt=gene_mut[4]
		pt=SeqRecord(Seq(pt,IUPAC.protein), id=mut, name=mut, description=mut)
		#problematic_acc = ['NM_001146344', 'NM_001002252']
		#if acc in problematic_acc:
		#	mutated_seq=seq_windows(pt, pos, ref, alt, n=psize, check=False)
		#else:
		mutated_seq=seq_windows(pt, pos, ref, alt, n=psize, check=check_f)
		start=pos-psize
		if start < 0:
			start=0
		if not mutated_seq:
			print ("[Error] the seq for %s is null" % acc)
			return None
		else:
			mutated_pep_seq=mutated_seq[start:pos+psize-1]
			return mutated_pep_seq
	else:
		return None


def entrez(acc,mut, psize=9, check_f=True):
	database='nuccore'
	Entrez.email='ti-cheng.chang@stjude.org'
	batchSize=100
	retmax=10**9

	### Esearch
	esearch_handle=Entrez.esearch(db=database, term=acc)
	esearch_result=Entrez.read(esearch_handle)
	### there is a bug in the biopython parser.py, that may raise the error of "local variable 'url' referenced before assignment". Need to modify it based on solution on web:
	### https://github.com/biopython/biopython/commit/0c4bd6aec730c78dd11c4d7d79a68eda45125beb
	esearch_handle.close()
	#print (esearch_result)
	#print (esearch_result["IdList"][0])
	mRNA_gi_number=esearch_result["IdList"][0]
	
	### Elink
	retrieved_protein_uids = []
	elink_handle=Entrez.elink(dbfrom="nuccore", db="protein", LinkName="nuccore_protein", id=mRNA_gi_number)
	elink_result=Entrez.read(elink_handle)
	elink_handle.close()
	for each_record in elink_result:
		#print(acc, each_record)
		try:
			protein_id=each_record["LinkSetDb"][0]["Link"][0]["Id"]
			retrieved_protein_uids.append(protein_id)
		except IndexError:
			print ("[Error] no protein link: %s" % (acc))
	#print (retrieved_protein_uids)
	
	### Epost
	#epost_handle = Entrez.epost(db="protein", id=",".join(retrieved_protein_uids))
	#epost_result = Entrez.read(epost_handle)
	#epost_handle.close()
	
	#webenv=epost_result["WebEnv"]
	#query_key= epost_result["QueryKey"]
	
	### Efetch
	count=len(retrieved_protein_uids)
	batch_size=20
	the_records=""
	mutated_pep_seq=""
	#for start in range(0,count, batch_size):
	if count==1:
		
		### multiple records:
		#end=min(count, start + batch_size)
		#print("Fetching records %i thru %i..." % (start + 1, end))
		#fetch_handle = Entrez.efetch(db="protein", rettype="fasta", retmode="text", retstart=start, retmax=batch_size, webenv=webenv, query_key=query_key)
		
		### single records:
		fetch_handle = Entrez.efetch(db="protein", rettype="fasta", retmode="text",id=retrieved_protein_uids)
		#seq=fetch_handle.read()
		for seq in SeqIO.parse(fetch_handle, "fasta"):
			print ("%s\t%s\t%s\t%s" %(acc, mut, seq.id, seq.description))
		gene_mut=re.compile("(.*)_(.)(\d+)(.)").split(mut.lstrip())
		gene=gene_mut[1]
		ref=gene_mut[2]
		pos=int(gene_mut[3])
		alt=gene_mut[4]
		#mutated_seq=seq_windows(seq, pos, ref, alt, n=psize)
		problematic_acc = ['NM_001146344', 'NM_001002252']
		if acc in problematic_acc:
			mutated_seq=seq_windows(seq, pos, ref, alt, n=psize, check=False)
		else:
			mutated_seq=seq_windows(seq, pos, ref, alt, n=psize, check=check_f)

		start=pos-psize
		if start < 0:
			start=0
		if not mutated_seq:
			print ("[Error] the seq for %s is null" % acc)
			return None
		else:
			mutated_pep_seq=mutated_seq[start:pos+psize-1]
			return mutated_pep_seq

		#print("%s" % seq.seq)
		#print("\t%s" % seq.seq[pos-psize:pos+psize-1])
		#print("\t%s" % mutated_pep_seq)
		fetch_handle.close()
		#the_records = the_records + seq
	elif count==0:
		print ("[Error] no protein match of the mRNA acc: %s" % (acc))
	else:
		print ("[Error] 1-to-many match of mRNA acc: %s\t%s" % (acc, ";".join(retrieved_protein_uids)))
	#print (the_records)
	return mutated_pep_seq

def seq_windows(seq_obj,mut_pos,ref, alt, n=9, check=True):
	ID=seq_obj.id
	refseq=seq_obj.seq
	target=refseq[mut_pos-1:mut_pos]
	if check:
		if target.strip()[0] == ref.strip():	
			mutseq=refseq[0:mut_pos-1]+alt+refseq[mut_pos:]
			return mutseq
			#for pos in range(mut_pos-n, mut_pos):
			#	subseq=mutseq[pos:pos+n]
			#	print("\t%s" % (subseq))
		else:
			print("[Error] the target residue is wrong: ID..%s\ttarget..%s\tref..%s" %(ID, target, ref))
	else:
		print("[Warn] no check: ID..%s\ttarget..%s\tref..%s" %(ID, target, ref))
		mutseq=refseq[0:mut_pos-1]+alt+refseq[mut_pos:]
		return mutseq
def pd_print_full(x):
	pd.set_option('display.max_rows', len(x))
	print(x)
	pd.reset_option('display.max_rows')

	
def netMHCcons(peptide_fn, allele="HLA-A02:01", MHC="netMHCcons", filter=0, size=9, nM=500):
	if os.path.isfile(MHC):
		out=peptide_fn + '.affinity.xls'
		#cmd=MHC + " -inptype 0 -f " + peptide_fn + ' -a ' + allele + ' -length ' + str(size) + '| grep -E -v \'^#|^$|^---|pos|Number\''
		### with default netMHCcons Excel output
		#cmd=MHC + " -inptype 0 -filter " + str(filter) + " -f " + peptide_fn + ' -a ' + str(allele) + ' -length ' + str(size) + ' -xls 1 -xlsfile ' + out + ' | grep -E -v \'^#|^$|^---|^\s+pos|Number\''
		### without default netMGCcons Excel output and never filter at the begining
		cmd=MHC + " -inptype 0 -filter " + str(0) + " -f " + peptide_fn + ' -a ' + str(allele) + ' -length ' + str(size) + ' | grep -E -v \'^#|^$|^---|^\s+pos|Number\''
		#cmd=MHC + " -inptype 0 -filter " + str(0) + " -f " + peptide_fn + ' -a ' + str(allele) + ' -length ' + str(size)		

		print("[cmd] %s" % cmd)
		status, output = commands.getstatusoutput(cmd)
		#print(output)
		MHC_out=pd.read_table(StringIO(output), header=None, sep='\s+', names=['Position','Allele', 'Peptide','ID', '1-log50k', 'nM', 'Rank', 'Binding'])
		#print(MHC_out.values)
		MHC_out_filtered=MHC_out
		if filter == 1:
			MHC_out_filtered = MHC_out[MHC_out['nM'] <= int(nM)]	
		#pd_print_full(MHC_out_filtered)
		return MHC_out_filtered
	else:
		print("[Error] netMHCcons not found ... quit\n")
		quit()

def clean_tmp(dir='.', ext='.snv_flanking.seq'):
	for File in os.listdir(dir):
		if File.endswith(ext):
			print("[warning] Remove all temp files from the current folder")
			os.remove(os.path.join(dir, File))
		
def all_netMHCcons(snv, out_fn, dir='.', ext='.fornetmhccons.tmp', Peptide_size=9, HLA_allele="HLA-A02:01", Report_filter=0, MHC_binary="netMHCcons", nM_cutoff=500):
	tmp_file_cnt=0
	for fas in os.listdir(dir):
		if fas.endswith(ext):
			print("[Processing] run netMHCcons for %s" % fas)
			out=netMHCcons(os.path.join(dir,fas), MHC=MHC_binary, allele=HLA_allele, filter=Report_filter, size=Peptide_size, nM=nM_cutoff)
			tmp_file_cnt=tmp_file_cnt+1
			#out_fn=os.path.join(dir,fas) + '.netmhccons_out'
			finaldf = merge_df(snv, out, out_fn, HLA_allele)
			return finaldf
	if tmp_file_cnt == 0:
		print ("No fasta files were found for running netMHCcons\n")
		quit()
			
	
def merge_df(SNV, epi_df, fn, allele):
	data_pd = SNV
	#print(list(data_pd))
	### correct name for the data frame for PCGP
	data_pd=PCGP_correct_name(data_pd)
	### rename columns
	#data_pd_sel = data_pd[['annovar_gene', 'Sample', 'Chr', 'Position_hg19', 'class', 'annovar_peptide_change', 'ReferenceAllele', 'MutantAllele', 'annovar_acc']]
	data_pd_sel = data_pd[['GeneName', 'Sample', 'Chr', 'Position_hg19', 'Class', 'AAChange', 'ReferenceAllele', 'MutantAllele', 'mRNA_acc']]
	data_pd_sel.columns = ['GeneName', 'Sample', 'Chr', 'Position_hg19', 'Class', 'AAChange', 'ReferenceAllele', 'MutantAllele', 'mRNA_acc']
	#data_pd_sel = data_pd_sel[data_pd_sel.Class=='missense']
	
	data_pd_sel['GeneName']=data_pd_sel['GeneName'].str.replace('.','_')
	print(epi_df)
	epi_df['GeneName'], epi_df['Mutation'] = epi_df['ID'].str.rsplit('_',1).str
	epi_df_sel=epi_df.drop('Binding',1)
	
	
	merged=data_pd_sel.merge(epi_df_sel, how='right', sort=True, left_on=['GeneName', 'AAChange'], right_on=['GeneName' ,'Mutation'])
	merged=merged.drop(['AAChange', 'Position'],1)
	#merged.sort(['Mutation', 'nM'], inplace=True, ascending=True)
	merged=merged.rename(columns={'ID':'Gene_variant'})
	merged['HLAtype']=allele
	merged.drop_duplicates()
	merged.to_csv(fn, sep='\t', na_rep=".", header=True, index=False)
	return merged

def PCGP_correct_name(pd_df, field_name='Sample'):
	### PCGP sample mutation file has no tumor extension at the end of file
	for index, row in pd_df.iterrows():
		if not '_' in row[field_name]:
			pd_df.ix[index,field_name]=str(row[field_name])+'_D'
	return pd_df

#### functions for check expression of mutations, from Michel Edmenson and Rob Carter
def get_prot_position_to_genomic_codon_positions_from_mRNA(mrna_exons_all_df, mrna_acc, chromosome, debug=True):
	"""
	Returns a dict mapping protein position to genomic positions in 1-based coordinates from a single mRNA accession.

	-mrna_exons_all_df: dataframe derived from refflat file containing exon locations along the genome
	-mrna_acc: A singe mRNA accession 
	"""
	genomic_positions_of_mrna = []
	mrna_exons_df = get_one_row_df(mrna_acc, chromosome, mrna_exons_all_df)
	if not isinstance(mrna_exons_df, pd.core.frame.DataFrame):
		sys.stderr.write("Cannot find unique refflat entry for " + mrna_acc + " on chromosome " + chromosome + ". Skipping...")
		return None
	
	prot_pos_to_genomic_dict = {}
	es =	mrna_exons_df.loc[mrna_acc, 'exonStarts']
	ee =	mrna_exons_df.loc[mrna_acc, 'exonEnds']
	if "," in es or "," in ee:
		exon_starts_list = [int(i) for	i in es.split(",")]
		exon_ends_list = [int(i) for	i in ee.split(",")]
	else:
		exon_starts_list = [int(es)]
		exon_ends_list = [int(ee)]
	for i in range(0,len(exon_starts_list)):
		genomic_positions_of_mrna += range(exon_starts_list[i], exon_ends_list[i] + 1)
	genomic_positions_of_coding_sequence_set = set(range(mrna_exons_df.loc[mrna_acc,'cdsStart'], mrna_exons_df.loc[mrna_acc,'cdsEnd'] + 1)) & set(genomic_positions_of_mrna)
	genomic_positions_of_coding_sequence_list = list(genomic_positions_of_coding_sequence_set)
	prot_length = len(genomic_positions_of_coding_sequence_list)/3
	if mrna_exons_df.loc[mrna_acc,'strand'] == '+':
		genomic_positions_of_coding_sequence_list.sort()
	elif mrna_exons_df.loc[mrna_acc,'strand'] == '-':
		genomic_positions_of_coding_sequence_list = sorted(genomic_positions_of_coding_sequence_list, reverse=True)
	else:
		sys.exit("No valid strand!")
	for i in range(0,prot_length):
		prot_pos_to_genomic_dict[i+1] = genomic_positions_of_coding_sequence_list[i*3:(i+1)*3]
	#debug and sys.stderr.write("##########################################\n")
	#debug and sys.stderr.write("Mapping from protein to genomic coordinates:\n")
	#debug and sys.stderr.write("\n".join([(str(i) + "\t" + ','.join([str(k) for k in j])) for i,j in prot_pos_to_genomic_dict.iteritems()]))
	#debug and sys.stderr.write("\n##########################################\n")
	#print (prot_pos_to_genomic_dict)
	return prot_pos_to_genomic_dict

def get_genomic_coverage(sample, bams_path, strand, fastaobj, chromosome, pep_genomic_positions, genomic_pos, ref_allele, alt_allele, aa_ref, aa_alt, debug=True):
	'''
	This function determines various coverage statisitcs for a given peptide sequence overlapping a variant in a sample by querying an associated RNAseq 
	file. It returns a tuple with the following elements: ref_count: count of reference base, alt_count: count of alternative base, total_count: ######,
	fully_covered: ######, pep_sequence: peptide sequence)

	If bams_path is None, which will happen when there is no RNASeq data corresponding to mutation data, the tuple will be populated with None


	'''
	bases_counts_dict = {}
	pep_sequence = get_translated_seq_string_from_discontiguous_genomic_positions(fastaobj, strand, chromosome, pep_genomic_positions, genomic_pos, ref_allele, alt_allele)
	genomic_pos_str = ",".join([str(p) for p in pep_genomic_positions])
	
	if not bams_path:
		return (None, None, None, None, pep_sequence, False)

	covered=True
	b = pysam.AlignmentFile(bams_path, 'rb')
	for pc in b.pileup(str(chromosome), min(pep_genomic_positions) -1, max(pep_genomic_positions)):
		if pc.pos+1 in pep_genomic_positions:
			if sum([(not pr.is_del) and (pr.indel == 0 ) for pr in pc.pileups]) < 1:
				covered=False
			elif pc.pos+1 == genomic_pos:
				bases_counts_dict = count_base_type_at_position([ref_allele, alt_allele], pc)
	if not bases_counts_dict.keys():
		bases_counts_dict = {ref_allele:0, alt_allele:0, 'total':0}
	#debug and sys.stderr.write(str(bases_counts_dict) + "\n")	
	return (bases_counts_dict[ref_allele], bases_counts_dict[alt_allele], bases_counts_dict['total'], covered, pep_sequence, True)




def get_one_row_df(mrna_acc, chrom, df):
	'''
	Returns a oneDataframe corresponding to mrna_acc on chromosome chrom.
	
	input:
	-mrna_acc: mrna_accession in refflat file
	-chrom: chromosome location of mrna
	
	returns: a 1-row df if there exists a unique entry corresponding to mrna_acc and chrom, otherwise None
	'''
	if(mrna_acc in df.index):
		sub_df = df[df.index.isin([mrna_acc]) & df.chrom.isin([chrom])]
		if sub_df.shape[0] != 1:
			#print sub_df
			return None
		else:
			return sub_df
	else:
		return None

	return (bases_counts_dict[ref_allele], bases_counts_dict[alt_allele], bases_counts_dict['total'], covered, pep_sequence, True)



def count_base_type_at_position(bases_list, pu_col, debug=True):
	'''
	Returns a dict of the number of reads of each base aligned to the given pileup column, which includes all reads overlapping 
	genomic region
	
	-bases_list: bases to consider for overlap counting
	-PileupColumn object associated with a particular site in the genome
	'''
	base_count_dict = {}
	base_count_dict['total'] = 0
	for i in bases_list:
		base_count_dict[i] = 0
	bases_at_pos_list = []
	for r in pu_col.pileups:
		if (not r.is_del) and (r.indel == 0):
			bases_at_pos_list.append(r.alignment.query_sequence[r.query_position])
	for b in bases_at_pos_list:
		base_count_dict['total'] += 1
		if b in bases_list:
			base_count_dict[b] += 1
	return base_count_dict


def get_dict_of_discontiguous_bases_or_exons(sorted_int_list):
	exons_dict = {}
	exons_count = 0
	last_entry = None
	exons_dict[exons_count] = [sorted_int_list[0]]
	last_entry = sorted_int_list[0]
	for i in sorted_int_list[1:]:
		if i-last_entry != 1:
			exons_count += 1
			exons_dict[exons_count] = [i]
		else:
			exons_dict[exons_count] += [i]
		last_entry = i
	return exons_dict


def get_fasta_sequence_string_from_exons_dict(fasta_obj, strand, ref_name, exons_dict, genomic_pos, ref_allele, alt_allele, debug=False):
	seq_string = ''
	#translated_ref_string = Utils.get_bam_ref_query_string_by_genome_build(ref_name)
	#print ",".join([strand, str(ref_name), str(exons_dict), str(genomic_pos), ref_allele, alt_allele])
	ref_allele_len = len(ref_allele)
	
	genomic_pos_list = range(genomic_pos, genomic_pos + ref_allele_len)
	exons_traversal_list = range(0, len(exons_dict.keys()))
	if strand == '-':
		exons_traversal_list = list(reversed(exons_traversal_list))
	if len(list(exons_traversal_list)) > 1:
		logging.debug("Peptides overlapping {}->{} mutation at pos {} in {} span more than one exon".format(ref_allele, alt_allele, genomic_pos, ref_name))
	pos_overlap_set_flag=False
	for i in exons_traversal_list:
		exons_list_orig = exons_dict[i]
		exons_list = copy.copy(exons_list_orig)
		pos_overlap_set = set(genomic_pos_list) & set(exons_list)
		seq_frag = fasta_obj.fetch(ref_name, min(exons_list)-1, max(exons_list))
	#print("[Genomic pos] %s" % genomic_pos)
	#print("[Seq frag] %s" % seq_frag)
	#print("[pos_overlap_set] %s" % pos_overlap_set) 
		if pos_overlap_set:
			pos_overlap_set_flag=True
			#if strand == '+':
			for _genomic_pos_ind in range(0,len(genomic_pos_list)):
				genomic_pos = genomic_pos_list[_genomic_pos_ind]
				if genomic_pos in exons_list:
					ind = exons_list.index(genomic_pos)
					if seq_frag[ind].upper() != ref_allele[_genomic_pos_ind].upper():
						sys.exit("+'ve strand: Reference base {} does NOT match {}!".format(ref_allele, seq_frag[ind:(ind+ref_allele_len)]))
					else:
						seq_frag_list = list(seq_frag)
						#debug and sys.stderr.write("changing base " + seq_frag[ind] + " from " + ref_allele + " to" + alt_allele + "\n")
						seq_frag_list[ind] = alt_allele[_genomic_pos_ind]
						seq_frag = "".join(seq_frag_list).upper()
		if strand == '-':
			seq_frag = str(Seq(seq_frag).reverse_complement())
		elif strand != '+':
			sys.exit("unrecognized strand! Exiting.")
		seq_string += seq_frag
	if not pos_overlap_set_flag:
	#seq_string=''
		print("[Error] the mutation position was not in the exon regions %s, %s, %s" % (genomic_pos, ref_allele, alt_allele))
	#print("[Error exon dict] %s" % exons_dict)	
	return seq_string

range(0,1)

def translate_dna_from_string(dna_string):
	seq_entry = Seq(dna_string)
	trans_string = seq_entry.translate()
	return str(trans_string)

def reversecomplement_dna_from_string(dna_string):
	seq_entry = Seq(dna_string)
	trans_string = seq_entry.reverse_complement()
	return str(trans_string)


def get_translated_seq_string_from_discontiguous_genomic_positions(fastaobj, strand, chromosome, sequence_list, genomic_pos, ref_allele, alt_allele, debug=True):
	if strand == '-':
		sequence_list.sort()
	#print ("Sequence_list:%s " % str(sequence_list))
	exons_dict = get_dict_of_discontiguous_bases_or_exons(sequence_list)
	#print ("Exons dict: %s" % str(exons_dict))
	#print ("Sending the following ref allele: %s; alt allele: %s" % (ref_allele, alt_allele))
	dnaseq_string = get_fasta_sequence_string_from_exons_dict(fastaobj, strand, chromosome, exons_dict, genomic_pos, ref_allele, alt_allele)
	#debug and sys.stderr.write("Translating " + dnaseq_string+'\n')
	#print("Translated_seq: %s" % translate_dna_from_string(dnaseq_string))
	return translate_dna_from_string(dnaseq_string)

def get_df_of_peptide_and_genomic_rna_coverage_for_sample(muts_df, fastaobj, exons_df, inbam_filename, peptide_length=9):
	#This will be populated and THEN converted to a dataframe
	coverage_to_df_list = []
	muts_sub_df = muts_df
	muts_sample_grouped = muts_sub_df.groupby(['sample', 'variant'])
	for (sample,var),group in muts_sample_grouped:
	
		#print 'peptide_length is a ' + str(peptide_length)
		prot_to_genomic_dict	= get_prot_position_to_genomic_codon_positions_from_mRNA(exons_df, group.loc[group.index[0],'mrna_accession'], group.loc[group.index[0], 'chromosome'])
		#If there was no mapping from protein position to genomic coordinates for any reason, skip this entry
		if not prot_to_genomic_dict:
			logging.error("No prot_to_genomic_dict for sample {} and var {}".format(sample, var))
		#print "Map from peptide position to genomic positions"
		#print str(prot_to_genomic_dict)
		#sys.stderr.write(str(group.index[0]))
		#sys.stderr.write(str(group))
		strand = exons_df.loc[group.loc[group.index[0],'mrna_accession'], 'strand']
		if len(strand) > 1:
			strand = strand[0]
		if prot_to_genomic_dict:
			#print (group.loc[group.index[0],'mrna_accession'])
			prot_length = len(prot_to_genomic_dict.keys())
			genomic_pos = group.loc[group.index[0],'pos']
			prot_pos = group.loc[group.index[0],'protein_pos']
			if prot_pos + peptide_length -1 <= prot_length:
				if prot_pos >= peptide_length:
				#print ("1 Prot position: " + str(prot_pos) + ", genomic position: " + str(genomic_pos))
					for i in range(prot_pos - peptide_length + 1, prot_pos + 1):
						genomic_positions = []
						for j in range(i, i+peptide_length):
							genomic_positions += prot_to_genomic_dict[j]
						#print ("Current genomic positions :" + str(genomic_positions))
						(ref_count, alt_count, total_count, fully_covered, pep_sequence, expression_data) = get_genomic_coverage(sample, inbam_filename, strand, fastaobj, group.loc[group.index[0],'chromosome'], genomic_positions, genomic_pos, group.loc[group.index[0],'reference_allele'], group.loc[group.index[0], 'non_reference_allele'], group.loc[group.index[0],'aa_ref'], group.loc[group.index[0], 'aa_alt'])
						#print("\n%s" % pep_sequence)
						coverage_to_df_list.append({'sample':sample, 'variant':var, 'peptide':pep_sequence, 'pep_start':i, 'pep_end':i+peptide_length - 1, 'ns_prot_pos':prot_pos, 'fully_covered':fully_covered, 'ref_count':ref_count, 'mut_pos_count':total_count, 'alt_count':alt_count, 'expression_data':expression_data})
				else:
				#print ("2 Prot position: " + str(prot_pos) + ", genomic position: " + str(genomic_pos))
					for i in range(1, prot_pos + 1):
						genomic_positions = []
						for j in range(i, i+peptide_length):
							genomic_positions += prot_to_genomic_dict[j]
						(ref_count, alt_count, total_count, fully_covered, pep_sequence, expression_data) = get_genomic_coverage(sample, inbam_filename, strand, fastaobj, group.loc[group.index[0],'chromosome'], genomic_positions, genomic_pos, group.loc[group.index[0],'reference_allele'], group.loc[group.index[0], 'non_reference_allele'], group.loc[group.index[0],'aa_ref'], group.loc[group.index[0], 'aa_alt'])
						coverage_to_df_list.append({'sample':sample, 'variant':var, 'peptide':pep_sequence, 'pep_start':i, 'pep_end':i+peptide_length - 1, 'ns_prot_pos':prot_pos, 'fully_covered':fully_covered, 'ref_count':ref_count, 'mut_pos_count':total_count, 'alt_count':alt_count, 'expression_data':expression_data})
			else:
				if prot_pos >= peptide_length:
				#print ("3 Prot position: " + str(prot_pos) + ", genomic position: " + str(genomic_pos))
					for i in range(prot_pos - peptide_length + 1, prot_length-peptide_length + 1):
						genomic_positions = []
						for j in range(i, i+peptide_length):
							genomic_positions += prot_to_genomic_dict[j]
						(ref_count, alt_count, total_count, fully_covered, pep_sequence, expression_data) = get_genomic_coverage(sample, inbam_filename, strand, fastaobj, group.loc[group.index[0],'chromosome'], genomic_positions, genomic_pos, group.loc[group.index[0],'reference_allele'], group.loc[group.index[0], 'non_reference_allele'], group.loc[group.index[0],'aa_ref'], group.loc[group.index[0], 'aa_alt'])
						coverage_to_df_list.append({'sample':sample, 'variant':var, 'peptide':pep_sequence, 'pep_start':i, 'pep_end':i+peptide_length - 1, 'ns_prot_pos':prot_pos, 'fully_covered':fully_covered, 'ref_count':ref_count, 'mut_pos_count':total_count, 'alt_count':alt_count, 'expression_data':expression_data})
				else:
				#print ("4 Prot position: " + str(prot_pos) + ", genomic position: " + str(genomic_pos))
					for i in range(1, prot_length-peptide_length + 1):
						genomic_positions = []
						for j in range(i, i+peptide_length):
							genomic_positions += prot_to_genomic_dict[j]
						(ref_count, alt_count, total_count, fully_covered, pep_sequence, expression_data) = get_genomic_coverage(sample, inbam_filename, strand, fastaobj, group.loc[group.index[0],'chromosome'], genomic_positions, genomic_pos, group.loc[group.index[0],'reference_allele'], group.loc[group.index[0], 'non_reference_allele'], group.loc[group.index[0],'aa_ref'], group.loc[group.index[0], 'aa_alt'])
						coverage_to_df_list.append({'sample':sample, 'variant':var, 'peptide':pep_sequence, 'pep_start':i, 'pep_end':i+peptide_length - 1, 'ns_prot_pos':prot_pos, 'fully_covered':fully_covered, 'ref_count':ref_count, 'mut_pos_count':total_count,	'alt_count':alt_count, 'expression_data':expression_data})
			
		else:
			pass
	#print "mojo"
	#print coverage_to_df_list
	#print "magoo"
	#print muts_sub_df
	#logging.info("Merging data sources. Here are the heads of the to-be-merged dfs: %s\n%s\n".format(pd.DataFrame(coverage_to_df_list).head().to_string(), muts_sub_df.head().to_string())) 
	merged_df = pd.merge(pd.DataFrame(coverage_to_df_list), muts_sub_df, on=['sample', 'variant'], sort=False, how='left')
	return merged_df


def mut_df_io(data_pd):
	data_pd['variant']=data_pd.apply(lambda data_pd: '%s' % data_pd['AAChange'] if '_' in data_pd['AAChange'] else '%s_%s' % (data_pd['GeneName'], data_pd['AAChange']), axis=1)
	#pd_print_full(data_pd)
	data_pd_sel = data_pd[data_pd.Class.str.contains('missense', case=False)]
	#print(list(data_pd_sel))
	
	### WGS SNV
	#data_pd_sel = data_pd_sel[['ProteinGI', 'Sample', 'variant', 'Chr', 'Position_hg19', 'ReferenceAllele', 'MutantAllele', 'mRNA_acc']]
	#data_pd_sel.columns=['','sample', 'variant','chromosome', 'pos', 'reference_allele', 'non_reference_allele', 'mrna_accession']
	
	### Sith
	data_pd_sel = data_pd_sel[['Sample', 'variant', 'Chr', 'Position_hg19', 'ReferenceAllele', 'MutantAllele', 'mRNA_acc']]
	data_pd_sel.columns=['sample', 'variant','chromosome', 'pos', 'reference_allele', 'non_reference_allele', 'mrna_accession']	

	valid_variants_series = data_pd_sel['variant'].apply(lambda x: bool(re.search("_[A-Z][0-9]+[A-Z]$", x)))
	if valid_variants_series.size != data_pd_sel.shape[0]:
		print("Tossing %s mutants because they have invalid entries. %s remain.", data_pd_sel.shape[0]-valid_variants_series.size, valid_variants_series.size)
		data_pd_sel = data_pd_sel[valid_variants_series]
	data_pd_sel['protein_pos'] = [int(re.search("\d+", i.split("_")[1]).group(0)) for i in list(data_pd_sel['variant'])]
	data_pd_sel['aa_ref'] = [re.search("^[A-Z]+", i.split("_")[1]).group(0) for i in list(data_pd_sel['variant'])]
	data_pd_sel['aa_alt'] = [re.search("[A-Z]+$", i.split("_")[1]).group(0) for i in list(data_pd_sel['variant'])]
	### remove version number for the mrna_accession
	data_pd_sel.mrna_accession.apply(str)	
	data_pd_sel['mrna_accession'] = data_pd_sel['mrna_accession'].map(lambda data_pd_sel: data_pd_sel.rstrip('\..*'))
	#pd_print_full(data_pd_sel)
	return(data_pd_sel)

def PCGP_mut_expression(muts_df, bam_build, RNA_bam=None):
	bam_build_to_ref_locations_dict = {
		'hg18':'/home/dnanexus/genome/hg18.fa',
		'hg19':'/home/dnanexus/genome/hg19.fa',
                #'hg19':'hg19.fa',
		'GRCh37-lite':'/home/dnanexus/genome/GRCh37-lite.fa',
                'HG19_Broad_variant':'/home/dnanexus/genome/Homo_sapiens_assembly19.fasta',
                'g1k-human-build37':'/home/dnanexus/genome/GRCh37-lite.fa',
                'NCBI36_WUGSC_variant':'/home/dnanexus/genome/NCBI36_WUGSC_variant.fa'
	}
	
	refflat_filename = data_io.get_refflat_filename(bam_build)
	reference_fa_filename = data_io.get_ref_fa_filename(bam_build, bam_build_to_ref_locations_dict)
	fastaobj = pysam.Fastafile(reference_fa_filename)
	refflat_df = data_io.read_refflat_as_df(refflat_filename, bam_build)
	
	### make sure the comptability of the dataframe for "get_df_of_peptide_and_genomic_rna_coverage_for_sample 
	print("Original")
	print(muts_df)
	sample_muts_df = data_io.ensure_refflat_muts_df_accession_compatibility(muts_df, refflat_df)
	print("Compatible")
	print(sample_muts_df)
	compatible_refflat_df = data_io.get_build_compatibile_refflat(refflat_df, bam_build)
	compatible_sample_muts_df = data_io.get_build_compatibile_muts_df(sample_muts_df, bam_build)
	
	print("Compatible sample")
	print(compatible_sample_muts_df)

	### extract exoerssion level for all
	expression_cov_df = get_df_of_peptide_and_genomic_rna_coverage_for_sample(compatible_sample_muts_df, fastaobj,	compatible_refflat_df, RNA_bam, peptide_length=int(args.size_of_peptide))
	#print(expression_cov_df)
	#print(list(expression_cov_df))
	expression_cov_df_sel = expression_cov_df[['sample', 'mrna_accession', 'variant', 'chromosome', 'pos', 'reference_allele', 'non_reference_allele', 'peptide', 'alt_count', 'ref_count', 'fully_covered', 'expression_data']]

	
	#print(list(expression_cov_df_sel))
	return(expression_cov_df_sel)

### ............... Main
def run(args):
	#clean_tmp()	
	data=f2ary(args.snv_input, head=args.header, reannotate=False)
	### get_genes(data,psize=args.size_of_peptide, header=args.header,acc_col=10, mut_col=2) ### Rob's sith file
	### current SNV output format
	fn_ls, data=get_genes(data,psize=int(args.size_of_peptide))

	### batch run
	#net_df=all_netMHCcons(data, args.output, MHC_binary=args.netMHCcons_bin, HLA_allele=args.HLA_alleles, Report_filter=args.filter_report, Peptide_size=args.size_of_peptide, nM_cutoff=args.nM_cutoff)
	### single run
	if len(fn_ls) == 1:
		print("[Single-sample run]")
		out=netMHCcons(fn_ls[0], MHC=args.netMHCcons_bin, allele=args.HLA_alleles, filter=args.filter_report, size=args.size_of_peptide, nM=args.nM_cutoff)
		#data['GeneName']=data['GeneName'].str.replace('.','_')
		net_df = merge_df(data, out, args.output, args.HLA_alleles)
	elif len(fn_ls) > 1:
		print("[Multi-sample run]")
		net_df=all_netMHCcons(data, args.output, MHC_binary=args.netMHCcons_bin, HLA_allele=args.HLA_alleles, Report_filter=args.filter_report, Peptide_size=args.size_of_peptide, nM_cutoff=args.nM_cutoff)
	else:
		print("[Error] No fasta file for netMHCcons was defined")
		quit()
	print("[step] create mutation dataframe")
	mutations_df=mut_df_io(data)
	if mutations_df.empty:
		print ("[Error] Mutations data was not imported properly and empty...exit\n")
		quit()
	#print("[step] check RNA bam version")
        #bam_build_version=Utils.get_genome_build(args.RNA_bam)

	if args.RNA_bam is not None and len(args.RNA_bam) >=1:
		print("[step] check RNA bam version")
        	bam_build_version=Utils.get_genome_build(args.RNA_bam)
		if len(bam_build_version)==0:
			print ("[Error] Bam build version unknown ...exit\n")
			quit()
		#print("RNAseq bam version: %s" % bam_build_version)
			
		expression_df=PCGP_mut_expression(mutations_df, RNA_bam=args.RNA_bam, bam_build=bam_build_version)
		#expression_df.to_csv('test.result',	sep="\t", na_rep=".", header=True, index=False)
		print(expression_df)

		if expression_df.empty:
			print ("[Error] Expression data was not imported properly and empty...exit\n")
			quit()
		
		print("Ready for merge")
		print(list(expression_df))
		print(list(net_df))
		#print(expression_df)
		
		### final merge
		merge=net_df.merge(expression_df, left_on=['Sample', 'Peptide'], right_on=['sample', 'peptide'], how='left')
		merge.to_csv(args.expression_out, sep="\t", na_rep=".", header=True, index=False)
	else:
		print("[step] no RNA bam")
		bam_build_version='hg19'
		expression_df=PCGP_mut_expression(mutations_df, bam_build_version, RNA_bam=None)
		#expression_df.to_csv('test.result',	sep="\t", na_rep=".", header=True, index=False)
		if expression_df.empty:
			print ("[Error] Expression data was not imported properly and empty...exit\n")
			quit()

		print("Ready for merge")
		print(list(expression_df))
		#print(list(net_df))
		#print(expression_df)

		### final merge
		merge=pd.merge(net_df, expression_df, left_on=['Sample', 'Peptide'], right_on=['sample', 'peptide'], how='left')
		merge.to_csv(args.expression_out, sep="\t", na_rep=".", header=True, index=False)

if __name__=='__main__':
	parser=argparse.ArgumentParser(description='peptide affinity prediction based on the SNP')
	parser.add_argument('-i', '--snv_input', help='Input SNV file name (only missense mutations considered)', required=True)
	parser.add_argument('-o', '--output', help='Output with user defined peptide length, the epitope with the corresponding epitope affinity 500 nM', default="stout")
	parser.add_argument('-n', '--netMHCcons_bin', help='Path to netMHCcons binary', default="netMHCcons")
	parser.add_argument('-s', '--size_of_peptide', help='Peptide size for prediction of affinity', default=9)
	parser.add_argument('-hd', '--header', help='Row to use as header', default=1)
	parser.add_argument('-a', '--HLA_alleles', help='List of HLAtype of the sample', default='HLA*A:02:01')
	parser.add_argument('-fl', '--filter_report', help='Filter report based on the nM (-n) value set(0/1)', default=0)
	parser.add_argument('-nM', '--nM_cutoff', help='nM cutoff (default 500), only the epitopes with less than the cutoff will be reported', default=500)
	parser.add_argument('-rB', '--RNA_bam', help='Bam file', default=None, required=False)
	parser.add_argument('-F', '--genome_fasta_file', help='Reference genome fasta file')	
	parser.add_argument('-R', '--reference_flat_file', help='Reference flat file with gene annotation')
	parser.add_argument('-rO', '--expression_out', help='Output with expression status annotaed', default="stout")
	parser.add_argument('-log', '--run_log', help='Debug message for the run', default='run.log')

	args=parser.parse_args()
	
	print("Input: %s" % args.snv_input)
	print("Output: %s" % args.output)
	print("netMHCcons: %s" % args.netMHCcons_bin)
	print("Peptide size: %s" % str(args.size_of_peptide))
	print("Tested HLA allele: %s" % str(args.HLA_alleles))
	print("Header in SNV file: %s" % args.header)
	print("RNAbam: %s" % args.RNA_bam)
	if args.filter_report==1:
		print("Report filter: %d; cutoff: %d" % (args.filter_report,	args.nM_cutoff))

	logging.basicConfig(filename = args.run_log, level= 'DEBUG')
	run(args)
	
