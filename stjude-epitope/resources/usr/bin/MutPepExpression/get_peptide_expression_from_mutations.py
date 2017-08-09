
# coding: utf-8

# In[1]:

#!/bin/env python


# In[2]:

from MutPepExpression import data_io
from MutPepExpression.args_handling import parse_args
#from MutPepExpression.args_handling import Utils
import pandas as pd
import numpy as np
import sys
import os.path
import pysam
import random
import re
from Bio.Seq import Seq
import copy
import string
import logging


# ### FUNCTION DEFINITIONS

# In[3]:

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
    es =  mrna_exons_df.loc[mrna_acc, 'exonStarts']
    ee =  mrna_exons_df.loc[mrna_acc, 'exonEnds']
    if "," in es or "," in ee:
      exon_starts_list = [int(i) for  i in es.split(",")]
      exon_ends_list = [int(i) for  i in ee.split(",")]
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
    debug and sys.stderr.write("##########################################\n")
    debug and sys.stderr.write("Mapping from protein to genomic coordinates:\n")
    debug and sys.stderr.write("\n".join([(str(i) + "\t" + ','.join([str(k) for k in j])) for i,j in prot_pos_to_genomic_dict.iteritems()]))
    debug and sys.stderr.write("\n##########################################\n")
    return prot_pos_to_genomic_dict

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
            print sub_df
            return None
        else:
            return sub_df
    else:
        return None


# In[4]:

def get_genomic_coverage(sample, bams_path, strand, fastaobj, chromosome, pep_genomic_positions, genomic_pos, ref_allele, alt_allele, debug=True):
  '''
  This function determines various coverage statisitcs for a given peptide sequence overlapping a variant in a sample by querying an associated RNAseq 
  file. It returns a tuple with the following elements: ref_count: count of reference base, alt_count: count of alternative base, total_count: ######,
  fully_covered: ######, pep_sequence: peptide sequence)
  
  
  '''
  b = pysam.Samfile(bams_path, 'rb')
  
  print "Bam path:"
  print bams_path
  covered=True
  bases_counts_dict = {}
  #-->
  print ",".join([sample, bams_path, strand, str(fastaobj), str(chromosome), str(pep_genomic_positions), str(genomic_pos), ref_allele, alt_allele])
  pep_sequence = get_translated_seq_string_from_discontiguous_genomic_positions(fastaobj, strand, chromosome, pep_genomic_positions, genomic_pos, ref_allele, alt_allele)
  #-->
  genomic_pos_str = ",".join([str(p) for p in pep_genomic_positions])
  debug and sys.stderr.write("{0}\t{1}\t{2}\t{3}\t{4}\t".format(pep_sequence, genomic_pos_str, genomic_pos, ref_allele, alt_allele))#Translated peptide from genomic positions (1-based coordinates) " + genomic_pos_str + ": " + str(pep_sequence) + "\n")
  print "looking for " + ",".join([str(chromosome), str(min(pep_genomic_positions) -1), str(max(pep_genomic_positions))])
  for pc in b.pileup(str(chromosome), min(pep_genomic_positions) -1, max(pep_genomic_positions)):
    if pc.pos+1 in pep_genomic_positions:
      if sum([(not pr.is_del) and (pr.indel == 0 ) for pr in pc.pileups]) < 1:
        covered=False
      elif pc.pos+1 == genomic_pos:
        bases_counts_dict = count_base_type_at_position([ref_allele, alt_allele], pc)
  if not bases_counts_dict.keys():
    bases_counts_dict = {ref_allele:0, alt_allele:0, 'total':0}
  debug and sys.stderr.write(str(bases_counts_dict) + "\n")  
  return (bases_counts_dict[ref_allele], bases_counts_dict[alt_allele], bases_counts_dict['total'], covered, pep_sequence)


# In[5]:

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
       bases_at_pos_list.append(r.alignment.seq[r.qpos])
  for b in bases_at_pos_list:
    base_count_dict['total'] += 1
    if b in bases_list:
      base_count_dict[b] += 1
  return base_count_dict


# In[6]:

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


# In[7]:

def get_fasta_sequence_string_from_exons_dict(fasta_obj, strand, ref_name, exons_dict, genomic_pos, ref_allele, alt_allele, debug=False):
  seq_string = ''
  #translated_ref_string = Utils.get_bam_ref_query_string_by_genome_build(ref_name)
  print ",".join([strand, str(ref_name), str(exons_dict), str(genomic_pos), ref_allele, alt_allele])
  for i in range(0, len(exons_dict.keys())):
    exons_list_orig = exons_dict[i]
    exons_list = copy.copy(exons_list_orig)
    if strand == '+':
      seq_frag = fasta_obj.fetch(ref_name, min(exons_list)-1, max(exons_list))
      if(genomic_pos in exons_list):
        ind = exons_list.index(genomic_pos)
        debug and sys.stderr.write("Comparing " +  seq_frag[ind] + " and " + ref_allele + "\n")
        if seq_frag[ind] != ref_allele:
          sys.exit("Reference base does NOT match!")
        else:
          seq_frag_list = list(seq_frag)
          debug and sys.stderr.write("changing base " + seq_frag[ind] + " from " + ref_allele + " to" + alt_allele + "\n")
          seq_frag_list[ind] = alt_allele
          seq_frag = "".join(seq_frag_list)
      seq_string += seq_frag
    elif strand == '-':
      #exons_list.sort()
      logging.debug("Grabbing %s %s-%s", ref_name, str(min(exons_list)-1), str(max(exons_list)))
      seq_frag_unrced_unmuted_str = fasta_obj.fetch(ref_name, min(exons_list)-1, max(exons_list))
      seq_frag=''
      if(genomic_pos in exons_list):
        ind = exons_list.index(genomic_pos)
        logging.debug("LOOKING FOR ref_allele:%s", ref_allele)
        logging.debug("LOOKING FOR SEQ_FRAG_UNRCED_UNMUTED:%s", seq_frag_unrced_unmuted_str)
        logging.debug("index: %s",str(ind))
        #debug and sys.stderr.write("Comparing " +  seq_frag_unrced_unmuted_str[ind] + " and " + ref_allele + "\n")
        if seq_frag_unrced_unmuted_str[ind] != ref_allele:
          sys.exit("Reference base does NOT match!")
        else:
          seq_frag_unrced_unmuted_str_list = list(seq_frag_unrced_unmuted_str)
          debug and sys.stderr.write("changing base " + seq_frag_unrced_unmuted_str_list[ind] + " from " + ref_allele + " to" + alt_allele + "\n")
          seq_frag_unrced_unmuted_str_list[ind] = alt_allele
          seq_frag = "".join(seq_frag_unrced_unmuted_str_list)
      else:
        seq_frag = seq_frag_unrced_unmuted_str
      seq_frag_rced_str = str(Seq(seq_frag).reverse_complement())
      seq_string += seq_frag_rced_str
    else:
      sys.exit('Wrong strand type')
  return seq_string


# In[8]:

def translate_dna_from_string(dna_string):
  seq_entry = Seq(dna_string)
  trans_string = seq_entry.translate()
  return str(trans_string)


# In[9]:

def get_translated_seq_string_from_discontiguous_genomic_positions(fastaobj, strand, chromosome, sequence_list, genomic_pos, ref_allele, alt_allele, debug=True):
    if strand == '-':
        sequence_list.sort()
    print "Sequence_list: " + str(sequence_list)
    exons_dict = get_dict_of_discontiguous_bases_or_exons(sequence_list)
    print "Exons dict: " + str(exons_dict)
    print "Sending the following ref allele: " + ref_allele
    dnaseq_string = get_fasta_sequence_string_from_exons_dict(fastaobj, strand, chromosome, exons_dict, genomic_pos, ref_allele, alt_allele)
    debug and sys.stderr.write("Translating " + dnaseq_string)
    return translate_dna_from_string(dnaseq_string)

def get_df_of_peptide_and_genomic_rna_coverage_for_sample(muts_df, fastaobj, exons_df, inbam_filename, peptide_length=9):
  #This will be populated and THEN converted to a dataframe
  coverage_to_df_list = []
  muts_sub_df = muts_df
  muts_sample_grouped = muts_sub_df.groupby(['sample', 'variant'])
  for (sample,var),group in muts_sample_grouped:
    
    print 'peptide_length is a ' + str(peptide_length)
    prot_to_genomic_dict  = get_prot_position_to_genomic_codon_positions_from_mRNA(exons_df, group.loc[group.index[0],'mrna_accession'], group.loc[group.index[0], 'chromosome'])
    #If there was no mapping from protein position to genomic coordinates for any reason, skip this entry
    if not prot_to_genomic_dict:
        continue
    #print "Map from peptide position to genomic positions"
    #print str(prot_to_genomic_dict)
    sys.stderr.write(str(group.index[0]))
    sys.stderr.write(str(group))
    strand = exons_df.loc[group.loc[group.index[0],'mrna_accession'], 'strand']
    if len(strand) > 1:
      strand = strand[0]
    if prot_to_genomic_dict:
      print group.loc[group.index[0],'mrna_accession']
      prot_length = len(prot_to_genomic_dict.keys())
      genomic_pos = group.loc[group.index[0],'pos']
      prot_pos = group.loc[group.index[0],'protein_pos']
      if prot_pos + peptide_length -1 <= prot_length:
        if prot_pos >= peptide_length:
          print "Prot position: " + str(prot_pos) + ", genomic position: " + str(genomic_pos)
          for i in range(prot_pos - peptide_length + 1, prot_pos + 1):
            genomic_positions = []
            for j in range(i, i+peptide_length):
              genomic_positions += prot_to_genomic_dict[j]
              print "Current genomic positions :" + str(genomic_positions)
            (ref_count, alt_count, total_count, fully_covered, pep_sequence) = get_genomic_coverage(sample, inbam_filename, strand, fastaobj, group.loc[group.index[0],'chromosome'], genomic_positions, genomic_pos, group.loc[group.index[0],'reference_allele'], group.loc[group.index[0], 'non_reference_allele'])
            coverage_to_df_list.append({'sample':sample, 'variant':var, 'peptide':pep_sequence, 'pep_start':i, 'pep_end':i+peptide_length - 1, 'ns_prot_pos':prot_pos, 'fully_covered':fully_covered, 'ref_count':ref_count, 'mut_pos_count':total_count, 'alt_count':alt_count})
        else:
          for i in range(1, prot_pos + 1):
            genomic_positions = []
            for j in range(i, i+peptide_length):
              genomic_positions += prot_to_genomic_dict[j]
            (ref_count, alt_count, total_count, fully_covered, pep_sequence) = get_genomic_coverage(sample, inbam_filename, strand, fastaobj, group.loc[group.index[0],'chromosome'], genomic_positions, genomic_pos, group.loc[group.index[0],'reference_allele'], group.loc[group.index[0], 'non_reference_allele'])
            coverage_to_df_list.append({'sample':sample, 'variant':var, 'peptide':pep_sequence, 'pep_start':i, 'pep_end':i+peptide_length - 1, 'ns_prot_pos':prot_pos, 'fully_covered':fully_covered, 'ref_count':ref_count, 'mut_pos_count':total_count, 'alt_count':alt_count})
      else:
        if prot_pos >= peptide_length:
          for i in range(prot_pos - peptide_length + 1, prot_length-peptide_length + 1):
            genomic_positions = []
            for j in range(i, i+peptide_length):
              genomic_positions += prot_to_genomic_dict[j]
            (ref_count, alt_count, total_count, fully_covered, pep_sequence) = get_genomic_coverage(sample, inbam_filename, strand, fastaobj, group.loc[group.index[0],'chromosome'], genomic_positions, genomic_pos, group.loc[group.index[0],'reference_allele'], group.loc[group.index[0], 'non_reference_allele'])
            coverage_to_df_list.append({'sample':sample, 'variant':var, 'peptide':pep_sequence, 'pep_start':i, 'pep_end':i+peptide_length - 1, 'ns_prot_pos':prot_pos, 'fully_covered':fully_covered, 'ref_count':ref_count, 'mut_pos_count':total_count, 'alt_count':alt_count})
        else:
          for i in range(1, prot_length-peptide_length + 1):
            genomic_positions = []
            for j in range(i, i+peptide_length):
              genomic_positions += prot_to_genomic_dict[j]
            (ref_count, alt_count, total_count, fully_covered, pep_sequence) = get_genomic_coverage(sample, inbam_filename, strand, fastaobj, group.loc[group.index[0],'chromosome'], genomic_positions, genomic_pos, group.loc[group.index[0],'reference_allele'], group.loc[group.index[0], 'non_reference_allele'])
            coverage_to_df_list.append({'sample':sample, 'variant':var, 'peptide':pep_sequence, 'pep_start':i, 'pep_end':i+peptide_length - 1, 'ns_prot_pos':prot_pos, 'fully_covered':fully_covered, 'ref_count':ref_count, 'mut_pos_count':total_count,  'alt_count':alt_count})
    else:
      pass
  print "mojo"
  print coverage_to_df_list
  print "magoo"
  print muts_sub_df
  merged_df = pd.merge(pd.DataFrame(coverage_to_df_list), muts_sub_df, on=['sample', 'variant'], sort=False, how='left')
  return merged_df

# ### Main

# In[11]:

def go(args):
  muts_df = data_io.read_missense_mutations_as_df(args.mutations_file)
  #exons_df = data_io.read_refflat_as_df(args.refflat_file, args.assembly_type)  
  fastaobj = pysam.Fastafile(args.fasta_file)
  genomic_cov_df = get_df_of_peptide_and_genomic_rna_coverage(muts_df, fastaobj, args.inbam_filename, args.pep_length)
  genomic_cov_df.to_csv(args.outfile_name, sep="\t", index=False)


# In[12]:

if __name__ == '__main__':
    args = parse_args()
    logging.basicConfig(filename = args.log_filename, level= 'DEBUG')
    go(args)


# ### Testing

# In[23]:

test=True
if test:
    sample_to_bam_df= pd.DataFrame({'path':[' /nfs_exports/genomes/1/projects/RNASEQ/TCGA_LungSkinOvarian/BucketRaw/SJLUAD/SJLUAD012528_D1-TCGA-05-4249-01A-01.bam']}, index=['SJLUAD012528_D1'])
    muts_df = data_io.read_missense_mutations_as_df("/nfs_exports/genomes/1/projects/WHOLEGENOME/TCGA_LungSkinOvarian/BucketIntermediate/HLATyping/delme_practice/SJLUAD012528_D1_muts.tsv")
    for (_sample,_bam_filename), group in sample_to_bam_df.groupby(['sample','bam']):
        bam_build = get_bam_build(_bam_filename)
        refflat_filename = get_refflat_filename(bam_build)
        reference_fa_filename = get_ref_fa_filename(bam_build)
        ref_fasta_obj = pysam.Fastafile(reference_fa_filename)
        
        test_refflat_df = data_io.read_refflat_as_df(refflat_filename, bam_build)
        sample_muts_df = muts_df[muts_df.sample == _sample]
        sample_muts_df = data_io.ensure_refflat_muts_df_accession_compatibility(sample_muts_df, test_refflat_df)
        compatible_refflat_df = data_io.get_build_compatibile_refflat(test_refflat_df, bam_build)
        compatible_sample_muts_df = data_io.get_build_compatibile_muts_df(sample_muts_df, bam_build)
        genomic_cov_df = get_df_of_peptide_and_genomic_rna_coverage_for_sample(compatible_sample_muts_df, ref_fasta_obj, compatible_refflat_df, _bam_filename, peptide_length=6)
    #print genomic_cov_df


# In[24]:

test_refflat_df


# In[25]:

test_muts_df.head()

