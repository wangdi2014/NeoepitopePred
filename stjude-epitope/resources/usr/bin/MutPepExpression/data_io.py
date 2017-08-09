
# coding: utf-8

# In[1]:

import pandas as pd
import numpy as np
from MutPepExpression import Utils
import re
import string
import sys
import logging


# In[2]:

def read_sample_to_bams_as_df(filename):
  sample_to_bams_df = pd.read_csv(filename, sep="\t", index_col=0, names=['loc'], dtype={'loc':str})
  return sample_to_bams_df


# In[3]:

def read_refflat_as_df(filename, assembly_type):
    '''
    Reads a refflat file and returns a Dataframe object. 
    
    filename = RefFlat filename
    assembly_type = string from the following list: ['hg18','hg19','GRCh37-lite']
    
    Note: Filename must contain header with the following column names, although order doesn't matter:
    bin,name,chrom,strand,txStart,txEnd,cdsStart,cdsEnd,exonCount,exonStarts,exonEnds
    The names column must correspond to the identifiers used in the mutation files (eg, refseq mRNA accessions)
    The assembly_type string is used to convert the chromosome to the correct string by adding or removing a 'chr' prefix 
    '''
    exons_df = pd.read_csv(filename, sep="\t", index_col = 'name', header=0, comment="#")#, names=['gene_name', 'chrom', 'strand', 'tx_start', 'tx_end', 'cds_start', 'cds_end', 'exon_count', 'exons_starts','exons_ends'], dtype={'gene_name':str, 'chrom':str, 'strand':str, 'tx_start':np.int64, 'tx_end':np.int64, 'cds_start':np.int64, 'cds_end':np.int64, 'exon_count':np.int64, 'exons_starts':str, 'exons_ends':str})
    exons_df['txStart'] = exons_df['txStart'] + 1
    exons_df['cdsStart'] = exons_df['cdsStart'] + 1
    exons_df['exonEnds'] = [string.rstrip(entry,",") for entry in list(exons_df['exonEnds'])]
    exons_df['exonStarts'] = [",".join([str(int(i) + 1) for i in string.rstrip(entry, ",").split(",")]) for entry in list(exons_df['exonStarts'])]
    #fix the name of the chromosomes according to the reference type
    #exons_df = fix_refflat_chromosome_column(exons_df, assembly_type)
    return exons_df


# In[5]:

def read_missense_mutations_as_df(filename):
    muts_df = pd.read_csv(filename, sep="\t", header=0, comment="#", dtype={'chromosome':str})
    #Remove variants with improper variant fileds
    #muts_df = pd.read_csv('/nfs_exports/genomes/1/projects/WHOLEGENOME/TCGA_LungSkinOvarian/BucketIntermediate/HLATyping/TCGA_MUTATION_DATA/all_diseases_level2_somatic_mutations_stjude_format.tsv', sep="\t", header=0, comment="#", dtype={'chromosome':str})
    valid_variants_series = muts_df['variant'].apply(lambda x: bool(re.search("_[A-Z][0-9]+[A-Z]$", x))
)
    if valid_variants_series.size != muts_df.shape[0]:
        logging.error("Tossing %s mutants because they have invalid entries. %s remain.", muts_df.shape[0]-valid_variants_series.size, valid_variants_series.size)
    muts_df = muts_df[valid_variants_series]
    #print muts_df[:3]
    #Kill the program if the mutation file chromosomes are not in the form 'chr[0-9]\+'
    muts_df['protein_pos'] = [int(re.search("\d+", i.split("_")[1]).group(0)) for i in list(muts_df['variant'])]
    muts_df['aa_ref'] = [re.search("^[A-Z]+", i.split("_")[1]).group(0) for i in list(muts_df['variant'])]
    muts_df['aa_alt'] = [re.search("[A-Z]+$", i.split("_")[1]).group(0) for i in list(muts_df['variant'])]
    #if is_valid_mutation_df(muts_df):
    return muts_df
    #else:
    #    sys.exit()


# In[7]:

def is_valid_mutation_df(mutation_df):
    '''
    Make sure the assembly column of the mutation file has values from the list:['hg18','hg19','GRCh37-lite', 'HG19_Broad_variant']
    '''
    #if len(mutation_df['assembly'].unique()) > 1: 
    #    sys.stderr.write("Assembly column of mutations_file should only have one value. It has the following values: " + ",".join([str(i) for i in list(mutation_df['assembly'].unique())]))
    #    return False
    if len(set(mutation_df['assembly']) | set(['hg18','hg19','GRCh37-lite', 'HG19_Broad_variant'])) != 4:
        sys.stderr.write("Assembly column of mutations_file should contain only one value from the following list: ['hg18','hg19','GRCh37-lite', 'HG19_Broad_variant']")
        return False
    else:
        return True


# In[48]:

def ensure_refflat_muts_df_accession_compatibility(mutation_df, refflat_df):
    tokeep_indices = mutation_df.mrna_accession.isin(refflat_df.index)
    logging.error("The following accessions from the sample mutations are not listed in the refflat file, so they are being dropped: %s", mutation_df[~tokeep_indices].to_string())
    return(mutation_df[tokeep_indices])


# In[ ]:

def get_build_compatibile_refflat(refflat_df, bam_build):
    refflat_df['chrom'] = Utils.get_bam_ref_query_string_by_genome_build(refflat_df['chrom'], bam_build)
    return(refflat_df)


# In[ ]:

def get_build_compatibile_muts_df(muts_df, bam_build):
    muts_df['chromosome'] = Utils.get_bam_ref_query_string_by_genome_build(muts_df['chromosome'], bam_build)
    return(muts_df)


# In[ ]:

def get_refflat_filename(bam_build):
    from pkg_resources import resource_filename
    if bam_build in ['GRCh37-lite', 'HG19_Broad_variant', 'hg19']:
        #print genome_build
        return(resource_filename('MutPepExpression', 'data/refflat_hg19_updated_08082015.tsv'))
    elif bam_build in ['hg18', 'NCBI36_WUGSC_variant']:
        return(resource_filename('MutPepExpression', 'data/refflat_hg18.txt'))
    else:
        sys.exit('Genome build "{}" is not recognized'.format(bam_build))


# In[ ]:

def get_ref_fa_filename(bam_build, build_to_fa_filename_dict):
    if bam_build in build_to_fa_filename_dict:
        return build_to_fa_filename_dict[bam_build]
    else:
        sys.exit('No record of a fasta reference for build {}'.format(bam_build))
    pass


# In[53]:

test=False
if test:
    logging.basicConfig()
    test_refflat_df = read_refflat_as_df("../examples/refflat_hg19_updated_08082015.tsv", 'hg19')
    #test_muts_df = read_missense_mutations_as_df("../examples/delme")
    #print test_refflat_df[:3]
    muts_df = read_missense_mutations_as_df("/nfs_exports/genomes/1/projects/WHOLEGENOME/TCGA_LungSkinOvarian/BucketIntermediate/HLATyping/TCGA_MUTATION_DATA/all_diseases_level2_somatic_mutations_stjude_format.tsv")
    #print muts_df.shape
    muts_mod_df = ensure_refflat_muts_df_compatibility(muts_df, test_refflat_df)
    #print muts_mod_df.shape
    #print muts_df[:3]
    


# In[ ]:



