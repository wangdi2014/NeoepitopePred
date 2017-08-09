
# coding: utf-8

# In[1]:

import sys
import re
import os
import logging
import pandas as pd
import pysam


# In[2]:

def get_bam_ref_query_string_by_genome_build(chromosome_series, genome_build):
    if genome_build in ['GRCh37-lite', 'HG19_Broad_variant']:
        #print genome_build
        return(chromosome_series.apply(lambda x: re.sub("^chr","",x)))
    elif genome_build in ['hg18', 'hg19']:
        return(chromosome_series.apply(lambda x: re.sub("(^chr)|(^)","chr",x)))
    else:
        sys.exit('Genome build "{}" is not recognized'.format(genome_build))


# In[3]:

def convert_hla_format(allele_str, output_format, digits = 4):
    temp_format = ''
    if digits != 4:
        sys.exit("only 4 digits is implemented")
    #if re.search("-", allele_str):
    temp_format = re.sub("HLA-", "", allele_str)
    #print temp_format
    temp_format = re.sub("^([A-Z])[_\*]?","\\1", temp_format)
    #print temp_format
    temp_format = re.sub("[_:]", "_", temp_format)
    #print temp_format
    if output_format == 'optitype':
        temp_format = re.sub("^([A-Z])","\\1*", temp_format)
        return re.sub("_",":", temp_format)
    if output_format == 'hlatyper':
        temp_format = re.sub("^([A-Z])","\\1_", temp_format)
        return re.sub("[\*:]","_", temp_format)
    if output_format == 'netmhcpan':
        temp_format = re.sub("^([A-Z])","HLA-\\1", temp_format)
        return re.sub("_",":", temp_format)
    else:
        sys.exit("Unknown hla format requested")
    


# In[4]:

def get_genome_build(bam_filename):
    '''
    Must return one of ['GRCh37-lite', 'HG19_Broad_variant', 'hg19'].
    If 'None' is passed as an argument, as will happen if no RNASeq data is present for a sample with mutation data, 
    'GRCh37-lite' will be returned.
    '''
    if not bam_filename:
        logging.info("No bam file. Returning GRCh37-lite as default")
        return 'GRCh37-lite'
    valid_genome_builds_dict = {'NCBI36_WUGSC_variant':'NCBI36_WUGSC_variant', 'g1k-human-build37':'GRCh37-lite', 'GRCh37-lite':'GRCh37-lite', 'HG19_Broad_variant':'HG19_Broad_variant', 'hg19':'hg19'}
    bam_obj = pysam.AlignmentFile(bam_filename, 'rb')
    bam_header_SQ_list = bam_obj.header['SQ']
    tag_list = []
    for _tag_dict in bam_header_SQ_list:
        for _tag in _tag_dict.keys():
            tag_list.append(_tag)
    tag_set = set(tag_list)
    if 'AS' in tag_set:
        as_list =  list(set([_entry['AS'] for _entry in bam_header_SQ_list]))
        if len(as_list) != 1:
            sys.exit('There should only be one "AS" field')
        else:
            genome_build = as_list[0]
            if genome_build in valid_genome_builds_dict:
                return valid_genome_builds_dict[genome_build]
            else:
                sys.exit("Genome build {} not recognized".format(genome_build))
    else:
        SN_to_ln_dict = {_entry['SN']:_entry['LN'] for _entry in bam_header_SQ_list if 'SN' in _entry}
        if any([re.match("chr", _ref_name) for _ref_name in SN_to_ln_dict.keys()]):
            if SN_to_ln_dict['chr1'] == 249250621:
                logging.warning("Returning genome build as hg19, although it could be a different variant of the GRCh37/hg19")
                return 'hg19'
            elif SN_to_ln_dict['chr1'] == 247249719:
                logging.warning("Returning genome build as hg18, although it could be a different variant of the GRCh36/hg18")
                return 'hg18'
            else:
                sys.exit("Unrecognized length for chromosome 1")
        else:
            if SN_to_ln_dict['1'] == 249250621:
                logging.warning("Returning genome build as GRCh37-lite, although it could be a different variant of the GRCh37/hg19")
                return 'GRCh37-lite'
            elif SN_to_ln_dict['1'] == 247249719:
            #    logging.warning("Returning genome build as hg18, although it could be a different variant of the GRCh36/hg18")
                logging.warning("Genome build is a variant of hg18 based on the length of chromosome 1, but does not contain the 'chr' prefix on the references. NCBI36_WUGSC_variant is being returned.")
                return('NCBI36_WUGSC_variant')
            else:
                sys.exit("Unrecognized length for chromosome 1")


# In[ ]:

#netmhcpan HLA-B40:98
#HLATyper A_02_02
#optitype A*01:01

