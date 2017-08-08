
# coding: utf-8

# In[2]:

import argparse
import re
from MutPepExpression import Utils


# In[6]:

def parse_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--ref_fa_filenames', default='', nargs='+', type=str, help='Change the default  locations of the fasta references. Format is space separated key-value pairs as follows: "bam_build1=location bam_build2=location2"')
    #parser.add_argument('refflat_file', help=
    #                    '''
    #                    RefFlat file containing refseq accessions. Must contain a header and have standard column names:
    #                   ''')
    parser.add_argument('mutations_file', help=
    '''
    TSV file of mutation data. The following columns are required:
    -sample:                eg) SJLUAD_G1
    -variant:               Must be in form ^[^_]+_[A-z][0-9]+[A-z]
                            eg) ZNF717_R44K
    -mrna_accession:        This must match the name column from the refflat file.
                            eg) NM_001005484
    -pos:                   eg) 12990152
    -chromosome:            eg) 3
    -class:                 Must be 'MISSENSE'. Program will stop if any other values are present
    -reference_allele:      eg) T
    -non_reference_allele:  eg) A
    -assembly:              Must be ONE of ['hg18','hg19','GRCh37-lite'] for all rows
    
    The chromosome name will be adjusted to match the assembly by adding or removing 'chr' prefixes. 
    If the assembly column has more than one value or is not one of the specified values, 
    the program will exit.
    
    ''')
    parser.add_argument('sample', type = str, help='Sample name. Eg) SJMEL0001_D')
    #parser.add_argument('fasta_file', help='')
    #parser.add_argument('assembly_type', choices=['hg18','hg19','GRCh37-lite'], help='')
    parser.add_argument('outfile_name', help='')
    parser.add_argument('log_filename', default= 'logfile.txt', help='Name of the logfile')
    parser.add_argument('pep_length', type=int, help='')
    parser.add_argument('hla_alleles', type = str, nargs='+', help='list of HLA alleles for sample')
    parser.add_argument('--inbam_filename', default = None, help='Name of RNASeq bam file. Leave blank if none available and epitopes will be determined from GRCH37-lite without expression data')
    parser.add_argument('--no_allele_validation', action = 'store_true', help='Do not require 1 or 2 alleles at each locus. Any combination of alleles will be checked')
    args = parser.parse_args()
    #Update fa locations
    args.ref_fa_filenames = update_ref_fa_filenames(args)
    args.hla_alleles_dict =  get_and_validate_hla_dict(args)
    if is_valid_args(args):
        return args
    else:
        sys.exit("Arguments are not valid")


# In[3]:

def update_ref_fa_filenames(args):
    bam_build_to_ref_locations_dict = {
        'hg18':'/nfs_exports/genomes/1/Homo_sapiens/hg18/FASTA/Hg18.fa',
        'hg19':'/nfs_exports/genomes/1/Homo_sapiens/hg19/FASTA/Hg19.fa',
        'GRCh37-lite':'/nfs_exports/genomes/1/Homo_sapiens/GRCh37-lite/FASTA/GRCh37-lite.fa',
        'HG19_Broad_variant':'/nfs_exports/genomes/1/Homo_sapiens/Broad_hg19/FASTA/Homo_sapiens_assembly19.fasta',
        'g1k-human-build37':'/nfs_exports/genomes/1/Homo_sapiens/GRCh37-lite/FASTA/GRCh37-lite.fa',
        'NCBI36_WUGSC_variant':'/nfs_exports/genomes/1/Homo_sapiens/NCBI36_WUGSC_variant/FASTA/NCBI36_WUGSC_variant.fa'
    }
    if args.ref_fa_filenames:
        sys.exit("Not implemented yet! Only the default locations are useable!")
    else:
        return bam_build_to_ref_locations_dict


# In[5]:

def get_and_validate_hla_dict(args):
    if args.no_allele_validation:
        locus_to_alleles_dict = {}
        for _allele in args.hla_alleles:
            new_format = Utils.convert_hla_format(_allele, 'netmhcpan')
            _locus = re.sub("HLA-([A-Z])\d.+$", "\\1", new_format)
            if _locus in locus_to_alleles_dict:
                locus_to_alleles_dict[_locus].append(new_format)
            else:
                locus_to_alleles_dict[_locus] = [new_format]
    else:
        locus_to_alleles_dict = {'A':[], 'B':[], 'C':[]}
        for _allele in args.hla_alleles:
            new_format = Utils.convert_hla_format(_allele, 'netmhcpan')
            locus_to_alleles_dict[re.sub("HLA-([A-Z])\d.+$", "\\1", new_format)].append(new_format)
        for _locus_list in locus_to_alleles_dict.values():
            if not len(_locus_list) in [1,2]:
                sys.exit('Must have at least one allele for each locus and no more than two')
    return locus_to_alleles_dict
        


# In[ ]:

def is_valid_args(args):
    return True

