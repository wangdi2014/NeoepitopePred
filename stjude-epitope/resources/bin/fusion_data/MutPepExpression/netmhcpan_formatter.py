
# coding: utf-8

# In[1]:

import tempfile
import pandas as pd
from Bio import SeqIO
from Bio import SeqRecord
from Bio import Seq
import subprocess
import re


# In[2]:

class epitope_predictor:
    def __init__(self, peptide_series, hla_types_dict):
        self.peptides = list(peptide_series)
        self.A_alleles = self.ensure_valid_format(hla_types_dict['A'])
        self.B_alleles = self.ensure_valid_format(hla_types_dict['B'])
        self.C_alleles = self.ensure_valid_format(hla_types_dict['C'])
        self.all_alleles_list = self.A_alleles + self.B_alleles + self.C_alleles
        self.supported_human_class1_4digit_alleles = self.get_supported_human_class1_4digit_alleles()
        if not self.check_that_alleles_are_supported():
            sys.exit("Allele(s) NOT supported")
    
    def check_that_alleles_are_supported(self):
        return True
    
    def ensure_valid_format(self, alleles_list):
        '''
        This function ensures that any HLA allele is in the format required by netmhcpan, that is:
        4-digit resolution matching the following regular expression: HLA-[ABC]\d+:\d+.
        An example output is: HLA-C17:03
        '''
        reformatted_alleles_list = []
        for _allele in alleles_list:
            _temp_allele = re.sub("(^HLA-)|(^)", "HLA-", _allele)
            if re.search("_", _temp_allele):
                _temp_allele = re.sub("([ABC])_", "\1", _temp_allele)
                #_temp_allele = re.sub("_", "*", _temp_allele)
                _temp_allele = re.sub("_", ":", _temp_allele)
            reformatted_alleles_list.append(_temp_allele)
                
        return reformatted_alleles_list
        
    def get_supported_human_class1_4digit_alleles(self):
        all_supported_4digit_alleles_list = subprocess.check_output(['ssh', 'rcarter@erus', "netmhcpan -listMHC"])
        return([_allele for _allele in all_supported_4digit_alleles_list if re.match('HLA-[ABC]', _allele)])
        
    def predict_epitopes(self):
        pep_fa_file = self.get_peptide_fasta_filename()
        print pep_fa_file
        #for _allele in self.A_alleles:
        preds_filename = self.run_netmhcpan(pep_fa_file)
        preds_df = self.parse_preds_output(preds_filename)
        return preds_df
    
    def run_netmhcpan(self, pep_fa_file):
        temp_preds_handle = tempfile.NamedTemporaryFile(delete=False)
        temp_preds_filename = temp_preds_handle.name
        alleles_string_for_netmhcpan = ",".join(self.all_alleles_list)
        print 'executing ' + " ".join(['netmhcpan', '-a', alleles_string_for_netmhcpan, '-f', pep_fa_file, '-xls', '-xlsfile', temp_preds_filename])
        if subprocess.check_call(['netmhcpan', '-a', alleles_string_for_netmhcpan, '-f', pep_fa_file, '-xls', '-xlsfile', temp_preds_filename]):
            sys.exit("Problem running netmhcpan")
        return temp_preds_filename
    
    def parse_preds_output(self, preds_filename):
        #print preds_filename
        netmhcpan_output_df = pd.read_csv(preds_filename, sep="\t", header=1)
        netmhcpan_output_fh = open(preds_filename, 'r')
        alleles = re.sub("\t+", "\t", netmhcpan_output_fh.readline()).split("\t")[1:-1]
        netmhcpan_output_fh.close()
        df_list = []
        for _index in range(0,len(alleles)):
            _allele = alleles[_index]
            _pos = 3 + 3*_index
            #print "{} and {}".format(_pos, _pos+3)
            temp_df = pd.concat([netmhcpan_output_df.iloc[:,1:2], netmhcpan_output_df.iloc[:,(_pos):(_pos + 3 )]], axis = 1)
            temp_df.columns = [re.sub("\.\d+$", "", _col) for _col in temp_df.columns]
            temp_df['allele'] = _allele
            df_list.append(temp_df)
            #print temp_df
        final_df = pd.concat(df_list)
        return final_df
        
    
    def get_peptide_fasta_filename(self):
        temp_fasta_handle = tempfile.NamedTemporaryFile(delete=False)
        temp_fasta_name = temp_fasta_handle.name
        counter = 0
        seq_rec_list = []
        for _pep in self.peptides:
            print type(_pep)
            seq_rec_list.append(SeqRecord.SeqRecord(id= 's' + str(counter), seq = Seq.Seq(_pep)))
            counter += 1
        SeqIO.write(seq_rec_list, temp_fasta_handle, 'fasta')
        return temp_fasta_name


# In[1]:

test = False
if test:
    predder = epitope_predictor(pd.Series(['ASLLQLWQLW','LQLLQLWWLS','LWLLQLWSHE','DITELQLWYF','TEYFTLQLWG']),{'C':['HLA-C12:02','HLA-C08:01'],'B': ['HLA-B56:10'], 'A':['HLA-A03:01']})
    preds = predder.predict_epitopes()
    print preds

