#!/usr/bin/python
import os
import argparse
import commands
import pandas as pd
import string

from StringIO import StringIO

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

def run(args):
	out=netMHCcons(args.peptide_fasta_file, MHC=args.netMHCcons_bin, allele=args.HLA_alleles, filter=args.filter_report, size=args.size_of_peptide, nM=args.nM_cutoff)
	fn=args.peptide_fasta_file
	out['file']=fn.split('.', 1)[0]
	out['HLAtype']=args.HLA_alleles
	out_dedup=out.drop_duplicates()
	out_dedup.to_csv(args.output, sep='\t', na_rep=".", header=True, index=False)


if __name__=='__main__':
	parser=argparse.ArgumentParser(description='peptide affinity prediction using netMHCcons')
	parser.add_argument('-i', '--peptide_fasta_file', help='Input peptide fasta sequence file', required=True)
	parser.add_argument('-o', '--output', help='Output filename', default="stout")
	parser.add_argument('-n', '--netMHCcons_bin', help='Path to netMHCcons binary', default="netMHCcons")
	parser.add_argument('-s', '--size_of_peptide', help='Peptide size for prediction of affinity', default=9)
	parser.add_argument('-a', '--HLA_alleles', help='List of HLAtype of the sample', default='HLA*A:02:01')
	parser.add_argument('-fl', '--filter_report', help='Filter report based on the nM (-n) value set(0/1)', default=0)
	parser.add_argument('-nM', '--nM_cutoff', help='nM cutoff (default 500), only the epitopes with less than the cutoff will be reported', default=500)

	args=parser.parse_args()
	
	print("Input: %s" % args.peptide_fasta_file)
	print("Output: %s" % args.output)
	print("netMHCcons: %s" % args.netMHCcons_bin)
	print("Peptide size: %s" % str(args.size_of_peptide))
	print("Tested HLA allele: %s" % str(args.HLA_alleles))
	if args.filter_report==1:
                print("Report filter: %d; cutoff: %d" % (args.filter_report,  args.nM_cutoff))

	run(args)		

	
