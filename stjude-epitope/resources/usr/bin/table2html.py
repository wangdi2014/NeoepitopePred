#!/usr/bin/env python


import sys
import argparse
import numpy as np
import pandas as pd
import re
import os
import HTML
import time

def html(args):
	fls=[]
        main_html = open("Epitope_affinity_prediction.html", "w")
        all_html_together=sig_class_style()
        for file in os.listdir(args.path):
            if file.endswith(args.file_extension):
                print "Processing " + file
                sn=re.sub("\..*","", file)
                fn=os.path.join(args.path, file)

		#makedir(args.path+'/html')
                makedir(args.path+'/xlsx')

		#out_fn=args.path + '/html/' + fn +'.html'
		out_ex=args.path + '/xlsx/' + fn +'.xlsx'
                out_fn=args.path + fn +'.html'

                out_fn=re.sub("\.\/\/","./", out_fn)
		out_ex=re.sub("\.\/\/","./", out_ex)		
                out_fn=re.sub("\/\.\/","/", out_fn)
                out_ex=re.sub("\/\.\/","/", out_ex)    

                infile=pd.read_csv(fn, sep="\t", header=0)
		
                dataclass=''
		if 'snv' in file:
                        dataclass='Missense mutation'
		elif 'fusion' in file:
                        dataclass='Gene fusion'
                infile['Class']=dataclass
		
		if 'file' in infile:
                        infile=infile.rename(columns = {'file':'Sample'})
		if 'ID' in infile:
			infile=infile.rename(columns = {'ID':'Gene_variant'})               
	
                ### construct output table
                out_fn_ln='<a href="' + out_fn + '"' + sn + '> HTML report</a>'
                out_fn_bk='<a href="#%s">Link</a>'; out_fn_bk=out_fn_bk % (sn)
                out_ex_ln='<a href="' + out_ex + '"' + sn + '> Excel report</a>'
                out_cnt=infile[infile['nM'] <= int(args.affinity_highlight)].count()['nM']
                #fls.append([sn, str(out_cnt), out_fn_bk, out_ex_ln])
                fls.append([sn, str(out_cnt), out_fn_bk])

		#print infile	

		#report=infile.pivot_table(index=["Sample", "HLAtype","Gene_variant"], values=["Peptide", "1-log50k","nM", "Rank"],dropna=False, aggfunc = lambda x: " ".join([str(y) for y in x]))
                report=infile.pivot_table(index=["Sample", "HLAtype","Gene_variant", "Class", "Peptide", "Allele"], 
                                          values=["1-log50k","nM", "Rank"],
                                          dropna=True)

		#print report
		# simple prinit
		cls=html_head2(sn) + sig_class_style()
                
                htmlpage=report.to_html(index=True, 
                                       formatters={'nM': lambda x: format_column(x, int(args.affinity_highlight))},
                                       escape=False)
		
                htmlpage_for_all=table_div(sn) + htmlpage + table_div_foot(sn)
                all_html_together += htmlpage_for_all
                
                htmlpage=cls+htmlpage
                
                #with open(out_fn, 'w') as f:
                #    f.write(htmlpage)
                #f.close
		
                # highlight print
                #with open(out_fn, 'w') as f:
		#    f.write(report.style.\
                #        applymap(highlight_vals, subset='nM').\
                #        set_table_attributes("border=1").\
                #        render()
                #    )
		
		### format output to excel
		writer = pd.ExcelWriter(out_ex, engine='xlsxwriter')
		#writer = pd.ExcelWriter(out_ex)
		report.to_excel(writer, sheet_name='Sheet1')
		workbook  = writer.book
                worksheet = writer.sheets['Sheet1']
		num_format = workbook.add_format({'num_format': '#,##0.00'})		
		# Set the column width and format
		worksheet.set_column('A:A', 14)
		worksheet.set_column('B:B', 60)
		worksheet.set_column('C:C', 20)
		worksheet.set_column('D:D', 16)
		worksheet.set_column('E:E', 16)
		worksheet.set_column('F:F', 16)  	
		writer.save()	
        ### create main table for link
        #print fls
	head=html_head()
        
        htmlcode=HTML.table(fls, header_row=HTML.TableRow(["Sample", "nM<="+str(args.affinity_highlight), "Table"], 
                                 bgcolor="#800000",
                                 attribs={'style': 'color:white'}
                                 ),
                            col_align=['left','center','left','left']
				
                            )
        foot=html_foot()        
	
        htmlcode=head + htmlcode + all_html_together + foot
        #print htmlcode
        main_html.write(htmlcode)
	main_html.close()
		

def highlight(s):
    print s.index.get_level_values('nM')
    is_mos = s.index.get_level_values('nM')<=20000
    return ['color: darkorange' if v else 'color: darkblue' for v in is_mos]

def highlight_vals(val, min=0, max=500, bcolor='green', fcolor='white'):
	if min < val < max:
		return 'background-color: %s' % bcolor
	else:
		return ''
def is_significant(value):
        return '<div class="significant">{}</div>'.format(value)
        #return '<span class="significant">{}</span>'.format(value)

def sig_class_style():
	style_class='''
        <style>
        div.significant {
            background-color: green;
            color: white;
        }

        div.insignificant {
            background-color: white;
            color: black;
        }
        </style>\n'''
        return style_class

def html_head():
        head="""
        <!DOCTYPE>
        <html>
        <body>
        <title> Epitope Affinity </title>
        <h1 style="color:Black"> Affinity prediction of epitope </h1>
        <hr color=#800000>
        <p id='top'>Please click the links below to access the results for each sample:</p>
        """
        return head

def html_head2(s='sample'):
        head="""
        <!DOCTYPE>
        <html>
        <body>
        <title>Epitope affinity prediction for %s</title>
        <h2 style="color:#800000" align="right"> %s</h2>
        <hr color=#800000>
        <div>
         <a href="#sample" align="right"">%s</a>
        </div>
        """
        head= head % (s, s, s)
        return head

def table_div(s='sample'):
        div="""
        <p id=%s style="color:#800000" align="right"> %s</p>
        <hr color=#800000>
        """
        div= div % (s,s)
        return div


def html_foot():
        foot="""
        <hr color=#800000>
        <p align="right"> The data was generated on %s.</p>
        </body>
        </html>"""
        foot = foot % (time.strftime("%c"))
	return foot

def table_div_foot(s='sample'):
       div="""
       <p style="text-align:center;"><a href="#" class="backtotop">Back to Top &uarr;</a></p>
       """
       return div


def format_column(value, cutoff):
        '''
        format function applied to all elements in a selected column
        '''
        result=str(value)
        if value <= cutoff:
            result=is_significant(result)
        return result

def makedir(s):
    if not os.path.exists(s):
        os.makedirs(s)


if __name__=='__main__':
        parser=argparse.ArgumentParser(description='convert table to html')
        parser.add_argument('-p', '--path', help='Path to files that need to be convered', required=True)
        parser.add_argument('-e', '--file_extension', help='File extension', required=True)
        parser.add_argument('-c', '--affinity_highlight', help='The IC50 lower than the cutoff will be highlighted', required=False, default=500)

        args=parser.parse_args()

        print("Path: %s" % args.path)
        print("File extension: %s" %args.file_extension)
        print("Highlight threshold: %s" % str(args.affinity_highlight))

        html(args)


