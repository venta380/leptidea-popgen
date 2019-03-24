import sys
import string
import numpy as np
import fileinput
import pandas
import personal_popgen
from itertools import izip
import argparse
import os


pwd='/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/bam_covrage_files/'
os.chdir(pwd)


def get_args():
    parser = argparse.ArgumentParser(description='''outputs sites covered 1X in all Induveduals in population pairs

        -i file with list of paths to population lists
        -o output_prefix in database
        -p which part of the scafolds scaf1 scaf2 scaf3
        -n number of induvedluas in the set
        ''')

    parser.add_argument('-i', '--input', type=str, required=True)
    parser.add_argument('-o', '--output', type=str, required=True)
    parser.add_argument('-p', '--part', type=str, required=True)
    parser.add_argument('-n', '--num', type=int, required=True)
    return parser.parse_args()

args = get_args()
coverage_list_pop1=args.input
coverage_list_pop1=[stuff.strip() for stuff in open(coverage_list_pop1)]
output=args.output
inds=int(args.num)
############################# making windows ################################
chromosome=pandas.read_csv("/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/windows.csv")[['CHROM','BIN_START','BIN_END']]

#############################################################################



 
 
def  assess_coverage_filter_7X_induv(list_1, inds):
        headera=['CHROM', 'POS', 'COV0']
        keys=['COV0']
        list_all=[]
        for i in range(0,len(list_1)):
                if i == 0:
                        n=pandas.read_table(list_1[i], header=None, engine='c')
                        n.columns = headera
                        list_all=[n]
                else:
                        n=pandas.read_table(list_1[i], header=None, engine='c')
                        headera=['CHROM', 'POS', 'COV'+str(i)]
                        n.columns = headera
                        list_all.append(n)
                        keys.append('COV'+str(i))
                        nextone=personal_popgen.join_bam_coverage(list_all)
                        del list_all
                        list_all=[nextone]
        del list_all
        temp=nextone[keys]
        nextone=nextone[temp[temp >= 1].count(axis=1) >= inds]
        del temp
        return nextone[['CHROM','POS']]

stats_dir="/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/bam_covrage_files/"
list_1=[]
for i in coverage_list_pop1:
        n=[stats_dir+"/"+args.part+"/"+stuff.strip()+".genomeCoverage.txt" for stuff in open(i)]
        list_1.append(assess_coverage_filter_7X_induv(n, inds))

print "loaded all files"


df_pass=reduce(lambda x, y: pandas.merge(x, y, on=['CHROM','POS'], how='inner'), list_1)
print "merged all"

del list_1
df_pass['BIN_START']=(np.floor(df_pass['POS']/10000)*10000)+1
df_pass['BIN_END']=df_pass['BIN_START']+(10000-1)




df_merge = pandas.merge(df_pass, chromosome, on=['CHROM','BIN_START','BIN_END'],how='inner')
df_merge = df_merge.query('BIN_START <= POS and BIN_END >= POS')

print "computed windows"
del df_pass

out=pandas.DataFrame({output+'_N_sites' :df_merge.groupby(['CHROM','BIN_START','BIN_END']).size()}).reset_index()
print "computed output"



out.to_csv('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/bam_covrage_files'+"/"+args.part+"/"+output+'_N_sites')

print "completed"

