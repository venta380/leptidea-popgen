import sys
import string
import numpy as np
import fileinput
import pandas
import personal_popgen
from itertools import izip

#coverage_list_pop1=str(sys.argv[1])
#coverage_list_pop2=str(sys.argv[2])
#windows=str(sys.argv[3])
#output=str(sys.argv[4])

coverage_list_pop1=str(sys.argv[1])
coverage_list_pop2=str(sys.argv[2])
output=str(sys.argv[3])
dir_c=str(sys.argv[4])



############################# making windows ################################

chromosome=pandas.read_csv("/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private/final_fst_dxy_db")[['CHROM','BIN_START','BIN_END']]


#############################################################################

stats_dir="/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/"
list_1=[stats_dir+dir_c+stuff.strip()+".genomeCoverage.txt" for stuff in open(coverage_list_pop1)]
list_2=[stats_dir+dir_c+stuff.strip()+".genomeCoverage.txt" for stuff in open(coverage_list_pop2)]



def  assess_coverage_filter_7X_induv(list_1):
        ref_files = map(open, list_1)
        for line_1 in ref_files[0]:
                scaffold=line_1.strip().split()[0]
                position=line_1.strip().split()[1]
                cov1=[int(line_1.strip().split()[2])]
                cov2=[int(i.readline().strip().split()[2]) for i in ref_files[1:]]
                cover_list=cov1+cov2
                count=0
                for l in cover_list:
                        if l >= 7:
                                count+=1
                if count >= 7:
                        DESSISION="P"
                else:
                        DESSISION="F"
                yield [scaffold, int(position), DESSISION]


gen_gen=assess_coverage_filter_7X_induv(list_1)
gen_gen2=assess_coverage_filter_7X_induv(list_2)






chromosome[output+'_N_sites']=0


for i, j in izip(gen_gen, gen_gen2):        
        scaffold=i[0]
        position=i[1]
        DESSISION1=i[2]
        DESSISION2=j[2]
        if DESSISION1=="P" and DESSISION2=="P":
                chromosome.ix[(chromosome['CHROM'] == scaffold) & (chromosome['BIN_START'].astype(int) <= position) & (chromosome['BIN_END'].astype(int) >= position),  output+'_N_sites']+=1





chromosome.to_csv('/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/'+output+'_N_sites')

