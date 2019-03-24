import sys
import string
import personal_popgen
import itertools
from itertools import izip
import numpy as np
import subprocess
import os
import pandas


#####divide coverage files into three files#############

chromosome=pandas.read_csv("/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private/final_fst_dxy_db")[['CHROM','BIN_START','BIN_END']]

CHROMS=list(set(list(chromosome['CHROM'])))



files=subprocess.check_output([str(sys.argv[1])],shell=True).split()

#subprocess.call(["mkdir scaf1"], shell=True) 
#subprocess.call(["mkdir scaf2"], shell=True) 
#subprocess.call(["mkdir scaf3"], shell=True) 
#subprocess.call(["mkdir scaf4"], shell=True) 
#subprocess.call(["mkdir scaf5"], shell=True) 
#subprocess.call(["mkdir scaf6"], shell=True) 
#subprocess.call(["mkdir scaf7"], shell=True) 
#subprocess.call(["mkdir scaf8"], shell=True) 
#subprocess.call(["mkdir scaf9"], shell=True) 
#subprocess.call(["mkdir scaf10"], shell=True) 
#subprocess.call(["mkdir scaf11"], shell=True)
#subprocess.call(["mkdir scaf12"], shell=True) 
#subprocess.call(["mkdir scaf13"], shell=True) 
#subprocess.call(["mkdir scaf14"], shell=True) 
#subprocess.call(["mkdir scaf15"], shell=True) 
#subprocess.call(["mkdir scaf16"], shell=True) 
#subprocess.call(["mkdir scaf17"], shell=True) 
#subprocess.call(["mkdir scaf18"], shell=True) 



def split_shit_up(CHROMS,n):
	scaf1 =''
	scaf2 =''
	scaf3 =''
	scaf4 =''
	scaf5 =''
	scaf6 =''
	scaf7 =''
	scaf8 =''
	scaf9 =''
	scaf10=''
	scaf11=''
	scaf12=''
	scaf13=''
	scaf14=''
	scaf15=''
	scaf16=''
	scaf17=''
	scaf18=''

	for s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15, s16, s17, s18 in izip(np.array_split(CHROMS,18)[0], np.array_split(CHROMS,18)[1], np.array_split(CHROMS,18)[2], np.array_split(CHROMS,18)[3], np.array_split(CHROMS,18)[4], np.array_split(CHROMS,18)[5], np.array_split(CHROMS,18)[6], np.array_split(CHROMS,18)[7], np.array_split(CHROMS,18)[8], np.array_split(CHROMS,18)[9], np.array_split(CHROMS,18)[10], np.array_split(CHROMS,18)[11], np.array_split(CHROMS,18)[12], np.array_split(CHROMS,18)[13], np.array_split(CHROMS,18)[14], np.array_split(CHROMS,18)[15], np.array_split(CHROMS,18)[16], np.array_split(CHROMS,18)[17]):
		scaf1=scaf1+s1+"|"
		scaf2=scaf2+s2+"|"
		scaf3=scaf3+s3+"|"
		scaf4=scaf4+s4+"|"
		scaf5=scaf5+s5+"|"
		scaf6=scaf6+s6+"|"
		scaf7=scaf7+s7+"|"
		scaf8=scaf8+s8+"|"
		scaf9=scaf9+s9+"|"
		scaf10=scaf10+s10+"|"
		scaf11=scaf11+s11+"|"
		scaf12=scaf12+s12+"|"
		scaf13=scaf13+s13+"|"
		scaf14=scaf14+s14+"|"
		scaf15=scaf15+s15+"|"
		scaf16=scaf16+s16+"|"
		scaf17=scaf17+s17+"|"
		scaf18=scaf18+s18+"|"
	scaf1=scaf1[:-1]
	scaf2=scaf2[:-1] 
	scaf3=scaf3[:-1] 
	scaf4=scaf4[:-1] 
	scaf5=scaf5[:-1] 
	scaf6=scaf6[:-1] 
	scaf7=scaf7[:-1] 
	scaf8=scaf8[:-1] 
	scaf9=scaf9[:-1] 
	scaf10=scaf10[:-1]
	scaf11=scaf11[:-1]
	scaf12=scaf12[:-1]
	scaf13=scaf13[:-1]
	scaf14=scaf14[:-1]
	scaf15=scaf15[:-1]
	scaf16=scaf16[:-1]
	scaf17=scaf17[:-1]
	scaf18=scaf18[:-1]
	com1="grep -w -E '"+scaf1+"' "+n +"> ./scaf1"+"/"+n
	com2="grep -w -E '"+scaf2+"' "+n +"> ./scaf2"+"/"+n
	com3="grep -w -E '"+scaf3+"' "+n +"> ./scaf3"+"/"+n
	com4="grep -w -E '"+scaf4+"' "+n +"> ./scaf4"+"/"+n
	com5="grep -w -E '"+scaf5+"' "+n +"> ./scaf5"+"/"+n
	com6="grep -w -E '"+scaf6+"' "+n +"> ./scaf6"+"/"+n
	com7="grep -w -E '"+scaf7+"' "+n +"> ./scaf7"+"/"+n
	com8="grep -w -E '"+scaf8+"' "+n +"> ./scaf8"+"/"+n
	com9="grep -w -E '"+scaf9+"' "+n +"> ./scaf9"+"/"+n
	com10="grep -w -E '"+scaf10+"' "+n +"> ./scaf10"+"/"+n
	com11="grep -w -e '"+scaf11+"' "+n +"> ./scaf11"+"/"+n
	com12="grep -w -e '"+scaf12+"' "+n +"> ./scaf12"+"/"+n
	com13="grep -w -e '"+scaf13+"' "+n +"> ./scaf13"+"/"+n
	com14="grep -w -e '"+scaf14+"' "+n +"> ./scaf14"+"/"+n
	com15="grep -w -e '"+scaf15+"' "+n +"> ./scaf15"+"/"+n
	com16="grep -w -e '"+scaf16+"' "+n +"> ./scaf16"+"/"+n
	com17="grep -w -e '"+scaf17+"' "+n +"> ./scaf17"+"/"+n
	com18="grep -w -e '"+scaf18+"' "+n +"> ./scaf18"+"/"+n
	print scaf1
	print scaf2
	print scaf3
	print scaf4
	print scaf5
	print scaf6
	print scaf7
	print scaf8
	print scaf9
	print scaf10
	print scaf11
	print scaf12
	print scaf13
	print scaf14
	print scaf15
	print scaf16
	print scaf17
	print scaf18
	subprocess.call([com1], shell=True)
	subprocess.call([com2], shell=True)
	subprocess.call([com3], shell=True)
	subprocess.call([com4], shell=True)
	subprocess.call([com5], shell=True)
	subprocess.call([com6], shell=True)
	subprocess.call([com7], shell=True)
	subprocess.call([com8], shell=True)
	subprocess.call([com9], shell=True)
	subprocess.call([com10], shell=True)
	subprocess.call([com11], shell=True)
	subprocess.call([com12], shell=True)
	subprocess.call([com13], shell=True)
	subprocess.call([com14], shell=True)
	subprocess.call([com15], shell=True)
	subprocess.call([com16], shell=True)
	subprocess.call([com17], shell=True)
	subprocess.call([com18], shell=True)

for n in files:
	split_shit_up(CHROMS,n)







