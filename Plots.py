import sys
import os
import string
import pandas
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import seaborn as sns
import itertools
import math
import time
import sys
import personal_popgen
from pandas.tools.plotting import scatter_matrix
import matplotlib.cm as cm

sns.set(font_scale=1.5)
sns.set_style("whitegrid", {'axes.grid' : False})

pwd='/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/'
os.chdir(pwd)
def join_raw_data_base(lists):
        for i in range(1,len(lists)):
                j=i-1
                if i == 1:
                        new_2=pandas.merge(lists[j], lists[i], on=['CHROM','BIN_START'], how='inner')
                else:
                        new_2=pandas.merge(nextone, lists[i], on=['CHROM','BIN_START'], how='inner')
                nextone=new_2
        return nextone



irish_juvernica_kaz_sin=pandas.read_csv("irish_juvernica_kaz_sin_dxy_DB")
irish_juvernica_kazak_juvernica=pandas.read_csv("irish_juvernica_kazak_juvernica_dxy_DB")
irish_juvernica_spanish_reali=pandas.read_csv("irish_juvernica_spanish_reali_dxy_DB")
irish_juvernica_spanish_sinapis=pandas.read_csv("irish_juvernica_spanish_sinapis_dxy_DB")
irish_juvernica_swe_sin_allele=pandas.read_csv("irish_juvernica_swe_sin_allele_dxy_DB")
kaz_sin_spanish_reali=pandas.read_csv("kaz_sin_spanish_reali_dxy_DB")
kaz_sin_spanish_sinapis=pandas.read_csv("kaz_sin_spanish_sinapis_dxy_DB")
kaz_sin_swe_sin_allele=pandas.read_csv("kaz_sin_swe_sin_allele_dxy_DB")
kazak_juvernica_kaz_sin=pandas.read_csv("kazak_juvernica_kaz_sin_dxy_DB")
kazak_juvernica_spanish_reali=pandas.read_csv("kazak_juvernica_spanish_reali_dxy_DB")
kazak_juvernica_spanish_sinapis=pandas.read_csv("kazak_juvernica_spanish_sinapis_dxy_DB")
kazak_juvernica_swe_sin_allele=pandas.read_csv("kazak_juvernica_swe_sin_allele_dxy_DB")
spanish_reali_spanish_sinapis=pandas.read_csv("spanish_reali_spanish_sinapis_dxy_DB")
spanish_reali_swe_sin_allele=pandas.read_csv("spanish_reali_swe_sin_allele_dxy_DB")
spanish_sinapis_swe_sin_allele=pandas.read_csv("spanish_sinapis_swe_sin_allele_dxy_DB")
reali_juvernica=pandas.read_csv("reali_juvernica_dxy_DB")
sinapis_juvernica=pandas.read_csv("sinapis_juvernica_dxy_DB")
sinapis_reali=pandas.read_csv("sinapis_reali_dxy_DB")

sinapis_my_PI_thetaW_tajD=pandas.read_csv("sinapis_my_PI_thetaW_tajD")
spanish_sinapis_my_PI_thetaW_tajD=pandas.read_csv("spanish_sinapis_my_PI_thetaW_tajD")
swe_sin_allele_my_PI_thetaW_tajD=pandas.read_csv("swe_sin_allele_my_PI_thetaW_tajD")
kaz_sin_my_PI_thetaW_tajD=pandas.read_csv("kaz_sin_my_PI_thetaW_tajD")
reali_my_PI_thetaW_tajD=pandas.read_csv("reali_my_PI_thetaW_tajD")
spanish_reali_my_PI_thetaW_tajD=pandas.read_csv("spanish_reali_my_PI_thetaW_tajD")
kazak_juvernica_my_PI_thetaW_tajD=pandas.read_csv("kazak_juvernica_my_PI_thetaW_tajD")
juvernica_my_PI_thetaW_tajD=pandas.read_csv("juvernica_my_PI_thetaW_tajD")
irish_juvernica_my_PI_thetaW_tajD=pandas.read_csv("irish_juvernica_my_PI_thetaW_tajD")


sinapis_Avg_Cov=pandas.read_csv("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/sinapis_Avg_Cov.csv", names=['CHROM','BIN_START','BIN_END','sinapis_Avg_COV'], sep=' ',header=None)
spanish_sinapis_Avg_Cov=pandas.read_csv("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/spanish_sinapis_Avg_Cov.csv", names=['CHROM','BIN_START','BIN_END','spanish_sinapis_Avg_COV'], sep=' ',header=None)
swe_sin_allele_Avg_Cov=pandas.read_csv("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/swe_sin_allele_Avg_Cov.csv", names=['CHROM','BIN_START','BIN_END','swe_sin_allele_Avg_COV'], sep=' ',header=None)
kaz_sin_Avg_Cov=pandas.read_csv("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/kaz_sin_Avg_Cov.csv", names=['CHROM','BIN_START','BIN_END','kaz_sin_Avg_COV'], sep=' ',header=None)
reali_Avg_Cov=pandas.read_csv("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/reali_Avg_Cov.csv", names=['CHROM','BIN_START','BIN_END','reali_Avg_COV'], sep=' ',header=None)
spanish_reali_Avg_Cov=pandas.read_csv("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/spanish_reali_Avg_Cov.csv", names=['CHROM','BIN_START','BIN_END','spanish_reali_Avg_COV'], sep=' ',header=None)
kazak_juvernica_Avg_Cov=pandas.read_csv("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/kazak_juvernica_Avg_Cov.csv", names=['CHROM','BIN_START','BIN_END','kazak_juvernica_Avg_COV'], sep=' ',header=None)
juvernica_Avg_Cov=pandas.read_csv("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/juvernica_Avg_Cov.csv", names=['CHROM','BIN_START','BIN_END','juvernica_Avg_COV'], sep=' ',header=None)
irish_juvernica_Avg_Cov=pandas.read_csv("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/irish_juvernica_Avg_Cov.csv", names=['CHROM','BIN_START','BIN_END','irish_juvernica_Avg_COV'], sep=' ',header=None)




irish_juvernica_TajD=pandas.read_table("irish_juvernica_TajD.Tajima.D")
irish_juvernica_TajD=irish_juvernica_TajD.rename(index=str, columns={'TajimaD': 'irish_juvernica_TajimaD' })

juvernica_TajD=pandas.read_table("juvernica_TajD.Tajima.D")
juvernica_TajD=juvernica_TajD.rename(index=str, columns={'TajimaD': 'juvernica_TajimaD' })

kazak_juvernica_TajD=pandas.read_table("kazak_juvernica_TajD.Tajima.D")
kazak_juvernica_TajD=kazak_juvernica_TajD.rename(index=str, columns={'TajimaD': 'kazak_juvernica_TajimaD' })

kaz_sin_TajD=pandas.read_table("kaz_sin_TajD.Tajima.D")
kaz_sin_TajD=kaz_sin_TajD.rename(index=str, columns={'TajimaD': 'kaz_sin_TajimaD' })


reali_TajD=pandas.read_table("reali_TajD.Tajima.D")
reali_TajD=reali_TajD.rename(index=str, columns={'TajimaD': 'reali_TajimaD' })


sinapis_TajD=pandas.read_table("sinapis_TajD.Tajima.D")
sinapis_TajD=sinapis_TajD.rename(index=str, columns={'TajimaD': 'sinapis_TajimaD' })

spanish_reali_TajD=pandas.read_table("spanish_reali_TajD.Tajima.D")
spanish_reali_TajD=spanish_reali_TajD.rename(index=str, columns={'TajimaD': 'spanish_reali_TajimaD' })


spanish_sinapis_TajD=pandas.read_table("spanish_sinapis_TajD.Tajima.D")
spanish_sinapis_TajD=spanish_sinapis_TajD.rename(index=str, columns={'TajimaD': 'spanish_sinapis_TajimaD' })

swe_sin_TajD=pandas.read_table("swe_sin_TajD.Tajima.D")
swe_sin_TajD=swe_sin_TajD.rename(index=str, columns={'TajimaD': 'swe_sin_TajimaD' })


TajimaD=join_raw_data_base([irish_juvernica_TajD,juvernica_TajD,kazak_juvernica_TajD,kaz_sin_TajD,reali_TajD,sinapis_TajD,spanish_reali_TajD,spanish_sinapis_TajD,swe_sin_TajD])
TajimaD=TajimaD[['CHROM','BIN_START',"irish_juvernica_TajimaD","juvernica_TajimaD","kazak_juvernica_TajimaD","kaz_sin_TajimaD","reali_TajimaD","sinapis_TajimaD","spanish_reali_TajimaD","spanish_sinapis_TajimaD","swe_sin_TajimaD"]]
TajimaD['BIN_START']=TajimaD['BIN_START']+1
TajimaD['BIN_END']=TajimaD['BIN_START']+9999
dxy=[irish_juvernica_kaz_sin,irish_juvernica_kazak_juvernica,irish_juvernica_spanish_reali,irish_juvernica_spanish_sinapis,irish_juvernica_swe_sin_allele,kaz_sin_spanish_reali,kaz_sin_spanish_sinapis,kaz_sin_swe_sin_allele,kazak_juvernica_kaz_sin,kazak_juvernica_spanish_reali,kazak_juvernica_spanish_sinapis,kazak_juvernica_swe_sin_allele,spanish_reali_spanish_sinapis,spanish_reali_swe_sin_allele,spanish_sinapis_swe_sin_allele,reali_juvernica,sinapis_juvernica,sinapis_reali,sinapis_my_PI_thetaW_tajD,spanish_sinapis_my_PI_thetaW_tajD,swe_sin_allele_my_PI_thetaW_tajD,kaz_sin_my_PI_thetaW_tajD,reali_my_PI_thetaW_tajD,spanish_reali_my_PI_thetaW_tajD,kazak_juvernica_my_PI_thetaW_tajD,juvernica_my_PI_thetaW_tajD,irish_juvernica_my_PI_thetaW_tajD,TajimaD]
merge_dxy=personal_popgen.join_raw_data_base(dxy)
Avg_cov=personal_popgen.join_raw_data_base([sinapis_Avg_Cov,spanish_sinapis_Avg_Cov,swe_sin_allele_Avg_Cov,kaz_sin_Avg_Cov,reali_Avg_Cov,spanish_reali_Avg_Cov,kazak_juvernica_Avg_Cov,juvernica_Avg_Cov,irish_juvernica_Avg_Cov])



###################################################################################

data1=pandas.read_csv('final_DB_pi.csv')
data3=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/FINAL_db_N_sites_1x_10_ind')
#angsd_data=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/angsd/Merged_all_pestPG.csv')
#angsd_data.BIN_START=angsd_data.BIN_START+1
lists_files=[data1,data3,merge_dxy,Avg_cov]

data=personal_popgen.join_raw_data_base(lists_files)


chromosome_file=[stuff.split() for stuff in open('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private/chromosomes.txt')]
chromosomes_dict={}
for j in chromosome_file:
	if j[1] not in chromosomes_dict.keys():
		chromosomes_dict[j[1]]=[]
		chromosomes_dict[j[1]].append(j[0])
	else:
		chromosomes_dict[j[1]].append(j[0])


################################abba_baba_input####################################
abba_baba= pandas.read_table('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/angsd/ABBABABA2.abbababa', header=None,)[[0,1,2,7,8]]
abba_baba['D']=(abba_baba[7]-abba_baba[8])/(abba_baba[7]+abba_baba[8]) 
abba_baba=abba_baba.rename(index=str, columns={0: 'CHROM', 1: 'BIN_START', 2: 'BIN_END', 7: 'ABBA', 8: 'BABA'})

abba_baba=abba_baba[(abba_baba['ABBA'] + abba_baba['BABA'] >= 10) ]
for i in chromosomes_dict.keys():
	ans=abba_baba[(abba_baba['CHROM'].isin(chromosomes_dict[i]))]['D'].mean()
	std=abba_baba[(abba_baba['CHROM'].isin(chromosomes_dict[i]))]['D'].std()
	print i
	print ans
	print std





#################################################DXY_correct############################################

data['DXY_0']=(data['kaz_sin_swe_sin_allele_dxy']*(data['kaz_sin_swe_sin_allele_fixed']+data['kaz_sin_swe_sin_allele_private_a']+data['kaz_sin_swe_sin_allele_private_b']+data['kaz_sin_swe_sin_allele_shared']))/data["./kaz_sin_swe_sin_allele_N_sites"]
data['DXY_1']=(data['spanish_sinapis_swe_sin_allele_dxy']*(data['spanish_sinapis_swe_sin_allele_fixed']+data['spanish_sinapis_swe_sin_allele_private_a']+data['spanish_sinapis_swe_sin_allele_private_b']+data['spanish_sinapis_swe_sin_allele_shared']))/data["./spanish_sinapis_swe_sin_allele_N_sites"]
data['DXY_2']=(data['kaz_sin_spanish_sinapis_dxy']*(data['kaz_sin_spanish_sinapis_fixed']+data['kaz_sin_spanish_sinapis_private_a']+data['kaz_sin_spanish_sinapis_private_b']+data['kaz_sin_spanish_sinapis_shared']))/data["./kaz_sin_spanish_sinapis_N_sites"]
data['DXY_3']=(data['irish_juvernica_kazak_juvernica_dxy']*(data['irish_juvernica_kazak_juvernica_fixed']+data['irish_juvernica_kazak_juvernica_private_a']+data['irish_juvernica_kazak_juvernica_private_b']+data['irish_juvernica_kazak_juvernica_shared']))/data["./irish_juvernica_kazak_juvernica_N_sites"]
data['DXY_4']=(data['irish_juvernica_kaz_sin_dxy']*(data['irish_juvernica_kaz_sin_fixed']+data['irish_juvernica_kaz_sin_private_a']+data['irish_juvernica_kaz_sin_private_b']+data['irish_juvernica_kaz_sin_shared']))/data["./irish_juvernica_kaz_sin_N_sites"]
data['DXY_5']=(data['irish_juvernica_spanish_reali_dxy']*(data['irish_juvernica_spanish_reali_fixed']+data['irish_juvernica_spanish_reali_private_a']+data['irish_juvernica_spanish_reali_private_b']+data['irish_juvernica_spanish_reali_shared']))/data["./irish_juvernica_spanish_reali_N_sites"]
data['DXY_6']=(data['irish_juvernica_spanish_sinapis_dxy']*(data['irish_juvernica_spanish_sinapis_fixed']+data['irish_juvernica_spanish_sinapis_private_a']+data['irish_juvernica_spanish_sinapis_private_b']+data['irish_juvernica_spanish_sinapis_shared']))/data["./irish_juvernica_spanish_sinapis_N_sites"]
data['DXY_7']=(data['irish_juvernica_swe_sin_allele_dxy']*(data['irish_juvernica_swe_sin_allele_fixed']+data['irish_juvernica_swe_sin_allele_private_a']+data['irish_juvernica_swe_sin_allele_private_b']+data['irish_juvernica_swe_sin_allele_shared']))/data["./irish_juvernica_swe_sin_allele_N_sites"]
data['DXY_8']=(data['kazak_juvernica_kaz_sin_dxy']*(data['kazak_juvernica_kaz_sin_fixed']+data['kazak_juvernica_kaz_sin_private_a']+data['kazak_juvernica_kaz_sin_private_b']+data['kazak_juvernica_kaz_sin_shared']))/data["./kazak_juvernica_kaz_sin_N_sites"]
data['DXY_9']=(data['kazak_juvernica_spanish_reali_dxy']*(data['kazak_juvernica_spanish_reali_fixed']+data['kazak_juvernica_spanish_reali_private_a']+data['kazak_juvernica_spanish_reali_private_b']+data['kazak_juvernica_spanish_reali_shared']))/data["./kazak_juvernica_spanish_reali_N_sites"]
data['DXY_10']=(data['kazak_juvernica_spanish_sinapis_dxy']*(data['kazak_juvernica_spanish_sinapis_fixed']+data['kazak_juvernica_spanish_sinapis_private_a']+data['kazak_juvernica_spanish_sinapis_private_b']+data['kazak_juvernica_spanish_sinapis_shared']))/data["./kazak_juvernica_spanish_sinapis_N_sites"]
data['DXY_11']=(data['kazak_juvernica_swe_sin_allele_dxy']*(data['kazak_juvernica_swe_sin_allele_fixed']+data['kazak_juvernica_swe_sin_allele_private_a']+data['kazak_juvernica_swe_sin_allele_private_b']+data['kazak_juvernica_swe_sin_allele_shared']))/data["./kazak_juvernica_swe_sin_allele_N_sites"]
data['DXY_12']=(data['kaz_sin_spanish_reali_dxy']*(data['kaz_sin_spanish_reali_fixed']+data['kaz_sin_spanish_reali_private_a']+data['kaz_sin_spanish_reali_private_b']+data['kaz_sin_spanish_reali_shared']))/data["./kaz_sin_spanish_reali_N_sites"]
data['DXY_13']=(data['spanish_reali_spanish_sinapis_dxy']*(data['spanish_reali_spanish_sinapis_fixed']+data['spanish_reali_spanish_sinapis_private_a']+data['spanish_reali_spanish_sinapis_private_b']+data['spanish_reali_spanish_sinapis_shared']))/data["./spanish_reali_spanish_sinapis_N_sites"]
data['DXY_14']=(data['spanish_reali_swe_sin_allele_dxy']*(data['spanish_reali_swe_sin_allele_fixed']+data['spanish_reali_swe_sin_allele_private_a']+data['spanish_reali_swe_sin_allele_private_b']+data['spanish_reali_swe_sin_allele_shared']))/data["./spanish_reali_swe_sin_allele_N_sites"]

data['DXY_15']=(data['sinapis_reali_dxy']*(data['sinapis_reali_fixed']+data['sinapis_reali_private_a']+data['sinapis_reali_private_b']+data['sinapis_reali_shared']))/data["spanish_reali_N_sites"]
data['DXY_16']=(data['reali_juvernica_dxy']*(data['reali_juvernica_fixed']+data['reali_juvernica_private_a']+data['reali_juvernica_private_b']+data['reali_juvernica_shared']))/data["juvernica_reali_N_sites"]
data['DXY_17']=(data['sinapis_juvernica_dxy']*(data['sinapis_juvernica_fixed']+data['sinapis_juvernica_private_a']+data['sinapis_juvernica_private_b']+data['sinapis_juvernica_shared']))/data["sinapis_juvernica_N_sites"]


data['kaz_sin_swe_sin_allele_Pi_a']=(data['kaz_sin_swe_sin_allele_Pi_a']*(data['kaz_sin_swe_sin_allele_fixed']+data['kaz_sin_swe_sin_allele_private_a']+data['kaz_sin_swe_sin_allele_private_b']+data['kaz_sin_swe_sin_allele_shared']))/data["./kaz_sin_swe_sin_allele_N_sites"]
data['spanish_sinapis_swe_sin_allele_Pi_a']=(data['spanish_sinapis_swe_sin_allele_Pi_a']*(data['spanish_sinapis_swe_sin_allele_fixed']+data['spanish_sinapis_swe_sin_allele_private_a']+data['spanish_sinapis_swe_sin_allele_private_b']+data['spanish_sinapis_swe_sin_allele_shared']))/data["./spanish_sinapis_swe_sin_allele_N_sites"]
data['kaz_sin_spanish_sinapis_Pi_a']=(data['kaz_sin_spanish_sinapis_Pi_a']*(data['kaz_sin_spanish_sinapis_fixed']+data['kaz_sin_spanish_sinapis_private_a']+data['kaz_sin_spanish_sinapis_private_b']+data['kaz_sin_spanish_sinapis_shared']))/data["./kaz_sin_spanish_sinapis_N_sites"]
data['irish_juvernica_kazak_juvernica_Pi_a']=(data['irish_juvernica_kazak_juvernica_Pi_a']*(data['irish_juvernica_kazak_juvernica_fixed']+data['irish_juvernica_kazak_juvernica_private_a']+data['irish_juvernica_kazak_juvernica_private_b']+data['irish_juvernica_kazak_juvernica_shared']))/data["./irish_juvernica_kazak_juvernica_N_sites"]

data['sinapis_reali_Pi_a']=(data['sinapis_reali_Pi_a']*(data['sinapis_reali_fixed']+data['sinapis_reali_private_a']+data['sinapis_reali_private_b']+data['sinapis_reali_shared']))/data["spanish_reali_N_sites"]
data['reali_juvernica_Pi_a']=(data['reali_juvernica_Pi_a']*(data['reali_juvernica_fixed']+data['reali_juvernica_private_a']+data['reali_juvernica_private_b']+data['reali_juvernica_shared']))/data["juvernica_reali_N_sites"]
data['sinapis_juvernica_Pi_a']=(data['sinapis_juvernica_Pi_a']*(data['sinapis_juvernica_fixed']+data['sinapis_juvernica_private_a']+data['sinapis_juvernica_private_b']+data['sinapis_juvernica_shared']))/data["sinapis_juvernica_N_sites"]

data['kaz_sin_swe_sin_allele_Pi_b']=(data['kaz_sin_swe_sin_allele_Pi_b']*(data['kaz_sin_swe_sin_allele_fixed']+data['kaz_sin_swe_sin_allele_private_a']+data['kaz_sin_swe_sin_allele_private_b']+data['kaz_sin_swe_sin_allele_shared']))/data["./kaz_sin_swe_sin_allele_N_sites"]
data['spanish_sinapis_swe_sin_allele_Pi_b']=(data['spanish_sinapis_swe_sin_allele_Pi_b']*(data['spanish_sinapis_swe_sin_allele_fixed']+data['spanish_sinapis_swe_sin_allele_private_a']+data['spanish_sinapis_swe_sin_allele_private_b']+data['spanish_sinapis_swe_sin_allele_shared']))/data["./spanish_sinapis_swe_sin_allele_N_sites"]
data['kaz_sin_spanish_sinapis_Pi_b']=(data['kaz_sin_spanish_sinapis_Pi_b']*(data['kaz_sin_spanish_sinapis_fixed']+data['kaz_sin_spanish_sinapis_private_a']+data['kaz_sin_spanish_sinapis_private_b']+data['kaz_sin_spanish_sinapis_shared']))/data["./kaz_sin_spanish_sinapis_N_sites"]
data['irish_juvernica_kazak_juvernica_Pi_b']=(data['irish_juvernica_kazak_juvernica_Pi_b']*(data['irish_juvernica_kazak_juvernica_fixed']+data['irish_juvernica_kazak_juvernica_private_a']+data['irish_juvernica_kazak_juvernica_private_b']+data['irish_juvernica_kazak_juvernica_shared']))/data["./irish_juvernica_kazak_juvernica_N_sites"]

data['sinapis_reali_Pi_b']=(data['sinapis_reali_Pi_b']*(data['sinapis_reali_fixed']+data['sinapis_reali_private_a']+data['sinapis_reali_private_b']+data['sinapis_reali_shared']))/data["spanish_reali_N_sites"]
data['reali_juvernica_Pi_b']=(data['reali_juvernica_Pi_b']*(data['reali_juvernica_fixed']+data['reali_juvernica_private_a']+data['reali_juvernica_private_b']+data['reali_juvernica_shared']))/data["juvernica_reali_N_sites"]
data['sinapis_juvernica_Pi_b']=(data['sinapis_juvernica_Pi_b']*(data['sinapis_juvernica_fixed']+data['sinapis_juvernica_private_a']+data['sinapis_juvernica_private_b']+data['sinapis_juvernica_shared']))/data["sinapis_juvernica_N_sites"]


data['kaz_sin_swe_sin_allele_Avg_Pi']=(data.kaz_sin_swe_sin_allele_Pi_b+data.kaz_sin_swe_sin_allele_Pi_a)/2.0
data['spanish_sinapis_swe_sin_allele_Avg_Pi']=(data.spanish_sinapis_swe_sin_allele_Pi_b+data.spanish_sinapis_swe_sin_allele_Pi_a)/2.0
data['kaz_sin_spanish_sinapis_Avg_Pi']=(data.kaz_sin_spanish_sinapis_Pi_b+data.kaz_sin_spanish_sinapis_Pi_a)/2.0
data['irish_juvernica_kazak_juvernica_Avg_Pi']=(data.irish_juvernica_kazak_juvernica_Pi_b+data.irish_juvernica_kazak_juvernica_Pi_a)/2.0

data['sinapis_reali_Avg_Pi']=(data.sinapis_reali_Pi_b+data.sinapis_reali_Pi_a)/2.0
data['reali_juvernica_Avg_Pi']=(data.reali_juvernica_Pi_b+data.reali_juvernica_Pi_a)/2.0
data['sinapis_juvernica_Avg_Pi']=(data.sinapis_juvernica_Pi_b+data.sinapis_juvernica_Pi_a)/2.0




###########################################################################################################



###################here mention elements you want to plot##########################
N_VARIANTS=["N_VARIANTS_0","N_VARIANTS_1","N_VARIANTS_2","N_VARIANTS_3","N_VARIANTS_4","N_VARIANTS_5","N_VARIANTS_6","N_VARIANTS_7","N_VARIANTS_8","N_VARIANTS_9","N_VARIANTS_10","N_VARIANTS_11","N_VARIANTS_12","N_VARIANTS_13","N_VARIANTS_14"]



fst_all_15=["MEAN_FST_0","MEAN_FST_1","MEAN_FST_2","MEAN_FST_3","MEAN_FST_4","MEAN_FST_5","MEAN_FST_6","MEAN_FST_7","MEAN_FST_8","MEAN_FST_9","MEAN_FST_10","MEAN_FST_11","MEAN_FST_12","MEAN_FST_13","MEAN_FST_14"]
dxy_all_15=["DXY_0","DXY_1","DXY_2","DXY_3","DXY_4","DXY_5","DXY_6","DXY_7","DXY_8","DXY_9","DXY_10","DXY_11","DXY_12","DXY_13","DXY_14"]
FIXED_all_15=["FIXED_0","FIXED_1","FIXED_2","FIXED_3","FIXED_4","FIXED_5","FIXED_6","FIXED_7","FIXED_8","FIXED_9","FIXED_10","FIXED_11","FIXED_12","FIXED_13","FIXED_14"]


pi=['irish_juvernica_final_pi','kazak_juvernica_final_pi','spanish_reali_final_pi','kaz_sin_final_pi','spanish_sinapis_final_pi','swe_sin_final_pi']
pi_colours=sns.color_palette(['#C8C800','#006400','#0000FF','#FF8C00','#FF0000','#C04000'])
fst=['MEAN_FST_0','MEAN_FST_1','MEAN_FST_2','MEAN_FST_3','MEAN_FST_8','MEAN_FST_13']

angsd_pi=['Pi_irish_juvernica','Pi_kaza_juvernica','Pi_spanish_reali','Pi_kaza_sinapis','Pi_sinapis_spain','Pi_sweden_sinapis']

my_pi=['irish_juvernica_my_PI','kazak_juvernica_my_PI','spanish_reali_my_PI','kaz_sin_my_PI','spanish_sinapis_my_PI','swe_sin_allele_my_PI']

my_taj_D=['irish_juvernica_TajimaD','kazak_juvernica_TajimaD','spanish_reali_TajimaD','kaz_sin_TajimaD','spanish_sinapis_TajimaD','swe_sin_TajimaD']


fst_sinapis=['MEAN_FST_0', 'MEAN_FST_1','MEAN_FST_2']
dxy_sinapis=['DXY_0', 'DXY_1','DXY_2']
pi_sinapis=['kaz_sin_final_PI','spanish_sinapis_final_PI','swe_sin_final_PI']
colour_sin=['#FF8C00','#FF0000','#C04000']

fst_sinapis_reali=['MEAN_FST_0', 'MEAN_FST_1','MEAN_FST_2','MEAN_FST_13']
dxy_sinapis_reali=['DXY_0', 'DXY_1','DXY_2','DXY_13']


fst_sinapis_juver=['MEAN_FST_0', 'MEAN_FST_1','MEAN_FST_2','MEAN_FST_3']
dxy_sinapis_juver=['DXY_0', 'DXY_1','DXY_2','DXY_3']
pi_sinapis_juver=['kaz_sin_final_PI','spanish_sinapis_final_PI','swe_sin_final_PI','irish_juvernica_final_PI','kazak_juvernica_final_PI']
colour_sin_juve=['#FF8C00','#FF0000','#C04000','#C8C800','#006400']


fst_simp=['MEAN_FST_8','MEAN_FST_13']
dxy_simp=['DXY_8','DXY_13']
pi_simp=['kazak_juvernica_final_PI','kaz_sin_final_PI','spanish_reali_final_PI','spanish_sinapis_final_PI']
colour_simp=['#006400','#FF8C00','#0000FF','#FF0000']


fst_all=fst
dxy_all=['DXY_0','DXY_1','DXY_2','DXY_3','DXY_8','DXY_13']
colour_all=['#C8C800','#006400','#FF8C00','#0000FF','#FF0000','#C04000']


fst_new=['MEAN_FST_1','MEAN_FST_3']
dxy_new=['DXY_1','DXY_3']
pi_new=['irish_juvernica_final_PI','kazak_juvernica_final_PI','spanish_sinapis_final_PI','swe_sin_final_PI']
colour_new=['#C8C800','#006400','#FF0000','#C04000']




fst_comp_species=['MEAN_FST_15','MEAN_FST_16','MEAN_FST_17']
dxy_comp_species=['DXY_15','DXY_16','DXY_17']
col_comp_species=sns.color_palette(['#56AACC','#B6AAE8','#E88CD9'])
pi_species=['juvernica_my_PI','reali_my_PI','sinapis_my_PI']
pi_colours_species=sns.color_palette(['#00ff00','#00e1ff','#a80000'])
taj_comp_species=['juvernica_TajimaD','reali_TajimaD','sinapis_TajimaD']


FIXED_0=["kaz_sin_swe_sin_allele_fixed"]
FIXED_1=["spanish_sinapis_swe_sin_allele_fixed"]
FIXED_2=["kaz_sin_spanish_sinapis_fixed"]
FIXED_3=["irish_juvernica_kazak_juvernica_fixed"]
FIXED_4=["irish_juvernica_kaz_sin_fixed"]
FIXED_5=["irish_juvernica_spanish_reali_fixed"]
FIXED_6=["irish_juvernica_spanish_sinapis_fixed"]
FIXED_7=["irish_juvernica_swe_sin_allele_fixed"]
FIXED_8=["kazak_juvernica_kaz_sin_fixed"]
FIXED_9=["kazak_juvernica_spanish_reali_fixed"]
FIXED_10=["kazak_juvernica_spanish_sinapis_fixed"]
FIXED_11=["kazak_juvernica_swe_sin_allele_fixed"]
FIXED_12=["kaz_sin_spanish_reali_fixed"]
FIXED_13=["spanish_reali_spanish_sinapis_fixed"]
FIXED_14=["spanish_reali_swe_sin_allele_fixed"]

FIXED_15=["sinapis_reali_fixed"]
FIXED_16=["reali_juvernica_fixed"]
FIXED_17=["sinapis_juvernica_fixed"]

fsts=["kaz_sin_swe_sin_allele_fst","spanish_sinapis_swe_sin_allele_fst","kaz_sin_spanish_sinapis_fst","irish_juvernica_kazak_juvernica_fst","irish_juvernica_kaz_sin_fst","irish_juvernica_spanish_reali_fst","irish_juvernica_spanish_sinapis_fst","irish_juvernica_swe_sin_allele_fst","kazak_juvernica_kaz_sin_fst","kazak_juvernica_spanish_reali_fst","kazak_juvernica_spanish_sinapis_fst","kazak_juvernica_swe_sin_allele_fst","kaz_sin_spanish_reali_fst","spanish_reali_spanish_sinapis_fst","spanish_reali_swe_sin_allele_fst"]


sites_filter=data[(data['irish_juvernica_N_sites'] >= 3000) & (data['kazak_juvernica_N_sites'] >= 3000) & (data['kaz_sin_N_sites'] >= 3000) & (data['spanish_reali_N_sites'] >= 3000) & (data['spanish_sinapis_N_sites'] >= 3000)& (data['swe_sin_allele_N_sites'] >= 3000)& (data['sinapis_N_sites'] >= 3000)& (data['reali_N_sites'] >= 3000)& (data['juvernica_N_sites'] >= 3000)]
sites_filter_0=sites_filter[(sites_filter['irish_juvernica_final_pi'] <= 1) & (sites_filter['kazak_juvernica_final_pi'] <= 1) & (sites_filter['kaz_sin_final_pi'] <= 1) & (sites_filter['spanish_reali_final_pi'] <= 1) & (sites_filter['spanish_sinapis_final_pi'] <= 1) & (sites_filter['swe_sin_final_pi'] <= 1)]
sites_filter_0=sites_filter_0[(sites_filter_0['irish_juvernica_final_pi'] <= 1) & (sites_filter['kazak_juvernica_final_pi'] <= 1) & (sites_filter['kaz_sin_final_pi'] <= 1) & (sites_filter['spanish_reali_final_pi'] <= 1) & (sites_filter['spanish_sinapis_final_pi'] <= 1) & (sites_filter['swe_sin_final_pi'] <= 1)]
sites_filter_0.MEAN_FST_0=sites_filter_0.MEAN_FST_0.clip(lower=0)
sites_filter_0.MEAN_FST_1=sites_filter_0.MEAN_FST_1.clip(lower=0)
sites_filter_0.MEAN_FST_2=sites_filter_0.MEAN_FST_2.clip(lower=0)
sites_filter_0.MEAN_FST_3=sites_filter_0.MEAN_FST_3.clip(lower=0)
sites_filter_0.MEAN_FST_4=sites_filter_0.MEAN_FST_4.clip(lower=0)
sites_filter_0.MEAN_FST_5=sites_filter_0.MEAN_FST_5.clip(lower=0)
sites_filter_0.MEAN_FST_6=sites_filter_0.MEAN_FST_6.clip(lower=0)
sites_filter_0.MEAN_FST_7=sites_filter_0.MEAN_FST_7.clip(lower=0)
sites_filter_0.MEAN_FST_8=sites_filter_0.MEAN_FST_8.clip(lower=0)
sites_filter_0.MEAN_FST_9=sites_filter_0.MEAN_FST_9.clip(lower=0)
sites_filter_0.MEAN_FST_10=sites_filter_0.MEAN_FST_10.clip(lower=0)
sites_filter_0.MEAN_FST_11=sites_filter_0.MEAN_FST_11.clip(lower=0)
sites_filter_0.MEAN_FST_12=sites_filter_0.MEAN_FST_12.clip(lower=0)
sites_filter_0.MEAN_FST_13=sites_filter_0.MEAN_FST_13.clip(lower=0)
sites_filter_0.MEAN_FST_14=sites_filter_0.MEAN_FST_14.clip(lower=0)

sites_filter_0.MEAN_FST_15=sites_filter_0.MEAN_FST_15.clip(lower=0)    
sites_filter_0.MEAN_FST_16=sites_filter_0.MEAN_FST_16.clip(lower=0)    
sites_filter_0.MEAN_FST_17=sites_filter_0.MEAN_FST_17.clip(lower=0)  

sites_filter_0_pi=sites_filter_0[pi]
sites_filter_0_dxy=sites_filter_0[dxy_all_15]
sites_filter_0_fst=sites_filter_0[fst_all_15]
sites_filter_0_vcftools_fst=sites_filter_0[(sites_filter_0['MEAN_FST_0'] >=0) & (sites_filter_0['MEAN_FST_1'] >=0) & (sites_filter_0['MEAN_FST_2'] >=0) & (sites_filter_0['MEAN_FST_3'] >=0) & (sites_filter_0['MEAN_FST_8'] >=0) & (sites_filter_0['MEAN_FST_13'] >=0)][fst]
#########################overlay_chrmosomes#############################

fasta=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034_nobackup/GENOME_ASSEMBLY/assembly_updates/v1.4/N.Backstrom_leptidea.scf.1.4.fasta')

chromosome = pandas.DataFrame(columns=['CHROM','BIN_START','BIN_END'])

for i in fasta.keys():
        temp={}
        if len(fasta[i]) >= 10000:
                list_start=range(1,len(fasta[i]),10000)
                list_start.pop()
                list_end=range(10000,len(fasta[i]),10000)
        else:
                list_start=[1]
                list_end=[10000]
        chrom=[i]*len(list_end)
        temp['CHROM']=pandas.Series(chrom)
        temp['BIN_START']=pandas.Series(list_start)
        temp['BIN_END']=pandas.Series(list_end)
        temp_df = pandas.DataFrame(temp)
        chromosome=personal_popgen.join_data_base([chromosome,temp_df])


list_1=[sites_filter_0, chromosome]
data_Z_FST_2=personal_popgen.join_data_base(list_1)
data['CHROM_ID']='.'

data2=sites_filter_0
data2['CHROM_ID']='.'
data2['auto_allo']='auto'
for n in chromosomes_dict.keys():
	for nj in chromosomes_dict[n]:
		data2.ix[(data2['CHROM'] == nj), 'CHROM_ID']=str(n)
                data.ix[(data['CHROM'] == nj), 'CHROM_ID']=str(n)

data2.ix[(data2['CHROM_ID'] == '21'), 'auto_allo']='Z'

z=data2[(data2['auto_allo'] == 'Z')]
a=data2[(data2['auto_allo'] != 'Z')]


a['Z_FST_0']=(a['MEAN_FST_0']-np.nanmean(a['MEAN_FST_0']))/np.nanstd(a['MEAN_FST_0'])
a['Z_FST_1']=(a['MEAN_FST_1']-np.nanmean(a['MEAN_FST_1']))/np.nanstd(a['MEAN_FST_1'])
a['Z_FST_2']=(a['MEAN_FST_2']-np.nanmean(a['MEAN_FST_2']))/np.nanstd(a['MEAN_FST_2'])
a['Z_FST_3']=(a['MEAN_FST_3']-np.nanmean(a['MEAN_FST_3']))/np.nanstd(a['MEAN_FST_3'])
a['Z_FST_4']=(a['MEAN_FST_4']-np.nanmean(a['MEAN_FST_4']))/np.nanstd(a['MEAN_FST_4'])
a['Z_FST_5']=(a['MEAN_FST_5']-np.nanmean(a['MEAN_FST_5']))/np.nanstd(a['MEAN_FST_5'])
a['Z_FST_6']=(a['MEAN_FST_6']-np.nanmean(a['MEAN_FST_6']))/np.nanstd(a['MEAN_FST_6'])
a['Z_FST_7']=(a['MEAN_FST_7']-np.nanmean(a['MEAN_FST_7']))/np.nanstd(a['MEAN_FST_7'])
a['Z_FST_8']=(a['MEAN_FST_8']-np.nanmean(a['MEAN_FST_8']))/np.nanstd(a['MEAN_FST_8'])
a['Z_FST_9']=(a['MEAN_FST_9']-np.nanmean(a['MEAN_FST_9']))/np.nanstd(a['MEAN_FST_9'])
a['Z_FST_10']=(a['MEAN_FST_10']-np.nanmean(a['MEAN_FST_10']))/np.nanstd(a['MEAN_FST_10'])
a['Z_FST_11']=(a['MEAN_FST_11']-np.nanmean(a['MEAN_FST_11']))/np.nanstd(a['MEAN_FST_11'])
a['Z_FST_12']=(a['MEAN_FST_12']-np.nanmean(a['MEAN_FST_12']))/np.nanstd(a['MEAN_FST_12'])
a['Z_FST_13']=(a['MEAN_FST_13']-np.nanmean(a['MEAN_FST_13']))/np.nanstd(a['MEAN_FST_13'])
a['Z_FST_14']=(a['MEAN_FST_14']-np.nanmean(a['MEAN_FST_14']))/np.nanstd(a['MEAN_FST_14'])
a['Z_FST_15']=(a['MEAN_FST_15']-np.nanmean(a['MEAN_FST_15']))/np.nanstd(a['MEAN_FST_15'])
a['Z_FST_16']=(a['MEAN_FST_16']-np.nanmean(a['MEAN_FST_16']))/np.nanstd(a['MEAN_FST_16'])
a['Z_FST_17']=(a['MEAN_FST_17']-np.nanmean(a['MEAN_FST_17']))/np.nanstd(a['MEAN_FST_17'])


z['Z_FST_0']=(z['MEAN_FST_0']-np.nanmean(z['MEAN_FST_0']))/np.nanstd(z['MEAN_FST_0'])
z['Z_FST_1']=(z['MEAN_FST_1']-np.nanmean(z['MEAN_FST_1']))/np.nanstd(z['MEAN_FST_1'])
z['Z_FST_2']=(z['MEAN_FST_2']-np.nanmean(z['MEAN_FST_2']))/np.nanstd(z['MEAN_FST_2'])
z['Z_FST_3']=(z['MEAN_FST_3']-np.nanmean(z['MEAN_FST_3']))/np.nanstd(z['MEAN_FST_3'])
z['Z_FST_4']=(z['MEAN_FST_4']-np.nanmean(z['MEAN_FST_4']))/np.nanstd(z['MEAN_FST_4'])
z['Z_FST_5']=(z['MEAN_FST_5']-np.nanmean(z['MEAN_FST_5']))/np.nanstd(z['MEAN_FST_5'])
z['Z_FST_6']=(z['MEAN_FST_6']-np.nanmean(z['MEAN_FST_6']))/np.nanstd(z['MEAN_FST_6'])
z['Z_FST_7']=(z['MEAN_FST_7']-np.nanmean(z['MEAN_FST_7']))/np.nanstd(z['MEAN_FST_7'])
z['Z_FST_8']=(z['MEAN_FST_8']-np.nanmean(z['MEAN_FST_8']))/np.nanstd(z['MEAN_FST_8'])
z['Z_FST_9']=(z['MEAN_FST_9']-np.nanmean(z['MEAN_FST_9']))/np.nanstd(z['MEAN_FST_9'])
z['Z_FST_10']=(z['MEAN_FST_10']-np.nanmean(z['MEAN_FST_10']))/np.nanstd(z['MEAN_FST_10'])
z['Z_FST_11']=(z['MEAN_FST_11']-np.nanmean(z['MEAN_FST_11']))/np.nanstd(z['MEAN_FST_11'])
z['Z_FST_12']=(z['MEAN_FST_12']-np.nanmean(z['MEAN_FST_12']))/np.nanstd(z['MEAN_FST_12'])
z['Z_FST_13']=(z['MEAN_FST_13']-np.nanmean(z['MEAN_FST_13']))/np.nanstd(z['MEAN_FST_13'])
z['Z_FST_14']=(z['MEAN_FST_14']-np.nanmean(z['MEAN_FST_14']))/np.nanstd(z['MEAN_FST_14'])
z['Z_FST_15']=(z['MEAN_FST_15']-np.nanmean(z['MEAN_FST_15']))/np.nanstd(z['MEAN_FST_15'])
z['Z_FST_16']=(z['MEAN_FST_16']-np.nanmean(z['MEAN_FST_16']))/np.nanstd(z['MEAN_FST_16'])
z['Z_FST_17']=(z['MEAN_FST_17']-np.nanmean(z['MEAN_FST_17']))/np.nanstd(z['MEAN_FST_17'])

az=a[['CHROM','BIN_START','BIN_END','Z_FST_0','Z_FST_1','Z_FST_2','Z_FST_3','Z_FST_4','Z_FST_5','Z_FST_6','Z_FST_7','Z_FST_8','Z_FST_9','Z_FST_10','Z_FST_11','Z_FST_12','Z_FST_13','Z_FST_14','Z_FST_15','Z_FST_16','Z_FST_17']]
zz=z[['CHROM','BIN_START','BIN_END','Z_FST_0','Z_FST_1','Z_FST_2','Z_FST_3','Z_FST_4','Z_FST_5','Z_FST_6','Z_FST_7','Z_FST_8','Z_FST_9','Z_FST_10','Z_FST_11','Z_FST_12','Z_FST_13','Z_FST_14','Z_FST_15','Z_FST_16','Z_FST_17']]


aazz=pandas.concat([az,zz], ignore_index=True)

data2=personal_popgen.join_raw_data_base([data2,aazz])


#####################################Ztransform TajimaD####################################
a['juvernica_TajimaD_Z_trans']=(a['juvernica_TajimaD']-np.nanmean(a['juvernica_TajimaD']))/np.nanstd(a['juvernica_TajimaD'])
a['kazak_juvernica_TajimaD_Z_trans']=(a['kazak_juvernica_TajimaD']-np.nanmean(a['kazak_juvernica_TajimaD']))/np.nanstd(a['kazak_juvernica_TajimaD'])
a['kaz_sin_TajimaD_Z_trans']=(a['kaz_sin_TajimaD']-np.nanmean(a['kaz_sin_TajimaD']))/np.nanstd(a['kaz_sin_TajimaD'])
a['reali_TajimaD_Z_trans']=(a['reali_TajimaD']-np.nanmean(a['reali_TajimaD']))/np.nanstd(a['reali_TajimaD'])
a['sinapis_TajimaD_Z_trans']=(a['sinapis_TajimaD']-np.nanmean(a['sinapis_TajimaD']))/np.nanstd(a['sinapis_TajimaD'])
a['spanish_reali_TajimaD_Z_trans']=(a['spanish_reali_TajimaD']-np.nanmean(a['spanish_reali_TajimaD']))/np.nanstd(a['spanish_reali_TajimaD'])
a['spanish_sinapis_TajimaD_Z_trans']=(a['spanish_sinapis_TajimaD']-np.nanmean(a['spanish_sinapis_TajimaD']))/np.nanstd(a['spanish_sinapis_TajimaD'])
a['swe_sin_TajimaD_Z_trans']=(a['swe_sin_TajimaD']-np.nanmean(a['swe_sin_TajimaD']))/np.nanstd(a['swe_sin_TajimaD'])

z['juvernica_TajimaD_Z_trans']=(z['juvernica_TajimaD']-np.nanmean(z['juvernica_TajimaD']))/np.nanstd(z['juvernica_TajimaD'])
z['kazak_juvernica_TajimaD_Z_trans']=(z['kazak_juvernica_TajimaD']-np.nanmean(z['kazak_juvernica_TajimaD']))/np.nanstd(z['kazak_juvernica_TajimaD'])
z['kaz_sin_TajimaD_Z_trans']=(z['kaz_sin_TajimaD']-np.nanmean(z['kaz_sin_TajimaD']))/np.nanstd(z['kaz_sin_TajimaD'])
z['reali_TajimaD_Z_trans']=(z['reali_TajimaD']-np.nanmean(z['reali_TajimaD']))/np.nanstd(z['reali_TajimaD'])
z['sinapis_TajimaD_Z_trans']=(z['sinapis_TajimaD']-np.nanmean(z['sinapis_TajimaD']))/np.nanstd(z['sinapis_TajimaD'])
z['spanish_reali_TajimaD_Z_trans']=(z['spanish_reali_TajimaD']-np.nanmean(z['spanish_reali_TajimaD']))/np.nanstd(z['spanish_reali_TajimaD'])
z['spanish_sinapis_TajimaD_Z_trans']=(z['spanish_sinapis_TajimaD']-np.nanmean(z['spanish_sinapis_TajimaD']))/np.nanstd(z['spanish_sinapis_TajimaD'])
z['swe_sin_TajimaD_Z_trans']=(z['swe_sin_TajimaD']-np.nanmean(z['swe_sin_TajimaD']))/np.nanstd(z['swe_sin_TajimaD'])

az=a[['CHROM','BIN_START','BIN_END','juvernica_TajimaD_Z_trans','kazak_juvernica_TajimaD_Z_trans','kaz_sin_TajimaD_Z_trans','reali_TajimaD_Z_trans','sinapis_TajimaD_Z_trans','spanish_reali_TajimaD_Z_trans','spanish_sinapis_TajimaD_Z_trans','swe_sin_TajimaD_Z_trans']]
zz=z[['CHROM','BIN_START','BIN_END','juvernica_TajimaD_Z_trans','kazak_juvernica_TajimaD_Z_trans','kaz_sin_TajimaD_Z_trans','reali_TajimaD_Z_trans','sinapis_TajimaD_Z_trans','spanish_reali_TajimaD_Z_trans','spanish_sinapis_TajimaD_Z_trans','swe_sin_TajimaD_Z_trans']]

aazz=pandas.concat([a,z], ignore_index=True)
data2=personal_popgen.join_raw_data_base([data2,aazz])


################outliers detections###################################
all_a_Z_FST_0 =a[a['Z_FST_0']  > a.Z_FST_0.dropna().quantile(0.99)][["kaz_sin_my_PI","swe_sin_allele_my_PI",'DXY_0','CHROM','BIN_START','BIN_END','auto_allo','kaz_sin_swe_sin_allele_Avg_Pi','Z_FST_0','MEAN_FST_0']]
all_a_Z_FST_1 =a[a['Z_FST_1']  > a.Z_FST_1.dropna().quantile(0.99)][["spanish_sinapis_my_PI","swe_sin_allele_my_PI",'DXY_1','CHROM','BIN_START','BIN_END','auto_allo', "spanish_sinapis_swe_sin_allele_Avg_Pi",'Z_FST_1','MEAN_FST_1']]
all_a_Z_FST_2 =a[a['Z_FST_2']  > a.Z_FST_2.dropna().quantile(0.99)][["kaz_sin_my_PI","spanish_sinapis_my_PI",'DXY_2','CHROM','BIN_START','BIN_END','auto_allo', "kaz_sin_spanish_sinapis_Avg_Pi",'Z_FST_2','MEAN_FST_2']]
all_a_Z_FST_3 =a[a['Z_FST_3']  > a.Z_FST_3.dropna().quantile(0.99)][["irish_juvernica_my_PI","kazak_juvernica_my_PI",'DXY_3','CHROM','BIN_START','BIN_END','auto_allo', "irish_juvernica_kazak_juvernica_Avg_Pi",'Z_FST_3','MEAN_FST_3']]
all_a_Z_FST_15=a[a['Z_FST_15'] > a.Z_FST_15.dropna().quantile(0.99)][["sinapis_my_PI","reali_my_PI",'DXY_15','CHROM','BIN_START','BIN_END','auto_allo',"sinapis_reali_Avg_Pi",'Z_FST_15','MEAN_FST_15']]
all_a_Z_FST_16=a[a['Z_FST_16'] > a.Z_FST_16.dropna().quantile(0.99)][["reali_my_PI","juvernica_my_PI",'DXY_16','CHROM','BIN_START','BIN_END','auto_allo',"reali_juvernica_Avg_Pi",'Z_FST_16','MEAN_FST_16']]
all_a_Z_FST_17=a[a['Z_FST_17'] > a.Z_FST_17.dropna().quantile(0.99)][["sinapis_my_PI","juvernica_my_PI",'DXY_17','CHROM','BIN_START','BIN_END','auto_allo',"sinapis_juvernica_Avg_Pi",'Z_FST_17','MEAN_FST_17']]
all_z_Z_FST_0 =z[z['Z_FST_0']  > z.Z_FST_0.dropna().quantile(0.99))][["kaz_sin_my_PI","swe_sin_allele_my_PI",'DXY_0','CHROM','BIN_START','BIN_END','auto_allo','kaz_sin_swe_sin_allele_Avg_Pi','Z_FST_0','MEAN_FST_0']]
all_z_Z_FST_1 =z[z['Z_FST_1']  > z.Z_FST_1.dropna().quantile(0.99))][["spanish_sinapis_my_PI","swe_sin_allele_my_PI",'DXY_1','CHROM','BIN_START','BIN_END','auto_allo', "spanish_sinapis_swe_sin_allele_Avg_Pi",'Z_FST_1','MEAN_FST_1']]
all_z_Z_FST_2 =z[z['Z_FST_2']  > z.Z_FST_2.dropna().quantile(0.99))][["kaz_sin_my_PI","spanish_sinapis_my_PI",'DXY_2','CHROM','BIN_START','BIN_END','auto_allo', "kaz_sin_spanish_sinapis_Avg_Pi",'Z_FST_2','MEAN_FST_2']]
all_z_Z_FST_3 =z[z['Z_FST_3']  > z.Z_FST_3.dropna().quantile(0.99))][["irish_juvernica_my_PI","kazak_juvernica_my_PI",'DXY_3','CHROM','BIN_START','BIN_END','auto_allo', "irish_juvernica_kazak_juvernica_Avg_Pi",'Z_FST_3','MEAN_FST_3']]
all_z_Z_FST_15=z[z['Z_FST_15'] > z.Z_FST_15.dropna().quantile(0.99))][["sinapis_my_PI","reali_my_PI",'DXY_15','CHROM','BIN_START','BIN_END','auto_allo',"sinapis_reali_Avg_Pi",'Z_FST_15','MEAN_FST_15']]
all_z_Z_FST_16=z[z['Z_FST_16'] > z.Z_FST_16.dropna().quantile(0.99))][["reali_my_PI","juvernica_my_PI",'DXY_16','CHROM','BIN_START','BIN_END','auto_allo',"reali_juvernica_Avg_Pi",'Z_FST_16','MEAN_FST_16']]
all_z_Z_FST_17=z[z['Z_FST_17'] > z.Z_FST_17.dropna().quantile(0.99))][["sinapis_my_PI","juvernica_my_PI",'DXY_17','CHROM','BIN_START','BIN_END','auto_allo',"sinapis_juvernica_Avg_Pi",'Z_FST_17','MEAN_FST_17']]



a_Z_FST_0 =a[(a['Z_FST_0']  > a.Z_FST_0.dropna().quantile(0.99)) &  (a['DXY_0']  > a.DXY_0.dropna().mean())][["kaz_sin_my_PI","swe_sin_allele_my_PI",'DXY_0','CHROM','BIN_START','BIN_END','auto_allo','kaz_sin_swe_sin_allele_Avg_Pi','Z_FST_0','MEAN_FST_0']]
a_Z_FST_1 =a[(a['Z_FST_1']  > a.Z_FST_1.dropna().quantile(0.99)) &  (a['DXY_1']  > a.DXY_1.dropna().mean())][["spanish_sinapis_my_PI","swe_sin_allele_my_PI",'DXY_1','CHROM','BIN_START','BIN_END','auto_allo', "spanish_sinapis_swe_sin_allele_Avg_Pi",'Z_FST_1','MEAN_FST_1']]
a_Z_FST_2 =a[(a['Z_FST_2']  > a.Z_FST_2.dropna().quantile(0.99)) &  (a['DXY_2']  > a.DXY_2.dropna().mean())][["kaz_sin_my_PI","spanish_sinapis_my_PI",'DXY_2','CHROM','BIN_START','BIN_END','auto_allo', "kaz_sin_spanish_sinapis_Avg_Pi",'Z_FST_2','MEAN_FST_2']]
a_Z_FST_3 =a[(a['Z_FST_3']  > a.Z_FST_3.dropna().quantile(0.99)) &  (a['DXY_3']  > a.DXY_3.dropna().mean())][["irish_juvernica_my_PI","kazak_juvernica_my_PI",'DXY_3','CHROM','BIN_START','BIN_END','auto_allo', "irish_juvernica_kazak_juvernica_Avg_Pi",'Z_FST_3','MEAN_FST_3']]
a_Z_FST_15=a[(a['Z_FST_15'] > a.Z_FST_15.dropna().quantile(0.99)) & (a['DXY_15'] > a.DXY_15.dropna().mean())][["sinapis_my_PI","reali_my_PI",'DXY_15','CHROM','BIN_START','BIN_END','auto_allo',"sinapis_reali_Avg_Pi",'Z_FST_15','MEAN_FST_15']]
a_Z_FST_16=a[(a['Z_FST_16'] > a.Z_FST_16.dropna().quantile(0.99)) & (a['DXY_16'] > a.DXY_16.dropna().mean())][["reali_my_PI","juvernica_my_PI",'DXY_16','CHROM','BIN_START','BIN_END','auto_allo',"reali_juvernica_Avg_Pi",'Z_FST_16','MEAN_FST_16']]
a_Z_FST_17=a[(a['Z_FST_17'] > a.Z_FST_17.dropna().quantile(0.99)) & (a['DXY_17'] > a.DXY_17.dropna().mean())][["sinapis_my_PI","juvernica_my_PI",'DXY_17','CHROM','BIN_START','BIN_END','auto_allo',"sinapis_juvernica_Avg_Pi",'Z_FST_17','MEAN_FST_17']]

z_Z_FST_0 =z[(z['Z_FST_0']  > z.Z_FST_0.dropna().quantile(0.99)) &   (z['DXY_0']  > z.DXY_0.dropna().mean())][["kaz_sin_my_PI","swe_sin_allele_my_PI",'DXY_0','CHROM','BIN_START','BIN_END','auto_allo','kaz_sin_swe_sin_allele_Avg_Pi','Z_FST_0','MEAN_FST_0']]
z_Z_FST_1 =z[(z['Z_FST_1']  > z.Z_FST_1.dropna().quantile(0.99)) &   (z['DXY_1']  > z.DXY_1.dropna().mean())][["spanish_sinapis_my_PI","swe_sin_allele_my_PI",'DXY_1','CHROM','BIN_START','BIN_END','auto_allo', "spanish_sinapis_swe_sin_allele_Avg_Pi",'Z_FST_1','MEAN_FST_1']]
z_Z_FST_2 =z[(z['Z_FST_2']  > z.Z_FST_2.dropna().quantile(0.99)) &   (z['DXY_2']  > z.DXY_2.dropna().mean())][["kaz_sin_my_PI","spanish_sinapis_my_PI",'DXY_2','CHROM','BIN_START','BIN_END','auto_allo', "kaz_sin_spanish_sinapis_Avg_Pi",'Z_FST_2','MEAN_FST_2']]
z_Z_FST_3 =z[(z['Z_FST_3']  > z.Z_FST_3.dropna().quantile(0.99)) &   (z['DXY_3']  > z.DXY_3.dropna().mean())][["irish_juvernica_my_PI","kazak_juvernica_my_PI",'DXY_3','CHROM','BIN_START','BIN_END','auto_allo', "irish_juvernica_kazak_juvernica_Avg_Pi",'Z_FST_3','MEAN_FST_3']]
z_Z_FST_15=z[(z['Z_FST_15'] > z.Z_FST_15.dropna().quantile(0.99))&   (z['DXY_15'] > z.DXY_15.dropna().mean())][["sinapis_my_PI","reali_my_PI",'DXY_15','CHROM','BIN_START','BIN_END','auto_allo',"sinapis_reali_Avg_Pi",'Z_FST_15','MEAN_FST_15']]
z_Z_FST_16=z[(z['Z_FST_16'] > z.Z_FST_16.dropna().quantile(0.99))&   (z['DXY_16'] > z.DXY_16.dropna().mean())][["reali_my_PI","juvernica_my_PI",'DXY_16','CHROM','BIN_START','BIN_END','auto_allo',"reali_juvernica_Avg_Pi",'Z_FST_16','MEAN_FST_16']]
z_Z_FST_17=z[(z['Z_FST_17'] > z.Z_FST_17.dropna().quantile(0.99))&   (z['DXY_17'] > z.DXY_17.dropna().mean())][["sinapis_my_PI","juvernica_my_PI",'DXY_17','CHROM','BIN_START','BIN_END','auto_allo',"sinapis_juvernica_Avg_Pi",'Z_FST_17','MEAN_FST_17']]

data_Z_FST_0 =pandas.concat([a_Z_FST_0 ,z_Z_FST_0 ], ignore_index=True)
data_Z_FST_1 =pandas.concat([a_Z_FST_1 ,z_Z_FST_1 ], ignore_index=True)
data_Z_FST_2 =pandas.concat([a_Z_FST_2 ,z_Z_FST_2 ], ignore_index=True)
data_Z_FST_3 =pandas.concat([a_Z_FST_3 ,z_Z_FST_3 ], ignore_index=True)
data_Z_FST_15=pandas.concat([a_Z_FST_15,z_Z_FST_15], ignore_index=True)
data_Z_FST_16=pandas.concat([a_Z_FST_16,z_Z_FST_16], ignore_index=True)
data_Z_FST_17=pandas.concat([a_Z_FST_17,z_Z_FST_17], ignore_index=True)






data_Z_FST_0.to_csv('Z_FST_0_all_windows')
data_Z_FST_1.to_csv('Z_FST_1_all_windows')
data_Z_FST_2.to_csv('Z_FST_2_all_windows')
data_Z_FST_3.to_csv('Z_FST_3_all_windows')
data_Z_FST_15.to_csv('Z_FST_15_all_windows')
data_Z_FST_16.to_csv('Z_FST_16_all_windows')
data_Z_FST_17.to_csv('Z_FST_17_all_windows')








data_Z_FST_0 =data_Z_FST_0[["kaz_sin_my_PI","swe_sin_allele_my_PI",'DXY_0','CHROM','BIN_START','BIN_END','auto_allo','MEAN_FST_0',"kaz_sin_swe_sin_allele_Avg_Pi"]]
data_Z_FST_1 =data_Z_FST_1[["spanish_sinapis_my_PI","swe_sin_allele_my_PI",'DXY_1','CHROM','BIN_START','BIN_END','auto_allo','MEAN_FST_1',"spanish_sinapis_swe_sin_allele_Avg_Pi"]]
data_Z_FST_2 =data_Z_FST_2[["kaz_sin_my_PI","spanish_sinapis_my_PI",'DXY_2','CHROM','BIN_START','BIN_END','auto_allo','MEAN_FST_2',"kaz_sin_spanish_sinapis_Avg_Pi"]]
data_Z_FST_3 =data_Z_FST_3[["irish_juvernica_my_PI","kazak_juvernica_my_PI",'DXY_3','CHROM','BIN_START','BIN_END','auto_allo','MEAN_FST_3',"irish_juvernica_kazak_juvernica_Avg_Pi"]]
data_Z_FST_15=data_Z_FST_15[["sinapis_my_PI","reali_my_PI",'DXY_15','CHROM','BIN_START','BIN_END','auto_allo','MEAN_FST_15',"sinapis_reali_Avg_Pi"]]
data_Z_FST_16=data_Z_FST_16[["reali_my_PI","juvernica_my_PI",'DXY_16','CHROM','BIN_START','BIN_END','auto_allo','MEAN_FST_16',"reali_juvernica_Avg_Pi"]]
data_Z_FST_17=data_Z_FST_17[["sinapis_my_PI","juvernica_my_PI",'DXY_17','CHROM','BIN_START','BIN_END','auto_allo','MEAN_FST_17',"sinapis_juvernica_Avg_Pi"]]




print len(data_Z_FST_0)
print len(data_Z_FST_1)
print len(data_Z_FST_2)
print len(data_Z_FST_3)
print len(data_Z_FST_15)
print len(data_Z_FST_16)
print len(data_Z_FST_17)



####################################
 
a_juvernica_TajimaD_outliers=a[a['juvernica_TajimaD_Z_trans']  > a.juvernica_TajimaD_Z_trans.dropna().quantile(0.99)][['CHROM','BIN_START','BIN_END','auto_allo']]
a_kazak_juvernica_TajimaD_outliers=a[a['kazak_juvernica_TajimaD_Z_trans']  > a.kazak_juvernica_TajimaD_Z_trans.dropna().quantile(0.99)][['CHROM','BIN_START','BIN_END','auto_allo']]
a_kaz_sin_TajimaD_outliers=a[a['kaz_sin_TajimaD_Z_trans']  > a.kaz_sin_TajimaD_Z_trans.dropna().quantile(0.99)][['CHROM','BIN_START','BIN_END','auto_allo']]
a_reali_TajimaD_outliers=a[a['reali_TajimaD_Z_trans']  > a.reali_TajimaD_Z_trans.dropna().quantile(0.99)][['CHROM','BIN_START','BIN_END','auto_allo']]
a_sinapis_TajimaD_outliers=a[a['sinapis_TajimaD_Z_trans']  > a.sinapis_TajimaD_Z_trans.dropna().quantile(0.99)][['CHROM','BIN_START','BIN_END','auto_allo']]
a_spanish_reali_TajimaD_outliers=a[a['spanish_reali_TajimaD_Z_trans']  > a.spanish_reali_TajimaD_Z_trans.dropna().quantile(0.99)][['CHROM','BIN_START','BIN_END','auto_allo']]
a_spanish_sinapis_TajimaD_outliers=a[a['spanish_sinapis_TajimaD_Z_trans']  > a.spanish_sinapis_TajimaD_Z_trans.dropna().quantile(0.99)][['CHROM','BIN_START','BIN_END','auto_allo']]
a_swe_sin_TajimaD_outliers=a[a['swe_sin_TajimaD_Z_trans']  > a.swe_sin_TajimaD_Z_trans.dropna().quantile(0.99)][['CHROM','BIN_START','BIN_END','auto_allo']]

z_juvernica_TajimaD_outliers=z[z['juvernica_TajimaD_Z_trans']  > z.juvernica_TajimaD_Z_trans.dropna().quantile(0.99)][['CHROM','BIN_START','BIN_END','auto_allo']]
z_kazak_juvernica_TajimaD_outliers=z[z['kazak_juvernica_TajimaD_Z_trans']  > z.kazak_juvernica_TajimaD_Z_trans.dropna().quantile(0.99)][['CHROM','BIN_START','BIN_END','auto_allo']]
z_kaz_sin_TajimaD_outliers=z[z['kaz_sin_TajimaD_Z_trans']  > z.kaz_sin_TajimaD_Z_trans.dropna().quantile(0.99)][['CHROM','BIN_START','BIN_END','auto_allo']]
z_reali_TajimaD_outliers=z[z['reali_TajimaD_Z_trans']  > z.reali_TajimaD_Z_trans.dropna().quantile(0.99)][['CHROM','BIN_START','BIN_END','auto_allo']]
z_sinapis_TajimaD_outliers=z[z['sinapis_TajimaD_Z_trans']  > z.sinapis_TajimaD_Z_trans.dropna().quantile(0.99)][['CHROM','BIN_START','BIN_END','auto_allo']]
z_spanish_reali_TajimaD_outliers=z[z['spanish_reali_TajimaD_Z_trans']  > z.spanish_reali_TajimaD_Z_trans.dropna().quantile(0.99)][['CHROM','BIN_START','BIN_END','auto_allo']]
z_spanish_sinapis_TajimaD_outliers=z[z['spanish_sinapis_TajimaD_Z_trans']  > z.spanish_sinapis_TajimaD_Z_trans.dropna().quantile(0.99)][['CHROM','BIN_START','BIN_END','auto_allo']]
z_swe_sin_TajimaD_outliers=z[z['swe_sin_TajimaD_Z_trans']  > z.swe_sin_TajimaD_Z_trans.dropna().quantile(0.99)][['CHROM','BIN_START','BIN_END','auto_allo']]



data_Z_FST_juvernica=pandas.concat([a_juvernica_TajimaD_outliers, z_juvernica_TajimaD_outliers], ignore_index=True)
data_Z_FST_kazak_juvernica=pandas.concat([a_kazak_juvernica_TajimaD_outliers, z_kazak_juvernica_TajimaD_outliers], ignore_index=True)
data_Z_FST_kaz_sin=pandas.concat([a_kaz_sin_TajimaD_outliers, z_kaz_sin_TajimaD_outliers], ignore_index=True)
data_Z_FST_reali=pandas.concat([a_reali_TajimaD_outliers, z_reali_TajimaD_outliers], ignore_index=True)
data_Z_FST_sinapis=pandas.concat([a_sinapis_TajimaD_outliers, z_sinapis_TajimaD_outliers], ignore_index=True)
data_Z_FST_spanish_reali=pandas.concat([a_spanish_reali_TajimaD_outliers, z_spanish_reali_TajimaD_outliers], ignore_index=True)
data_Z_FST_spanish_sinapis=pandas.concat([a_spanish_sinapis_TajimaD_outliers, z_spanish_sinapis_TajimaD_outliers], ignore_index=True)
data_Z_FST_swe_sin=pandas.concat([a_swe_sin_TajimaD_outliers, z_swe_sin_TajimaD_outliers], ignore_index=True)






def genes_in_outliers(data):
        head_stuff=['scaffold','source','feature','start','end','.','orentation','..','infor']
        annotation=pandas.read_table("/proj/uppstore2017185/b2014034/NBIS_annotation_leptidea/gff/gene-builds/leptidea_sinapis_rc1.gff",skiprows=1, names=head_stuff)
        annotation=annotation[['scaffold','source','feature','start','end','infor','orentation']]
        annotation=annotation[annotation.feature=='gene']
        list_genes=[]
        for window_index, window_1 in data.iterrows():
                window=annotation[(annotation['scaffold'] == window_1.CHROM)]
                window=window[(window['start']>=int(window_1.BIN_START)) & (window['end']<=int(window_1.BIN_END))]
                window=window.append(window[(window['start']<=int(window_1.BIN_START)) & (window['end']>=int(window_1.BIN_END))])
                #print window.empty
                if window.empty != 'False':
                        #print window
                        #print window_1.BIN_START
                        #print window_1.BIN_END
                        for i_index, i in window.iterrows():
                                list_genes.append(str(i['infor']))
        list_genes_=[i.split(';')[1][5:] for i in list_genes]
        #f=open('data_Z_FST_16.fasta', 'w+')
        file_name='ZFST_'+list(data)[2].split('_')[1]+'_genes_in_outliers.csv'
        d = {'a': list_genes_}
        df = pandas.DataFrame(data=d)
        df.to_csv(file_name, header=False, index=False,sep = " ")
        return list_genes




####annotation loaded ###############
head_stuff=['scaffold','source','feature','start','end','.','orentation','..','infor']
annotation=pandas.read_table("/proj/uppstore2017185/b2014034/NBIS_annotation_leptidea/gff/gene-builds/leptidea_sinapis_rc1.gff",skiprows=1, names=head_stuff)
annotation_gene=annotation[(annotation.feature=='mRNA')]
annotation_gene['Transcript']=annotation_gene.infor.str[3:31]
annotation_gene['Gene']=annotation_gene.infor.str[39:67]
annotation_gene=annotation_gene[['Transcript','Gene']]

DM_genes=pandas.read_csv("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/gene_ortholog_flybase/overlap_new.csv", names=['DM_gene','LS_gene'], sep=' ')


Z_genes_Z_FST_0 = [str(annotation_gene[(annotation_gene.Gene==i.split(';')[0][3:])]['Transcript']).split()[1] for i in genes_in_outliers(data_Z_FST_0[(data_Z_FST_0.auto_allo=='Z')]) ]
Z_genes_Z_FST_1 = [str(annotation_gene[(annotation_gene.Gene==i.split(';')[0][3:])]['Transcript']).split()[1] for i in genes_in_outliers(data_Z_FST_1[(data_Z_FST_1.auto_allo=='Z')]) ]
Z_genes_Z_FST_2 = [str(annotation_gene[(annotation_gene.Gene==i.split(';')[0][3:])]['Transcript']).split()[1] for i in genes_in_outliers(data_Z_FST_2[(data_Z_FST_2.auto_allo=='Z')]) ]
Z_genes_Z_FST_3 = [str(annotation_gene[(annotation_gene.Gene==i.split(';')[0][3:])]['Transcript']).split()[1] for i in genes_in_outliers(data_Z_FST_3[(data_Z_FST_3.auto_allo=='Z')]) ]
Z_genes_Z_FST_15= [str(annotation_gene[(annotation_gene.Gene==i.split(';')[0][3:])]['Transcript']).split()[1] for i in genes_in_outliers(data_Z_FST_15[(data_Z_FST_15.auto_allo=='Z')]) ]
Z_genes_Z_FST_16= [str(annotation_gene[(annotation_gene.Gene==i.split(';')[0][3:])]['Transcript']).split()[1] for i in genes_in_outliers(data_Z_FST_16[(data_Z_FST_16.auto_allo=='Z')]) ]
Z_genes_Z_FST_17= [str(annotation_gene[(annotation_gene.Gene==i.split(';')[0][3:])]['Transcript']).split()[1] for i in genes_in_outliers(data_Z_FST_17[(data_Z_FST_17.auto_allo=='Z')]) ]

A_genes_Z_FST_0 = [str(annotation_gene[(annotation_gene.Gene==i.split(';')[0][3:])]['Transcript']).split()[1] for i in genes_in_outliers(data_Z_FST_0[(data_Z_FST_0.auto_allo=='auto')]) ]
A_genes_Z_FST_1 = [str(annotation_gene[(annotation_gene.Gene==i.split(';')[0][3:])]['Transcript']).split()[1] for i in genes_in_outliers(data_Z_FST_1[(data_Z_FST_1.auto_allo=='auto')]) ]
A_genes_Z_FST_2 = [str(annotation_gene[(annotation_gene.Gene==i.split(';')[0][3:])]['Transcript']).split()[1] for i in genes_in_outliers(data_Z_FST_2[(data_Z_FST_2.auto_allo=='auto')]) ]
A_genes_Z_FST_3 = [str(annotation_gene[(annotation_gene.Gene==i.split(';')[0][3:])]['Transcript']).split()[1] for i in genes_in_outliers(data_Z_FST_3[(data_Z_FST_3.auto_allo=='auto')]) ]
A_genes_Z_FST_15= [str(annotation_gene[(annotation_gene.Gene==i.split(';')[0][3:])]['Transcript']).split()[1] for i in genes_in_outliers(data_Z_FST_15[(data_Z_FST_15.auto_allo=='auto')]) ]
A_genes_Z_FST_16= [str(annotation_gene[(annotation_gene.Gene==i.split(';')[0][3:])]['Transcript']).split()[1] for i in genes_in_outliers(data_Z_FST_16[(data_Z_FST_16.auto_allo=='auto')]) ]
A_genes_Z_FST_17= [str(annotation_gene[(annotation_gene.Gene==i.split(';')[0][3:])]['Transcript']).split()[1] for i in genes_in_outliers(data_Z_FST_17[(data_Z_FST_17.auto_allo=='auto')]) ]


genes_Z_FST_0 = [str(annotation_gene[(annotation_gene.Gene==i.split(';')[0][3:])]['Transcript']).split()[1] for i in genes_in_outliers(data_Z_FST_0) ]
genes_Z_FST_1 = [str(annotation_gene[(annotation_gene.Gene==i.split(';')[0][3:])]['Transcript']).split()[1] for i in genes_in_outliers(data_Z_FST_1) ]
genes_Z_FST_2 = [str(annotation_gene[(annotation_gene.Gene==i.split(';')[0][3:])]['Transcript']).split()[1] for i in genes_in_outliers(data_Z_FST_2) ]
genes_Z_FST_3 = [str(annotation_gene[(annotation_gene.Gene==i.split(';')[0][3:])]['Transcript']).split()[1] for i in genes_in_outliers(data_Z_FST_3) ]
genes_Z_FST_15= [str(annotation_gene[(annotation_gene.Gene==i.split(';')[0][3:])]['Transcript']).split()[1] for i in genes_in_outliers(data_Z_FST_15) ]
genes_Z_FST_16= [str(annotation_gene[(annotation_gene.Gene==i.split(';')[0][3:])]['Transcript']).split()[1] for i in genes_in_outliers(data_Z_FST_16) ]
genes_Z_FST_17= [str(annotation_gene[(annotation_gene.Gene==i.split(';')[0][3:])]['Transcript']).split()[1] for i in genes_in_outliers(data_Z_FST_17) ]

genes_Z_FST_0=DM_genes[(DM_genes.LS_gene.isin(genes_Z_FST_0))]['DM_gene']
genes_Z_FST_1=DM_genes[(DM_genes.LS_gene.isin(genes_Z_FST_1))]['DM_gene']
genes_Z_FST_2=DM_genes[(DM_genes.LS_gene.isin(genes_Z_FST_2))]['DM_gene']
genes_Z_FST_3=DM_genes[(DM_genes.LS_gene.isin(genes_Z_FST_3))]['DM_gene']
genes_Z_FST_15=DM_genes[(DM_genes.LS_gene.isin(genes_Z_FST_15))]['DM_gene']
genes_Z_FST_16=DM_genes[(DM_genes.LS_gene.isin(genes_Z_FST_16))]['DM_gene']
genes_Z_FST_17=DM_genes[(DM_genes.LS_gene.isin(genes_Z_FST_17))]['DM_gene']




##########################get specific genes########################

specific_0=set(genes_Z_FST_0).difference(list(set(genes_Z_FST_1+genes_Z_FST_2)))
specific_1=set(genes_Z_FST_1).difference(list(set(genes_Z_FST_0+genes_Z_FST_2)))
specific_2=set(genes_Z_FST_2).difference(list(set(genes_Z_FST_0+genes_Z_FST_1)))



specific_15=set(genes_Z_FST_15).difference(list(set(genes_Z_FST_16+genes_Z_FST_17)))
specific_16=set(genes_Z_FST_16).difference(list(set(genes_Z_FST_15+genes_Z_FST_17)))
specific_17=set(genes_Z_FST_17).difference(list(set(genes_Z_FST_16+genes_Z_FST_15)))

##########################get specific genes########################
#print corr_pvalue(data[['DXY_0','kaz_sin_swe_sin_allele_Avg_Pi']])
#print corr_pvalue(data[['DXY_1','spanish_sinapis_swe_sin_allele_Avg_Pi']])
#print corr_pvalue(data[['DXY_2','kaz_sin_spanish_sinapis_Avg_Pi']])
#print corr_pvalue(data[['DXY_3','irish_juvernica_kazak_juvernica_Avg_Pi']])
#print corr_pvalue(data[['DXY_15','sinapis_reali_Avg_Pi']])
#print corr_pvalue(data[['DXY_16','reali_juvernica_Avg_Pi']])
#print corr_pvalue(data[['DXY_17','sinapis_juvernica_Avg_Pi']])
#
#
#print corr_pvalue(data[['DXY_0' ,'kaz_sin_swe_sin_allele_Avg_Pi','MEAN_FST_0', 'kaz_sin_Avg_COV', 'swe_sin_allele_Avg_COV' ]])
#print corr_pvalue(data[['DXY_1' ,'spanish_sinapis_swe_sin_allele_Avg_Pi','MEAN_FST_1', 'spanish_sinapis_Avg_COV', 'swe_sin_allele_Avg_COV' ]])
#print corr_pvalue(data[['DXY_2' ,'kaz_sin_spanish_sinapis_Avg_Pi','MEAN_FST_2', 'kaz_sin_Avg_COV', 'spanish_sinapis_Avg_COV' ]])
#print corr_pvalue(data[['DXY_3' ,'irish_juvernica_kazak_juvernica_Avg_Pi','MEAN_FST_3', 'irish_juvernica_Avg_COV', 'kazak_juvernica_Avg_COV' ]])
#print corr_pvalue(data[['DXY_15','sinapis_reali_Avg_Pi','MEAN_FST_15', 'sinapis_Avg_COV','reali_Avg_COV']])
#print corr_pvalue(data[['DXY_16','reali_juvernica_Avg_Pi','MEAN_FST_16', 'reali_Avg_COV', 'juvernica_Avg_COV']])
#print corr_pvalue(data[['DXY_17','sinapis_juvernica_Avg_Pi','MEAN_FST_17', 'sinapis_Avg_COV','juvernica_Avg_COV']])











################### some plots ##################################
sns.boxplot(sites_filter_0_pi,showfliers=False,palette=pi_colours).set(ylim=(0,0.0105))
sns.boxplot(sites_filter_0_dxy,showfliers=False).set(ylim=(0,0.5))
sns.boxplot(sites_filter_0_fst,showfliers=False).set(ylim=(0,0.6))
sns.boxplot(sites_filter_0[pi_sinapis],showfliers=False,palette=colour_sin).set(ylim=(0,0.0105))


sns.boxplot(data2[my_taj_D],showfliers=False,palette=pi_colours)
plt.axhline(y=0.0, linestyle='--')

sns.kdeplot(sites_filter_0_pi, shade=True);

sns.kdeplot(sites_filter_0[fst_comp_species],palette=pi_colours_species,showfliers=False)



#sns.distplot(sites_filter_0_pi, hist=False, color=pi_colours, kde_kws={"shade": True})


def dist_plots(DataFrame, pi_colours, sites_filter_0):
        for i,j in zip(DataFrame, pi_colours):
                sns.kdeplot(sites_filter_0[i],color=j, shade=True) 

dist_plots(my_pi, colour_new, sites_filter_0)
plt.xlim([0.0, 0.010])

sites_filter_0[['MEAN_FST_15','MEAN_FST_16','MEAN_FST_17']].plot.kde(shade=True,color=species_colour)
sites_filter_0[['MEAN_FST_15','MEAN_FST_16','MEAN_FST_17']].plot.kde()

dist_plots(['MEAN_FST_15','MEAN_FST_16','MEAN_FST_17'], pi_colours_species, sites_filter_0)




fig, axes = plt.subplots(nrows=2,ncols=2, sharex=True, sharey=True)
data.plot(ax=axes[0,0],kind='scatter', x='DXY_15', y='sinapis_reali_Avg_Pi', ylim=(0,0.025), xlim=(0,0.025), color='#66c2a5', alpha=0.2)
data_Z_FST_15.plot(ax=axes[0,0],kind='scatter', x='DXY_15', y='sinapis_reali_Avg_Pi', ylim=(0,0.025), xlim=(0,0.025), color='#9b9393',marker='^')
axes[0,0].plot([0, 1], [0, 1], color='#66c2a5')
data.plot(ax=axes[1,0],kind='scatter', x='DXY_16', y='reali_juvernica_Avg_Pi', ylim=(0,0.025), xlim=(0,0.025), color='#fa8e63', alpha=0.2)
data_Z_FST_16.plot(ax=axes[1,0],kind='scatter', x='DXY_16', y='reali_juvernica_Avg_Pi', ylim=(0,0.025), xlim=(0,0.025), color='#9b9393',marker='^')
axes[1,0].plot([0, 1], [0, 1], color='#fa8e63')
data.plot(ax=axes[1,1],kind='scatter', x='DXY_17', y='sinapis_juvernica_Avg_Pi', ylim=(0,0.025), xlim=(0,0.025), color='#8da0cb', alpha=0.2)
data_Z_FST_17.plot(ax=axes[1,1],kind='scatter', x='DXY_17', y='sinapis_juvernica_Avg_Pi', ylim=(0,0.025), xlim=(0,0.025), color='#9b9393',marker='^')
axes[1,1].plot([0, 1], [0, 1], color='#8da0cb')




fig, axes = plt.subplots(nrows=2,ncols=2, sharex=True, sharey=True)
data.plot(ax=axes[0,0],kind='scatter', x='DXY_0', y='kaz_sin_swe_sin_allele_Avg_Pi', ylim=(0,0.025), xlim=(0,0.025), color='#36ada4', alpha=0.2)
data_Z_FST_0.plot(ax=axes[0,0],kind='scatter', x='DXY_0', y='kaz_sin_swe_sin_allele_Avg_Pi', ylim=(0,0.025), xlim=(0,0.025), color='#9b9393',marker='^')
axes[0,0].plot([0, 1], [0, 1], color='#36ada4')
data.plot(ax=axes[0,1],kind='scatter', x='DXY_1', y='spanish_sinapis_swe_sin_allele_Avg_Pi', ylim=(0,0.025), xlim=(0,0.025), color='#39a7d0', alpha=0.2)
data_Z_FST_1.plot(ax=axes[0,1],kind='scatter', x='DXY_1', y='spanish_sinapis_swe_sin_allele_Avg_Pi', ylim=(0,0.025), xlim=(0,0.025), color='#9b9393', marker='^')
axes[0,1].plot([0, 1], [0, 1], color='#39a7d0')
data.plot(ax=axes[1,0],kind='scatter', x='DXY_2', y='kaz_sin_spanish_sinapis_Avg_Pi', ylim=(0,0.025), xlim=(0,0.025), color='#a48cf4', alpha=0.2)
data_Z_FST_2.plot(ax=axes[1,0],kind='scatter', x='DXY_2', y='kaz_sin_spanish_sinapis_Avg_Pi', ylim=(0,0.025), xlim=(0,0.025),  color='#9b9393', marker='^')
axes[1,0].plot([0, 1], [0, 1], color='#a48cf4')
data.plot(ax=axes[1,1],kind='scatter', x='DXY_3', y='irish_juvernica_kazak_juvernica_Avg_Pi', ylim=(0,0.025), xlim=(0,0.025), color='#f561dd', alpha=0.2)
data_Z_FST_3.plot(ax=axes[1,1],kind='scatter', x='DXY_3', y='irish_juvernica_kazak_juvernica_Avg_Pi', ylim=(0,0.025), xlim=(0,0.025),  color='#9b9393', marker='^')
axes[1,1].plot([0, 1], [0, 1], color='#f561dd')


violin_fst=pandas.concat([data2.MEAN_FST_0,data2.MEAN_FST_1,data2.MEAN_FST_2,data2.MEAN_FST_3,data2.MEAN_FST_15,data2.MEAN_FST_16,data2.MEAN_FST_17], ignore_index=True)
violin_auto=pandas.concat([data2.auto_allo,data2.auto_allo,data2.auto_allo,data2.auto_allo,data2.auto_allo,data2.auto_allo,data2.auto_allo], ignore_index=True)
d = {'': ['MEAN_FST_0']*len(data2.auto_allo)+['MEAN_FST_1']*len(data2.auto_allo)+['MEAN_FST_2']*len(data2.auto_allo)+['MEAN_FST_3']*len(data2.auto_allo)+['MEAN_FST_15']*len(data2.auto_allo)+['MEAN_FST_16']*len(data2.auto_allo)+['MEAN_FST_17']*len(data2.auto_allo)}
violin_ID= pandas.DataFrame(data=d)

violin_df = pandas.DataFrame()
violin_df['FST']=violin_fst
violin_df['auto_allo']=violin_auto
violin_df['ID']=violin_ID

ax = sns.violinplot(x="ID", y="FST", hue="auto_allo",data=violin_df, palette="Set2", split=True, inner="quartile")


##################################################################### 

####################get correlation coeficients#############################

def corr_pvalue(df):
        from scipy.stats import pearsonr
        import numpy as np
        import pandas as pd
        numeric_df = df.dropna()._get_numeric_data()
        cols = numeric_df.columns
        mat = numeric_df.values
        arr = np.zeros((len(cols),len(cols)), dtype=object)
        for xi, x in enumerate(mat.T):
                for yi, y in enumerate(mat.T[xi:]):
                    arr[xi, yi+xi] = map(lambda _: round(_,3), pearsonr(x,y))
                    arr[yi+xi, xi] = arr[xi, yi+xi]
        
        return pandas.DataFrame(arr, index=cols, columns=cols)


def corr_plot(df):
        axes= scatter_matrix(df,diagonal='none')
        corr = df.corr(method='pearson').as_matrix()
        mask = np.zeros_like(corr, dtype=np.bool)
        mask[np.tril_indices_from(mask)] = True
        #for i, j in zip(*plt.np.tril_indices_from(axes, k=1)):
        #    axes[i, j].annotate("%.3f" %corr[i,j], (0.1, 0.8), xycoords='axes fraction', ha='center', va='center')
        
        for i, j in zip(*plt.np.triu_indices_from(axes, k=1)):
            axes[i, j].set_visible(False)
        
        for i, j in zip(*plt.np.diag_indices_from(axes)):
            axes[i, j].set_visible(False)

        for i, j in zip(*plt.np.tril_indices_from(axes, k=1)):
            axes[i, j].set(xlim=(0,0.025))
            axes[i, j].set(ylim=(0,0.025))
        plt.show()


species_pi=["sinapis_my_PI","reali_my_PI","juvernica_my_PI" ]
POP_pi=['swe_sin_allele_my_PI','kaz_sin_my_PI','spanish_sinapis_my_PI','spanish_reali_my_PI','kazak_juvernica_my_PI','irish_juvernica_my_PI']


corr_pvalue(sites_filter_0[species_pi])
corr_pvalue(sites_filter_0[POP_pi])




corr_plot(sites_filter_0[species_pi])
corr_plot(sites_filter_0[POP_pi])


corr = sites_filter_0[POP_pi].corr(method='pearson').as_matrix()

sns.boxplot(sites_filter_0_fst[fst_sinapis_reali],showfliers=False).set(ylim=(0,0.6))
sns.boxplot(sites_filter_0_fst[fst_sinapis_juver],showfliers=False).set(ylim=(0,0.6))

sns.boxplot(sites_filter_0[fst_all_15],showfliers=False).set(ylim=(0,0.6))

#sites_filter_1=sites_filter[(sites_filter['CHROM'].isin(chr_1)) & (sites_filter['irish_juvernica_final_PI'] <= 1) & (sites_filter['kazak_juvernica_final_PI'] <= 1) & (sites_filter['kaz_sin_final_PI'] <= 1) & (sites_filter['spanish_reali_final_PI'] <= 1) & (sites_filter['spanish_sinapis_final_PI'] <= 1) & (sites_filter['swe_sin_final_PI'] <= 1)]

fig, axes = plt.subplots(nrows=3, ncols=2)
fig.subplots_adjust(hspace=0.5)


sns.regplot(x=sites_filter_0['kaz_sin_swe_sin_allele_dxy'],y=sites_filter_0['spanish_sinapis_swe_sin_allele_dxy'], ci=None, ax=axes[0,0]).set(xlabel="Dxy ,Kazakhstan sinapis Swedish Sinapis", ylabel="Dxy ,Spanish sinapis vs Swedish sinapis")
sns.regplot(x=sites_filter_0['kaz_sin_spanish_sinapis_dxy'],y=sites_filter_0['spanish_sinapis_swe_sin_allele_dxy'], ci=None, ax=axes[0,1]).set(xlabel="Dxy ,Kazakhstan sinapis Spanish Sinapis", ylabel="Dxy ,Spanish sinapis vs Swedish sinapis")
sns.regplot(x=sites_filter_0['kaz_sin_swe_sin_allele_dxy'], y=sites_filter_0['spanish_sinapis_swe_sin_allele_dxy'], ci=None, ax=axes[1,0]).set(xlabel="Dxy ,Kazakhstan sinapis Swedish Sinapis", ylabel="Dxy ,Spanish sinapis vs Swedish sinapis")


sns.regplot(x=sites_filter_0['irish_juvernica_kazak_juvernica_dxy'], y=sites_filter_0['spanish_sinapis_swe_sin_allele_dxy'], ax=axes[1,1], ci=None).set(xlabel="Dxy ,Irish juvernica Kazakhstan juvernica", ylabel="Dxy ,Spanish sinapis vs Swedish sinapis")
sns.regplot(x=sites_filter_0['irish_juvernica_kazak_juvernica_dxy'], y=sites_filter_0['kaz_sin_spanish_sinapis_dxy'],ax=axes[2,0], ci=None).set(xlabel="Dxy ,Irish juvernica Kazakhstan juvernica", ylabel="Dxy ,Kazakhstan sinapis vs Spanish sinapis")


sns.regplot(x=sites_filter_0['MEAN_FST_0'], y=sites_filter_0['MEAN_FST_3'], color='#212121', ci=None)
sns.regplot(x=sites_filter_0['MEAN_FST_1'], y=sites_filter_0['MEAN_FST_3'], color='#212121', ci=None)
sns.regplot(x=sites_filter_0['MEAN_FST_2'], y=sites_filter_0['MEAN_FST_3'], color='#212121', ci=None)
sns.regplot(x=sites_filter_0['MEAN_FST_8'], y=sites_filter_0['MEAN_FST_13'], color='#212121', ci=None)


sns.jointplot("irish_juvernica_rho", "irish_juvernica_final_PI", data=sites_filter_0, kind="reg")
sns.jointplot("kazak_juvernica_rho", "kazak_juvernica_final_PI", data=sites_filter_0, kind="reg")
sns.jointplot("kazak_sinapis_rho", "kaz_sin_final_PI", data=sites_filter_0, kind="reg")
sns.jointplot("spanish_reali_rho", "spanish_reali_final_PI", data=sites_filter_0, kind="reg")
sns.jointplot("spanish_sinapis_rho", "spanish_sinapis_final_PI", data=sites_filter_0, kind="reg")
sns.jointplot("swedish_sinapis_rho", "swe_sin_final_PI", data=sites_filter_0, kind="reg")

sns.boxplot(sites_filter_0[["irish_juvernica_rho","kazak_juvernica_rho","spanish_reali_rho","kazak_sinapis_rho","spanish_sinapis_rho","swedish_sinapis_rho"]],showfliers=False,palette=pi_colours).set(ylim=(0,0.00003))


sns.jointplot("MEAN_FST_0", "DXY_0", data=sites_filter_0, kind="reg")
sns.jointplot("MEAN_FST_1", "DXY_1", data=sites_filter_0, kind="reg")
sns.jointplot("MEAN_FST_2", "DXY_2", data=sites_filter_0, kind="reg")
sns.jointplot("MEAN_FST_3", "DXY_3", data=sites_filter_0, kind="reg")


sns.jointplot("MEAN_FST_15", "DXY_15", data=sites_filter_0, kind="reg")
sns.jointplot("MEAN_FST_16", "DXY_16", data=sites_filter_0, kind="reg")
sns.jointplot("MEAN_FST_17", "DXY_17", data=sites_filter_0, kind="reg")






sns.jointplot("kaz_sin_spanish_sinapis_fixed", "spanish_sinapis_swe_sin_allele_fixed", data=sites_filter_0, kind="reg")


ten_perc_alop=sites_filter_0[(sites_filter_0['MEAN_FST_0'] >= sites_filter_0.MEAN_FST_0.dropna().quantile(0.90)) & (sites_filter_0['DXY_0'] >= sites_filter_0.DXY_0.dropna().quantile(0.90)) & (sites_filter_0['MEAN_FST_1'] >= sites_filter_0.MEAN_FST_1.dropna().quantile(0.90)) & (sites_filter_0['DXY_1'] >= sites_filter_0.DXY_1.dropna().quantile(0.90)) & (sites_filter_0['MEAN_FST_2'] >= sites_filter_0.MEAN_FST_2.dropna().quantile(0.90)) & (sites_filter_0['DXY_2'] >= sites_filter_0.DXY_2.dropna().quantile(0.90)) ]

f = open('fasta_to10.fasta', 'a')
for i_index , i in ten_perc_alop.iterrows():
	sequence=str(fasta[i.CHROM][int(i.BIN_START):int(i.BIN_END)])
	f.write(">"+i.CHROM+"_"+str(int(i.BIN_START))+"_"+str(int(i.BIN_END))+"\n")
	f.write(sequence+"\n")

#############ploting########change chromosome to plot#####################
def chrom_plot(sites_filter_0,chromosome, pi, fst, dxy,fixed, colour):
	pi_selection=pi
	fst_selection=fst
	dxy_selection=dxy
	fst_set=[]
	pi_set=[]
	dxy_set=[]
	fixed_set=[]
	for i in chromosome:
		pi_list=sites_filter_0[(sites_filter_0['CHROM']==i)][pi_selection]
		fst_list=sites_filter_0[(sites_filter_0['CHROM']==i)][fst_selection]
		dxy_list=sites_filter_0[(sites_filter_0['CHROM']==i)][dxy_selection]
		fixed_1=sites_filter_0[(sites_filter_0['CHROM']==i)][fixed]
		if len(pi_list)>=1:
			pi_=pi_list.rolling(window=10)
			pi_set.append(pi_)
			#pi_.mean().plot(fontsize=10,legend=True)
		if len(fst_list)>=1:
			fst_=fst_list.rolling(window=10)
			fst_set.append(fst_)
		if len(dxy_list)>=1:
			dxy_=dxy_list.rolling(window=10)
			dxy_set.append(dxy_)
		if len(fixed_1)>=1:
			fixed_=fixed_1.rolling(window=1)
			fixed_set.append(fixed_)
	
	fixed_set2=pandas.DataFrame(columns=fixed)
	for f in fixed_set:
		list_f=[fixed_set2,f.mean()]
		fixed_set2=pandas.concat(list_f,ignore_index=True)

	
	fst_set2=pandas.DataFrame(columns=fst_selection)
	for j in fst_set:
		list_j=[fst_set2,j.mean()]
		fst_set2=pandas.concat(list_j,ignore_index=True)
	
	
	
	dxy_set2=pandas.DataFrame(columns=dxy_selection)
	for n in dxy_set:
		list_n=[dxy_set2,n.mean()]
		dxy_set2=pandas.concat(list_n,ignore_index=True)
	
	
	
	pi_set2=pandas.DataFrame(columns=pi_selection)
	for p in pi_set:
		list_p=[pi_set2,p.mean()]
		pi_set2=pandas.concat(list_p,ignore_index=True)
	
	
	fig, axes = plt.subplots(nrows=4, sharex=True)
	fig.subplots_adjust(hspace=0.1)
	fixed_set2.plot(ax=axes[0],legend=False).set_ylabel('fixed', fontsize=15)
	fst_set2.plot(ax=axes[1],legend=False).set_ylabel('Fst', fontsize=15)
	dxy_set2.plot(ax=axes[2], legend=False).set_ylabel('Dxy', fontsize=15)
	pi_set2.plot(ax=axes[3], legend=False ,color=colour).set_ylabel(r'$\pi$', fontsize=15)
	ticks = axes[3].get_xticks()/100
	axes[3].set_xticklabels(ticks)
	axes[3].set_xlabel("Mega bases")
	return fig





chrom_plot(data,chromosomes_dict['21'], pi_sinapis, fst_sinapis, dxy_sinapis,FIXED_0+FIXED_1+FIXED_2, colour_sin)
plt.tight_layout(pad=0.4, h_pad=1.0)

chrom_plot(data2,['scaffold_2'], pi_sinapis, fst_sinapis, dxy_sinapis,FIXED_0+FIXED_1+FIXED_2, colour_sin)



chrom_plot(data2,chromosomes_dict['21'], pi_sinapis_juver, fst_sinapis_juver, dxy_sinapis_juver,FIXED_0+FIXED_1+FIXED_2+FIXED_3, colour_sin_juve)
plt.tight_layout(pad=0.4, h_pad=1.0)






chrom_plot(sites_filter_0,chromosomes_dict['1'], pi, fst_all_15, dxy_all_15, colour_all)

chrom_plot(sites_filter_0,chromosomes_dict['1'], pi_new, fst_new, dxy_new, colour_new)
######################################################################################
def manhattan_plot(sites_filter_0, fst, dxy,chromosomes_dict):
	chromosome=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21']
	fst_selection=fst+['CHROM']
	dxy_selection=dxy+['CHROM']
	fst_set=[]
	dxy_set=[]
	fixed_set=[]
	chrom={}
	for j in chromosome:
		chrom[j]=0
		for i in chromosomes_dict[str(j)]:
			chrom[j]+=len(sites_filter_0[(sites_filter_0['CHROM']==i)])
			fst_list=sites_filter_0[(sites_filter_0['CHROM']==i)][fst_selection]
			dxy_list=sites_filter_0[(sites_filter_0['CHROM']==i)][dxy_selection]
			if len(fst_list)>=5:
				fst_=fst_list.rolling(window=20)
				fst_set.append(fst_)
			if len(dxy_list)>=5:
				dxy_=dxy_list.rolling(window=20)
				dxy_set.append(dxy_)
	
	
	fst_set2=pandas.DataFrame(columns=fst_selection)
	for j in fst_set:
		list_j=[fst_set2,j.mean()]
		fst_set2=pandas.concat(list_j,ignore_index=True)
	
		
	dxy_set2=pandas.DataFrame(columns=dxy_selection)
	for n in dxy_set:
		list_n=[dxy_set2,n.mean()]
		dxy_set2=pandas.concat(list_n,ignore_index=True)

	fig, axes = plt.subplots(nrows=2, sharex=True)
	fig.subplots_adjust(hspace=0.1)
	fst_set2.plot(ax=axes[0],legend=False).set_ylabel('Fst', fontsize=15)
	dxy_set2.plot(ax=axes[1],legend=False).set_ylabel('Dxy', fontsize=15)
	li=[]
	for c in chromosome:
		if c == '1':
			start=0
			end=int(chrom[str(c)])
		else:
			start=end
			end=start+int(chrom[str(c)])
		li.append([start, end])
	for l in li[0::2]:
		axes[0].axvspan(l[0], l[1], facecolor='#efefef', alpha=0.5)
		axes[1].axvspan(l[0], l[1], facecolor='#efefef', alpha=0.5)
	ticks = axes[1].get_xticks()/100
	axes[1].set_xticklabels(ticks)
	axes[1].set_xlabel("Mega bases")
	return fig



manhattan_plot(data2, fst_all_15, dxy_all_15,chromosomes_dict)

def manhattan_plot_fst_dxy_pi(sites_filter_0, fst, dxy,pi,taj,fixed,chromosomes_dict,colour):
	species_colour=sns.color_palette("Set2", 8)[:3]
        chromosome=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21']
	pi_selection=pi
	fst_selection=fst+['CHROM']
	dxy_selection=dxy+['CHROM']
        taj_selection=taj
	fst_set=[]
	dxy_set=[]
	fixed_set=[]
	pi_set=[]
        taj_set=[]
	chrom={}
	for j in chromosome:
		chrom[j]=0
		for i in chromosomes_dict[str(j)]:
			pi_list=sites_filter_0[(sites_filter_0['CHROM']==i)][pi_selection]
			fixed_1=sites_filter_0[(sites_filter_0['CHROM']==i)][fixed]
			chrom[j]+=len(sites_filter_0[(sites_filter_0['CHROM']==i)])
			fst_list=sites_filter_0[(sites_filter_0['CHROM']==i)][fst_selection]
			dxy_list=sites_filter_0[(sites_filter_0['CHROM']==i)][dxy_selection]
                        taj_list=sites_filter_0[(sites_filter_0['CHROM']==i)][taj_selection]
			if len(fst_list)>=5:
				fst_=fst_list.rolling(window=20)
				fst_set.append(fst_)
			if len(dxy_list)>=5:
				dxy_=dxy_list.rolling(window=20)
				dxy_set.append(dxy_)
			if len(pi_list)>=5:
				pi_=pi_list.rolling(window=20)
				pi_set.append(pi_)
				#pi_.mean().plot(fontsize=10,legend=True)
			if len(fixed_1)>=1:
				fixed_=fixed_1.rolling(window=1)
				fixed_set.append(fixed_)
                        if len(taj_list)>=5:
                                taj_=taj_list.rolling(window=20)
                                taj_set.append(taj_)
                                #pi_.mean().plot(fontsize=10,legend=True)
	
	fst_set2=pandas.DataFrame(columns=fst_selection)
	for j in fst_set:
		list_j=[fst_set2,j.mean()]
		fst_set2=pandas.concat(list_j,ignore_index=True)
	
		
	dxy_set2=pandas.DataFrame(columns=dxy_selection)
	for n in dxy_set:
		list_n=[dxy_set2,n.mean()]
		dxy_set2=pandas.concat(list_n,ignore_index=True)

	fixed_set2=pandas.DataFrame(columns=fixed)
	for f in fixed_set:
		list_f=[fixed_set2,f.mean()]
		fixed_set2=pandas.concat(list_f,ignore_index=True)

	pi_set2=pandas.DataFrame(columns=pi_selection)
	for p in pi_set:
		list_p=[pi_set2,p.mean()]
		pi_set2=pandas.concat(list_p,ignore_index=True)

        taj_set2=pandas.DataFrame(columns=taj_selection)
        for p in taj_set:
                list_p=[taj_set2,p.mean()]
                taj_set2=pandas.concat(list_p,ignore_index=True)

	fig, axes = plt.subplots(nrows=5, sharex=True)
	fig.subplots_adjust(hspace=0.1)
	fixed_set2.plot(ax=axes[0],legend=False,color=species_colour).set_ylabel('fixed', fontsize=15)
	fst_set2.plot(ax=axes[1],legend=False,color=species_colour).set_ylabel('Fst', fontsize=15)
	dxy_set2.plot(ax=axes[2], legend=False,color=species_colour).set_ylabel('Dxy', fontsize=15)
	pi_set2.plot(ax=axes[3], legend=False ,color=colour).set_ylabel(r'$\pi$', fontsize=15)
        taj_set2.plot(ax=axes[4], legend=False ,color=colour).set_ylabel(r'$\pi$', fontsize=15)
	li=[]
	for c in chromosome:
		if c == '1':
			start=0
			end=int(chrom[str(c)])
		else:
			start=end
			end=start+int(chrom[str(c)])
		li.append([start, end])
	for l in li[0::2]:
		axes[0].axvspan(l[0], l[1], facecolor='#dbd9d9', alpha=0.5)
		axes[1].axvspan(l[0], l[1], facecolor='#dbd9d9', alpha=0.5)
		axes[2].axvspan(l[0], l[1], facecolor='#dbd9d9', alpha=0.5)
		axes[3].axvspan(l[0], l[1], facecolor='#dbd9d9', alpha=0.5)
                axes[4].axvspan(l[0], l[1], facecolor='#dbd9d9', alpha=0.5)
	ticks = axes[4].get_xticks()/100
	axes[4].set_xticklabels(ticks)
	axes[4].set_xlabel("Mega bases")
	return fig


taj_sinapis=['kaz_sin_TajimaD', 'spanish_sinapis_TajimaD', 'swe_sin_TajimaD']
pi_sinapis=['kaz_sin_my_PI', 'spanish_sinapis_my_PI', 'swe_sin_allele_my_PI']


manhattan_plot_fst_dxy_pi(data2, fst_sinapis ,dxy_sinapis ,pi_sinapis,taj_sinapis,FIXED_0+FIXED_1+FIXED_2, chromosomes_dict,colour_sin)

manhattan_plot_fst_dxy_pi(data2, fst_sinapis_juver ,dxy_sinapis_juver ,pi_sinapis_juver,FIXED_0+FIXED_1+FIXED_2+FIXED_3, chromosomes_dict,colour_sin_juve)

manhattan_plot_fst_dxy_pi(data2, ['MEAN_FST_4','MEAN_FST_8'] ,['DXY_4','DXY_8'] ,['kaz_sin_final_PI','irish_juvernica_final_PI','kazak_juvernica_final_PI'],FIXED_4+FIXED_8, chromosomes_dict,['#FF8C00', '#C8C800', '#006400'])


manhattan_plot_fst_dxy_pi(data2,  ['Z_FST_15', 'Z_FST_16', 'Z_FST_17'] ,dxy_comp_species ,pi_species,FIXED_15+FIXED_16+FIXED_17, chromosomes_dict,pi_colours_species)


manhattan_plot_fst_dxy_pi(data2,  fst_all_15 ,dxy_all_15 ,my_pi, my_taj_D ,FIXED_15+FIXED_16+FIXED_17, chromosomes_dict,pi_colours_species)

my_taj_D
my_pi
fst_all_15
dxy_all_15


fst_comp_species=['MEAN_FST_15','MEAN_FST_16','MEAN_FST_17']
dxy_comp_species=['DXY_15','DXY_16','DXY_17']

col_comp_species=sns.color_palette(['#56AACC','#B6AAE8','#E88CD9'])


def manhattan_plot_fst_sym(sites_filter_0, fst,a, colour):
        #sites_filter_0 is the dataframe
        #fst is the name of the coloum in dataframe to plot
        #'a' is the value for the line 
        chromosome=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21']
        fst_selection=fst+['CHROM']
        fst_set=[]
        chrom={}
        for j in chromosome:
                chrom[j]=0
                for i in chromosomes_dict[str(j)]:
                        chrom[j]+=len(sites_filter_0[(sites_filter_0['CHROM']==i)])
                        fst_list=sites_filter_0[(sites_filter_0['CHROM']==i)][fst_selection]
                        if len(fst_list)>=5:
                                fst_=fst_list.rolling(window=20)
                                fst_set.append(fst_)
 
        fst_set2=pandas.DataFrame(columns=fst_selection)
        for j in fst_set:
                list_j=[fst_set2,j.mean()]
                fst_set2=pandas.concat(list_j,ignore_index=True)
         
 
        fig, axes = plt.subplots(nrows=2, sharex=True)
        fst_set2.plot(ax=axes[0],color=[colour],legend=False).set_ylabel('Fst', fontsize=15)
        axes[0].axhline(y=a[0], xmin=0.0, xmax=0.96097,color=colour, drawstyle= 'steps')
        axes[0].axhline(y=a[1], xmin=0.96097, xmax=1.0,color=colour, drawstyle= 'steps')
        li=[]
        for c in chromosome:
                if c == '1':
                        start=0
                        end=int(chrom[str(c)])
                else:
                        start=end
                        end=start+int(chrom[str(c)])
                li.append([start, end])
        for l in li[0::2]:
                axes[0].axvspan(l[0], l[1], facecolor='#dbd9d9', alpha=0.5)
        ticks = axes[0].get_xticks()/100
        axes[0].set_xticklabels(ticks)
        axes[0].set_xlabel("Mega bases")
        return fig
 



##############################################roh##################PLOT##############

chromosome=chromosomes_dict['6']
#chromosome=['scaffold_2']
rho_selection=['RHO']
roh_data1=pandas.read_table("irish_juvernica.rho")
roh_data1.BIN_START=roh_data1.BIN_START+1
roh_data1.RHO=roh_data1.RHO/10000
l=[sites_filter,roh_data1]
rho_df=personal_popgen.join_raw_data_base(l)

rho_set=[]
for i in chromosome:
	rho_list=rho_df[(rho_df['CHROM']==i)][rho_selection]
	if len(rho_list)>=1:
		rho_=rho_list.rolling(window=20)
		rho_set.append(rho_)
		#pi_.mean().plot(fontsize=10,legend=True)

rho_set2=pandas.DataFrame(columns=rho_selection)
for j in rho_set:
	list_j=[rho_set2,j.mean()]
	rho_set2=pandas.concat(list_j,ignore_index=True)

rho_set2.plot()


###################################fixed_private_shared_same##########################################

labels = 'Fixed', 'Private A', 'Private B', 'shared'
list_sites=["sinapis_reali_","reali_juvernica_","sinapis_juvernica_","kaz_sin_swe_sin_allele_","spanish_sinapis_swe_sin_allele_","kaz_sin_spanish_sinapis_","irish_juvernica_kazak_juvernica_","irish_juvernica_kaz_sin_","irish_juvernica_spanish_reali_","irish_juvernica_spanish_sinapis_","irish_juvernica_swe_sin_allele_","kazak_juvernica_kaz_sin_","kazak_juvernica_spanish_reali_","kazak_juvernica_spanish_sinapis_","kazak_juvernica_swe_sin_allele_","kaz_sin_spanish_reali_","spanish_reali_spanish_sinapis_","spanish_reali_swe_sin_allele_"]



for i in list_sites:
	j=[i+'fixed', i+'private_a', i+'private_b', i+'shared']
	n=list(data[j].sum())
	print i
	print list(data[j].sum())
	print sum(list(data[j].sum()))
	sizes=[(j/np.nansum(n)) *100   for j in n]
	explode = (0, 0, 0, 0)
	#fig1, ax1 = plt.subplots()
	#ax1.pie(sizes, explode=explode, labels=labels, autopct='%1.1f%%', shadow=True, startangle=90)
	#ax1.axis('equal')
	#ax1.set_title(i[:-1], bbox={'facecolor':'0.8', 'pad':5})


plt.show()


labels = 'Fixed', 'Private A', 'Private B', 'shared'
list_sites=["sinapis_reali_","reali_juvernica_","sinapis_juvernica_"]



for i in list_sites:
    j=[i+'fixed', i+'private_a', i+'private_b', i+'shared']
    n=list(data[j].sum())
    print i
    print list(data[j].sum())
    print sum(list(data[j].sum()))
    sizes=[(j/np.nansum(n)) *100   for j in n]
    explode = (0, 0, 0, 0)
    fig1, ax1 = plt.subplots()
    ax1.pie(sizes, explode=explode, labels=labels, autopct='%1.1f%%', shadow=True, startangle=90)
    ax1.axis('equal')
    ax1.set_title(i[:-1], bbox={'facecolor':'0.8', 'pad':5})


plt.show()
###############################rough#####################################
ploter=sites_filter_1[sinapis_pi]
ploter=ploter.reset_index(drop=True)
r=ploter.rolling(window=10)
r.mean().plot(fontsize=10,legend=True, color=['#FF8C00','#C04000','#FF0000'])
plt.show()

ploter=sites_filter_1[sinapis_fst]
ploter=ploter[(ploter['MEAN_FST_0'] >= 0)&(ploter['MEAN_FST_1'] >= 0)&(ploter['MEAN_FST_2'] >= 0)]

ploter=ploter.rename(columns={'MEAN_FST_0': 'Kazak_sinapis_vs_Swedish_sinapis', 'MEAN_FST_1': 'Spanish_sinapis_vs_Swedish_sinapis', 'MEAN_FST_2': 'Kazak_sinapis_vs_Spanish_sinapis'})
ploter=ploter.reset_index(drop=True)
r=ploter.rolling(window=10)
r.mean().plot(fontsize=10,legend=True, color=['#8000ff','#0040ff','#a65959'])
plt.show()

"L. sinapis (Kaz)" "#FF8C00"
"L. sinapis (Spa)" "#FF0000"
"L. sinapis (Swe)" "#C04000"
"L. juvernica (Kaz)" "#006400"
"L. reali (Spa)" "#0000FF"
"L. juvernica (Ire)" "#C8C800"

['#C8C800','#006400','#FF8C00','#0000FF','#FF0000','#C04000']



all_1=[]
for i in list_sites:
	all_1.append(i+'fixed')
	all_1.append(i+'private_a')
	all_1.append(i+'private_b')
	all_1.append(i+'shared')

data_fixed_private_shared=data[['CHROM','BIN_START','BIN_END']+all_1]

data_fixed_private_shared.to_csv('data_fixed_private_shared')





a[(a['Z_FST_3'] >= a.Z_FST_3.dropna().quantile(0.99)) & (a['Z_FST_2'] >= a.Z_FST_2.dropna().quantile(0.99))] 

z[(z['Z_FST_3'] >= z.Z_FST_3.dropna().quantile(0.99)) & (z['Z_FST_2'] >= z.Z_FST_2.dropna().quantile(0.99))] 



a[(a['Z_FST_3'] >= a.Z_FST_3.dropna().quantile(0.99))].[]
z[(z['Z_FST_3'] >= z.Z_FST_3.dropna().quantile(0.99))].[]


a_col=species_colour=sns.color_palette("Set2", 8)[:3][0]
b_col=species_colour=sns.color_palette("Set2", 8)[:3][1]
c_col=species_colour=sns.color_palette("Set2", 8)[:3][2]

d_col=sns.color_palette("husl", 8)[-4:][0]
e_col=sns.color_palette("husl", 8)[-4:][1]
f_col=sns.color_palette("husl", 8)[-4:][2]
g_col=sns.color_palette("husl", 8)[-4:][3]


manhattan_plot_fst_sym(data2, ['Z_FST_15'],[a.Z_FST_15.dropna().quantile(0.99),z.Z_FST_15.dropna().quantile(0.99)],a_col)
manhattan_plot_fst_sym(data2, ['Z_FST_16'],[a.Z_FST_16.dropna().quantile(0.99),z.Z_FST_16.dropna().quantile(0.99)],b_col)
manhattan_plot_fst_sym(data2, ['Z_FST_17'],[a.Z_FST_17.dropna().quantile(0.99),z.Z_FST_17.dropna().quantile(0.99)],c_col)


manhattan_plot_fst_sym(data2, ['Z_FST_0'],[a.Z_FST_0.dropna().quantile(0.99),z.Z_FST_0.dropna().quantile(0.99)],d_col)
manhattan_plot_fst_sym(data2, ['Z_FST_1'],[a.Z_FST_1.dropna().quantile(0.99),z.Z_FST_1.dropna().quantile(0.99)],e_col)
manhattan_plot_fst_sym(data2, ['Z_FST_2'],[a.Z_FST_2.dropna().quantile(0.99),z.Z_FST_2.dropna().quantile(0.99)],f_col)
manhattan_plot_fst_sym(data2, ['Z_FST_3'],[a.Z_FST_3.dropna().quantile(0.99),z.Z_FST_3.dropna().quantile(0.99)],g_col)

sns.palplot([a_col,b_col,c_col,d_col,e_col,f_col,g_col])



def manhattan_plot_fst_dxy_pi(sites_filter_0, fst, dxy,pi,taj,fixed,chromosomes_dict,colour,species_colour):
        #species_colour=sns.color_palette("husl", 8)[-4:][0:3]
        chromosome=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21']
        pi_selection=pi
        fst_selection=fst+['CHROM']
        dxy_selection=dxy+['CHROM']
        taj_selection=taj
        fst_set=[]
        dxy_set=[]
        fixed_set=[]
        pi_set=[]
        taj_set=[]
        chrom={}
        for j in chromosome:
                chrom[j]=0
                for i in chromosomes_dict[str(j)]:
                        pi_list=sites_filter_0[(sites_filter_0['CHROM']==i)][pi_selection]
                        fixed_1=sites_filter_0[(sites_filter_0['CHROM']==i)][fixed]
                        chrom[j]+=len(sites_filter_0[(sites_filter_0['CHROM']==i)])
                        fst_list=sites_filter_0[(sites_filter_0['CHROM']==i)][fst_selection]
                        dxy_list=sites_filter_0[(sites_filter_0['CHROM']==i)][dxy_selection]
                        taj_list=sites_filter_0[(sites_filter_0['CHROM']==i)][taj_selection]
                        if len(fst_list)>=5:
                                fst_=fst_list.rolling(window=20)
                                fst_set.append(fst_)
                        if len(dxy_list)>=5:
                                dxy_=dxy_list.rolling(window=20)
                                dxy_set.append(dxy_)
                        if len(pi_list)>=5:
                                pi_=pi_list.rolling(window=20)
                                pi_set.append(pi_)
                                #pi_.mean().plot(fontsize=10,legend=True)
                        if len(fixed_1)>=1:
                                fixed_=fixed_1.rolling(window=1)
                                fixed_set.append(fixed_)
                        if len(taj_list)>=5:
                                taj_=taj_list.rolling(window=20)
                                taj_set.append(taj_)
                                #pi_.mean().plot(fontsize=10,legend=True)
        
        fst_set2=pandas.DataFrame(columns=fst_selection)
        for j in fst_set:
                list_j=[fst_set2,j.mean()]
                fst_set2=pandas.concat(list_j,ignore_index=True)
        
                
        dxy_set2=pandas.DataFrame(columns=dxy_selection)
        for n in dxy_set:
                list_n=[dxy_set2,n.mean()]
                dxy_set2=pandas.concat(list_n,ignore_index=True)

        fixed_set2=pandas.DataFrame(columns=fixed)
        for f in fixed_set:
                list_f=[fixed_set2,f.mean()]
                fixed_set2=pandas.concat(list_f,ignore_index=True)

        pi_set2=pandas.DataFrame(columns=pi_selection)
        for p in pi_set:
                list_p=[pi_set2,p.mean()]
                pi_set2=pandas.concat(list_p,ignore_index=True)

        taj_set2=pandas.DataFrame(columns=taj_selection)
        for p in taj_set:
                list_p=[taj_set2,p.mean()]
                taj_set2=pandas.concat(list_p,ignore_index=True)

        fig, axes = plt.subplots(nrows=5, sharex=True)
        fig.subplots_adjust(hspace=0.1)
        fixed_set2.plot(ax=axes[0],legend=False,color=species_colour).set_ylabel('fixed', fontsize=15)
        fst_set2.plot(ax=axes[1],legend=False,color=species_colour).set_ylabel('Fst', fontsize=15)
        dxy_set2.plot(ax=axes[2], legend=False,color=species_colour).set_ylabel('Dxy', fontsize=15)
        pi_set2.plot(ax=axes[3], legend=False ,color=colour).set_ylabel(r'$\pi$', fontsize=15)
        taj_set2.plot(ax=axes[4], legend=False ,color=colour).set_ylabel(r'$\pi$', fontsize=15)
        li=[]
        for c in chromosome:
                if c == '1':
                        start=0
                        end=int(chrom[str(c)])
                else:
                        start=end
                        end=start+int(chrom[str(c)])
                li.append([start, end])
        for l in li[0::2]:
                axes[0].axvspan(l[0], l[1], facecolor='#dbd9d9', alpha=0.5)
                axes[1].axvspan(l[0], l[1], facecolor='#dbd9d9', alpha=0.5)
                axes[2].axvspan(l[0], l[1], facecolor='#dbd9d9', alpha=0.5)
                axes[3].axvspan(l[0], l[1], facecolor='#dbd9d9', alpha=0.5)
                axes[4].axvspan(l[0], l[1], facecolor='#dbd9d9', alpha=0.5)
        ticks = axes[4].get_xticks()/100
        axes[4].set_xticklabels(ticks)
        axes[4].set_xlabel("Mega bases")
        return fig




taj_sinapis=['kaz_sin_TajimaD', 'spanish_sinapis_TajimaD', 'swe_sin_TajimaD']
pi_sinapis=['kaz_sin_my_PI', 'spanish_sinapis_my_PI', 'swe_sin_allele_my_PI']


manhattan_plot_fst_dxy_pi(data2, fst_sinapis ,dxy_sinapis ,pi_sinapis,taj_sinapis,FIXED_0+FIXED_1+FIXED_2, chromosomes_dict,colour_sin, sns.color_palette("husl", 8)[-4:][0:3])

manhattan_plot_fst_dxy_pi(data2, ['MEAN_FST_3'] ,['DXY_0'] ,['irish_juvernica_my_PI','kazak_juvernica_my_PI'],['irish_juvernica_TajimaD','kazak_juvernica_TajimaD'],FIXED_3, chromosomes_dict,['#C8C800', '#006400'],'#F461DD')



