import sys
import os
import string
import pandas
import personal_popgen
import itertools
import numpy as np




pwd='/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/'
os.chdir(pwd)



#load FST windows
file1= "/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/kazak_sinapis_Swedish_sinapis_10000.windowed.weir.fst"    #MEAN_FST_0
file2= "/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/spanish_sinapis_Swedish_sinapis_10000.windowed.weir.fst"  #MEAN_FST_1
file3= "/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/kazak_sinapis_spanish_sinapis_10000.windowed.weir.fst"    #MEAN_FST_2
file4= "/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/irish_juvernica_kazak_juvernica_10000.windowed.weir.fst"  #MEAN_FST_3
file5= "/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/irish_juvernica_kazak_sinapis_10000.windowed.weir.fst"    #MEAN_FST_4
file6= "/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/irish_juvernica_spanish_reali_10000.windowed.weir.fst"    #MEAN_FST_5
file7= "/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/irish_juvernica_spanish_sinapis_10000.windowed.weir.fst"  #MEAN_FST_6
file8= "/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/irish_juvernica_Swedish_sinapis_10000.windowed.weir.fst"  #MEAN_FST_7
file9= "/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/kazak_juvernica_kazak_sinapis_10000.windowed.weir.fst"    #MEAN_FST_8
file10="/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/kazak_juvernica_spanish_reali_10000.windowed.weir.fst"    #MEAN_FST_9
file11="/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/kazak_juvernica_spanish_sinapis_10000.windowed.weir.fst"  #MEAN_FST_10
file12="/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/kazak_juvernica_Swedish_sinapis_10000.windowed.weir.fst"  #MEAN_FST_11
file13="/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/kazak_sinapis_spanish_reali_10000.windowed.weir.fst"      #MEAN_FST_12
file14="/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/spanish_reali_spanish_sinapis_10000.windowed.weir.fst"    #MEAN_FST_13
file15="/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/spanish_reali_Swedish_sinapis_10000.windowed.weir.fst"    #MEAN_FST_14

file16="/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/sinapis_reali_10000.windowed.weir.fst"			          #MEAN_FST_15
file17="/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/juvernica_reali_10000.windowed.weir.fst"			      #MEAN_FST_16
file18="/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/sinapis_juvernica_10000.windowed.weir.fst"			      #MEAN_FST_17




# load CSV file with sites coverd in each window#load FST windows
sites_kazak_sinapis   ="/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/kazak_sinapis_10000.csv"
sites_Swedish_sinapis ="/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/Swedish_sinapis_10000.csv"
sites_spanish_sinapis ="/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/spanish_sinapis_10000.csv"
sites_irish_juvernica ="/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/irish_juv_10000.csv"
sites_kazak_juvernica ="/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/kazak_juv_10000.csv"
sites_spanish_reali   ="/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/spanish_reali_10000.csv"

#
#
#

# load site pi
irish_juvernica_site_PI ="/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/irish_juvernica_site_PI.sites.pi.sites.pi"
kazak_juvernica_site_PI ="/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/kazak_juvernica_site_PI.sites.pi.sites.pi"
kaz_sin_site_PI         ="/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/kaz_sin_site_PI.sites.pi.sites.pi"
spanish_reali_site_PI   ="/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/spanish_reali_site_PI.sites.pi.sites.pi"
spanish_sinapis_site_PI ="/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/spanish_sinapis_site_PI.sites.pi.sites.pi"
swe_sin_site_PI         ="/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/swe_sin_site_PI.sites.pi.sites.pi"

sinapis_site_PI         ="/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/sinapis.PI.sites.sites.pi.sites.pi"
reali_site_PI           ="/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/reali.PI.sites.sites.pi.sites.pi"
juvernica_site_PI       ="/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/juvernica.PI.sites.sites.pi.sites.pi"



# load allele freq files
irish_juvernica_allele_frquencey ="/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/irish_juvernica_freq.frq"
kazak_juvernica_allele_frquencey ="/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/kazak_juvernica_freq.frq"
kaz_sin_allele_frquencey         ="/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/kaz_sin_freq.frq"
spanish_reali_allele_frquencey   ="/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/spanish_reali_freq.frq"
spanish_sinapis_allele_frquencey ="/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/spanish_sinapis_freq.frq"
swe_sin_allele_frquencey         ="/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/swe_sin_allele_freq.frq"

sinapis_frequencey="/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/sinapis.frq"
reali_frequencey="/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/reali.frq"
juvernica_frequencey="/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/juvernica.frq"




test1 =pandas.read_table(file1)
test2 =pandas.read_table(file2)
test3 =pandas.read_table(file3)
test4 =pandas.read_table(file4)
test5 =pandas.read_table(file5)
test6 =pandas.read_table(file6)
test7 =pandas.read_table(file7)
test8 =pandas.read_table(file8)
test9 =pandas.read_table(file9)
test10=pandas.read_table(file10)
test11=pandas.read_table(file11)
test12=pandas.read_table(file12)
test13=pandas.read_table(file13)
test14=pandas.read_table(file14)
test15=pandas.read_table(file15)

test16=pandas.read_table(file16)
test17=pandas.read_table(file17)
test18=pandas.read_table(file18)



irish_juvernica_site_PI_table =pandas.read_table(irish_juvernica_site_PI)
kazak_juvernica_site_PI_table =pandas.read_table(kazak_juvernica_site_PI)
kaz_sin_site_PI_table         =pandas.read_table(kaz_sin_site_PI)
spanish_reali_site_PI_table   =pandas.read_table(spanish_reali_site_PI)
spanish_sinapis_site_PI_table =pandas.read_table(spanish_sinapis_site_PI)
swe_sin_site_PI_table         =pandas.read_table(swe_sin_site_PI)
sinapis_site_PI_table         =pandas.read_table(sinapis_site_PI)
reali_site_PI_table           =pandas.read_table(reali_site_PI)
juvernica_site_PI_table       =pandas.read_table(juvernica_site_PI)


sites_kazak_sinapis_table  =pandas.read_csv(sites_kazak_sinapis)
sites_Swedish_sinapis_table=pandas.read_csv(sites_Swedish_sinapis)
sites_spanish_sinapis_table=pandas.read_csv(sites_spanish_sinapis)
sites_irish_juvernica_table=pandas.read_csv(sites_irish_juvernica)
sites_kazak_juvernica_table=pandas.read_csv(sites_kazak_juvernica)
sites_spanish_reali_table  =pandas.read_csv(sites_spanish_reali)

sites=pandas.read_csv('/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/FINAL_db_N_sites_1x_10_ind')

lists_1=[test1,test2,test3,test4,test5,test6,test7,test8,test9,test10,test11,test12,test13,test14,test15,test16,test17,test18,sites]

nextone=personal_popgen.join_data_base(lists_1)
primary_keys=nextone[['CHROM','BIN_START','BIN_END']]



###########in process###########
def callcualte_windowed_pi_from_site_pi(population_pi, dataframe, sites_list, output):
	win_all_sites=dataframe[['CHROM','BIN_START','BIN_END',sites_list]]
	win_sites=population_pi
	win_sites['BIN_START']=(np.floor(win_sites['POS']/10000)*10000)+1
	win_sites['BIN_END']=win_sites['BIN_START']+(10000-1)
	df_merge = pandas.merge(win_sites, win_all_sites, on=['CHROM','BIN_START','BIN_END'],how='inner')
	df_merge = df_merge.query('BIN_START <= POS and BIN_END >= POS')
	out=pandas.DataFrame({output+'_final_pi' :df_merge.groupby(['CHROM','BIN_START','BIN_END',sites_list], as_index=True)['PI'].sum()}).reset_index()
	out[output+'_final_pi']=out[output+'_final_pi']/out[sites_list]
	return out


#callcualte_windowed_pi_from_site_pi(irish_juvernica_site_PI_table, nextone, 'irish_juvernica_N_sites', 'irish_juvernica_PI')

#site_pi=[irish_juvernica_site_PI_table,kazak_juvernica_site_PI_table,kaz_sin_site_PI_table,spanish_reali_site_PI_table,spanish_sinapis_site_PI_table,swe_sin_site_PI_table,sinapis_site_PI_table,reali_site_PI_table,juvernica_site_PI_table]
#final_window_pi=['irish_juvernica','kazak_juvernica','kaz_sin','spanish_reali','spanish_sinapis','swe_sin',"sinapis_site","reali_site","juvernica_site"]
#final_site_pi=['irish_juvernica_N_sites','kazak_juvernica_N_sites','kaz_sin_N_sites','spanish_reali_N_sites','spanish_sinapis_N_sites','swe_sin_allele_N_sites',"sinapis_N_sites","reali_N_sites","juvernica_N_sites"]
#
#new_stuff=[]
#for q, a, s in zip(site_pi, final_window_pi, final_site_pi):
#	#print [str(a)]
#	print " "+ a + " "+ s + " "
#	new=callcualte_windowed_pi_from_site_pi(q, nextone, s, a)
#	new_stuff.append(new)
#
#nextone=personal_popgen.join_data_base(new_stuff)
#
nextone=0

#sites_filter[[i+'_final_pi' for i in final_window_pi]].mean()

#sites_filter=nextone[(nextone['irish_juvernica_N_sites'] >= 3000) & (nextone['kazak_juvernica_N_sites'] >= 3000) & (nextone['kaz_sin_N_sites'] >= 3000) & (nextone['spanish_reali_N_sites'] >= 3000) & (nextone['spanish_sinapis_N_sites'] >= 3000)& (nextone['swe_sin_allele_N_sites'] >= 3000)& (nextone['sinapis_N_sites'] >= 3000)& (nextone['reali_N_sites'] >= 3000)& (nextone['juvernica_N_sites'] >= 3000)]

##population script 
freq_list=["irish_juvernica","kazak_juvernica","kaz_sin","spanish_reali","spanish_sinapis","swe_sin_allele"]
dxy_comb={str(i[0]+'_'+i[1]+'_dxy'): [pwd+(i[0]+'_freq.frq'), (pwd+i[1]+'_freq.frq')] for i in list(itertools.combinations(freq_list, 2))}
columns={str(i[0]+'_'+i[1]+'_dxy'): [(i[0]), (i[1])] for i in list(itertools.combinations(freq_list, 2))}
for i in dxy_comb.keys():
	filea=dxy_comb[i][0]
	fileb=dxy_comb[i][1]
	output_col=i
	fstcol=i[:-3]+'fst'
	fixedcol=i[:-3]+'fixed'
	privatea=i[:-3]+'private_a'
	privateb=i[:-3]+'private_b'
	sharedcol=i[:-3]+'shared'
	samecol=i[:-3]+'same_sites'
	primary_keys[str(output_col)]=0.0
	primary_keys[str(fstcol)]=0.0
	primary_keys[str(fixedcol)]=0.0
	primary_keys[str(privatea)]=0.0
	primary_keys[str(privateb)]=0.0
	primary_keys[str(sharedcol)]=0.0
	primary_keys[str(samecol)]=0.0
	for j in personal_popgen.dxy_window_function(filea,fileb,primary_keys):
		#output_list [window_index, np.nanmean(avg_dxy), np.nanmean(avg_fst),fixed_diff, private_pop1, private_pop2, shared, fixed_same]
		primary_keys[str(output_col)]=primary_keys[str(output_col)].set_value(j[0], value=j[1])
		primary_keys[str(fstcol)]=primary_keys[str(fstcol)].set_value(j[0], value=j[2])
		primary_keys[str(fixedcol)]=primary_keys[str(fixedcol)].set_value(j[0], value=j[3])
		primary_keys[str(privatea)]=primary_keys[str(privatea)].set_value(j[0], value=j[4])
		primary_keys[str(privateb)]=primary_keys[str(privateb)].set_value(j[0], value=j[5])
		primary_keys[str(sharedcol)]=primary_keys[str(sharedcol)].set_value(j[0], value=j[6])
		primary_keys[str(samecol)]=primary_keys[str(samecol)].set_value(j[0], value=j[7])
		print j[0]
		sys.stdout.flush()





primary_keys.to_csv(pwd+"final_fst_dxy_FIXED_db", )

