#! /bin/bash -l
#SBATCH -A b2014034
#SBATCH -p node -C mem512GB
#SBATCH -n 16
#SBATCH -t 100:00:00
#SBATCH -J coverage_filter_PI
#SBATCH --mail-user venkat.talla@ebc.uu.se
#SBATCH --mail-type=ALL

#####rm -rf sinapis_juvernica
#####rm -rf juvernica_reali
#####rm -rf sinapis_reali
#####rm -rf sinapis
#####rm -rf reali
#####rm -rf juvernica
#####
#####echo "./lists/swe_sin_allele" >> sinapis_juvernica
#####echo "./lists/spanish_sinapis" >> sinapis_juvernica
#####echo "./lists/kaz_sin" >> sinapis_juvernica
#####echo "./lists/irish_juvernica" >> sinapis_juvernica
#####echo "./lists/kazak_juvernica" >> sinapis_juvernica
#####
#####
#####echo "./lists/reali.txt" >> juvernica_reali
#####echo "./lists/irish_juvernica" >> juvernica_reali
#####echo "./lists/kazak_juvernica" >> juvernica_reali
#####
#####
#####
#####echo "./lists/reali.txt" >> sinapis_reali
#####echo "./lists/swe_sin_allele" >> sinapis_reali
#####echo "./lists/spanish_sinapis" >> sinapis_reali
#####echo "./lists/kaz_sin" >> sinapis_reali
#####
#####
#####echo "./lists/swe_sin_allele" >> sinapis
#####echo "./lists/spanish_sinapis" >> sinapis
#####echo "./lists/kaz_sin" >> sinapis
#####
#####echo "./lists/reali.txt" >> reali
#####echo "./lists/irish_juvernica" >> juvernica
#####echo "./lists/kazak_juvernica" >> juvernica
#####echo "./lists/irish_juvernica" >> irish_juvernica
#####echo "./lists/kazak_juvernica" >> kazak_juvernica
#####echo "./lists/spanish_reali"   >> spanish_reali
#####echo "./lists/swe_sin_allele"  >> swe_sin_allele
#####echo "./lists/kaz_sin" 		   >> kaz_sin
#####echo "./lists/spanish_sinapis" >> spanish_sinapis



lists='/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/lists/'


for i in scaf1 scaf2 scaf3 scaf4 scaf5 scaf6 scaf7 scaf8 scaf9 scaf10
do
echo '/'$i'/'
python get_N_sites.py   './lists/irish_juvernica'  './lists/kaz_sin'         './irish_juvernica_kaz_sin'         '/'$i'/' &
python get_N_sites.py   './lists/irish_juvernica'  './lists/kazak_juvernica' './irish_juvernica_kazak_juvernica' '/'$i'/' &
python get_N_sites.py   './lists/irish_juvernica'  './lists/spanish_reali'   './irish_juvernica_spanish_reali'   '/'$i'/' &
python get_N_sites.py   './lists/irish_juvernica'  './lists/spanish_sinapis' './irish_juvernica_spanish_sinapis' '/'$i'/' &
python get_N_sites.py   './lists/irish_juvernica'  './lists/swe_sin_allele'  './irish_juvernica_swe_sin_allele'  '/'$i'/' &
python get_N_sites.py   './lists/kaz_sin'          './lists/spanish_reali'   './kaz_sin_spanish_reali'           '/'$i'/' &
python get_N_sites.py   './lists/kaz_sin'          './lists/spanish_sinapis' './kaz_sin_spanish_sinapis'         '/'$i'/' &
python get_N_sites.py   './lists/kaz_sin'          './lists/swe_sin_allele'  './kaz_sin_swe_sin_allele'          '/'$i'/' &
wait
python get_N_sites.py   './lists/kazak_juvernica'  './lists/kaz_sin'         './kazak_juvernica_kaz_sin'         '/'$i'/' &
python get_N_sites.py   './lists/kazak_juvernica'  './lists/spanish_reali'   './kazak_juvernica_spanish_reali'   '/'$i'/' &
python get_N_sites.py   './lists/kazak_juvernica'  './lists/spanish_sinapis' './kazak_juvernica_spanish_sinapis' '/'$i'/' &
python get_N_sites.py   './lists/kazak_juvernica'  './lists/swe_sin_allele'  './kazak_juvernica_swe_sin_allele'  '/'$i'/' &
python get_N_sites.py   './lists/spanish_reali'    './lists/spanish_sinapis' './spanish_reali_spanish_sinapis'   '/'$i'/' &
python get_N_sites.py   './lists/spanish_reali'    './lists/swe_sin_allele'  './spanish_reali_swe_sin_allele'    '/'$i'/' &
python get_N_sites.py   './lists/spanish_sinapis'  './lists/swe_sin_allele'  './spanish_sinapis_swe_sin_allele'  '/'$i'/' &
wait
python get_N_sites_V2.0.py -i sinapis_juvernica -o  sinapis_juvernica -p $i &
python get_N_sites_V2.0.py -i juvernica_reali   -o  juvernica_reali   -p $i &
python get_N_sites_V2.0.py -i sinapis_reali     -o  sinapis_reali     -p $i &
python get_N_sites_V2.0.py -i sinapis           -o  sinapis        	  -p $i &
python get_N_sites_V2.0.py -i reali             -o  reali        	  -p $i &
python get_N_sites_V2.0.py -i juvernica         -o  juvernica 		  -p $i &
wait
python get_N_sites_V2.0.py -i irish_juvernica -o irish_juvernica  -p $i &
python get_N_sites_V2.0.py -i kazak_juvernica -o kazak_juvernica  -p $i &
python get_N_sites_V2.0.py -i spanish_reali   -o spanish_reali    -p $i &
python get_N_sites_V2.0.py -i swe_sin_allele  -o swe_sin_allele   -p $i &
python get_N_sites_V2.0.py -i kaz_sin 		  -o kaz_sin          -p $i &
python get_N_sites_V2.0.py -i spanish_sinapis -o spanish_sinapis  -p $i &
wait
done



#rm -rf scaf*


