# CHIP inclusive
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90165001-GCST90166000/GCST90165267/GCST90165267_buildGRCh38.tsv.gz .

# CHIP exclusive
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90165001-GCST90166000/GCST90165261/GCST90165261_buildGRCh38.tsv.gz . 1>>myrgc.log 2>>myrgc.err &

# AFR: 375 African ancestry cases, 8,177 African ancestry controls
# Clonal haematopoiesis of indeterminate potential (with or without mosaic chromosomal alteration)
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90165001-GCST90166000/GCST90165263/GCST90165263_buildGRCh38.tsv.gz . 1>>myrgc.log 2>>myrgc.err &

# TET2: Clonal haematopoiesis of indeterminate potential (TET2 mutation)
# 3,918 European ancestry cases, 342,869 European ancestry controls
# 1,888 European ancestry cases, 135,106 European ancestry controls
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90165001-GCST90166000/GCST90165281/GCST90165281_buildGRCh38.tsv.gz . 1>>myrgc.log 2>>myrgc.err &

#DNMT3A 
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90165001-GCST90166000/GCST90165271/GCST90165271_buildGRCh38.tsv.gz .  1>>myrgc.log 2>>myrgc.err &

# ASXL1: Clonal haematopoiesis of indeterminate potential (ASXL1 mutation)
# 2,007 European ancestry cases, 342,869 European ancestry controls
# 1,007 European ancestry cases, 135,106 European ancestry controls
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90165001-GCST90166000/GCST90165259/GCST90165259_buildGRCh38.tsv.gz . 1>>myrgc.log 2>>myrgc.err &


## 
# prepare summary

for files in $(ls *_buildGRCh38.tsv.gz | awk '{print $NF}'); do  zcat ${files} | awk 'NR==1{print $0}(NR>1 && $12<=5e-6 && $13>=0.001 && $13<=0.999){print $0}' > rgc_p5e6.$(basename ${files} ".gz").tsv; done &
