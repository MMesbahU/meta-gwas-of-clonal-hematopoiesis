#!/bin/bash


## 
# /medpop/esp2/mesbah/tools/GWAMA_v2.2.2/GWAMA 

# /medpop/esp2/mesbah/tools/GWAMA_v2.2.2/GWAMA --filelist /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/02.Metal_MetaAnalysis/01.run_GWAMA/test_EURonly.txt -qt --name_marker SNPID --name_n N --name_ea Alt --name_nea Ref --name_eaf AAF --name_beta REGENIE_BETA --name_se REGENIE_SE --output /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/02.Metal_MetaAnalysis/01.run_GWAMA/testGwama

## 
# /medpop/esp2/mesbah/tools/GWAMA_v2.2.2/GWAMA -i /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/02.Metal_MetaAnalysis/01.run_GWAMA/input.meta4.chr1.eur.txt -qt --name_marker SNP --name_n N --name_ea ALT --name_nea REF --name_eaf AAF --name_beta Beta --name_se SE --indel_alleles -o /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/02.Metal_MetaAnalysis/01.run_GWAMA/output.meta4.chr1.eur

## 
/medpop/esp2/mesbah/tools/MR-MEGA_v0.2/MR-MEGA  -i /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/02.Metal_MetaAnalysis/01.run_GWAMA/mrmega.input -o /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/02.Metal_MetaAnalysis/01.run_GWAMA/output.MR_mega --qt --pc 1 --no_std_names --name_pos POS --name_chr CHR --name_n N --name_se SE --name_beta Beta --name_eaf AAF --name_ea ALT --name_nea REF --name_marker SNP

