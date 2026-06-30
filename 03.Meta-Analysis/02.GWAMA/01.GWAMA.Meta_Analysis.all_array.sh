#!/bin/bash

source /broad/software/scripts/useuse

# use Tabix
################
# MCPS vs all
# qsub -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/tmpdir -t 1-$(wc -l /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/input_files/meta5_vs_MCPSlist_of_21summary.list |awk '{print $1}') -pe smp 4 -binding linear:4 -l h_vmem=30G -l h_rt=10:00:00 -N mcps_vs_others /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/01.GWAMA.Meta_Analysis.all_array.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/input_files/meta5_vs_MCPSlist_of_21summary.list /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out /medpop/esp2/mesbah/tools/GWAMA_v2.2.2/GWAMA
####
###### qsub -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/tmpdir -t 1-$(wc -l /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/list_of_96_gwama_inputfiles.list |awk '{print $1}') -pe smp 4 -binding linear:4 -l h_vmem=30G -l h_rt=10:00:00 -N gwama_metaAnalyses /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/01.GWAMA.Meta_Analysis.all_array.sh /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/list_of_96_gwama_inputfiles.list /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out /medpop/esp2/mesbah/tools/GWAMA_v2.2.2/GWAMA

###########################
GWAMA_file_list=${1} # /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/list_of_96_gwama_inputfiles.list

outDir=${2}

## Run GWAMA
GWAMA=${3} # /medpop/esp2/mesbah/tools/GWAMA_v2.2.2/GWAMA
## input files
gwas_list=$(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${GWAMA_file_list} )

output_prefix=${outDir}/GWAMA.chr1_22.$(basename $(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${GWAMA_file_list} ) ".txt")
############################
# Ensure input file exists
if [[ ! -f "${gwas_list}" ]]; then
    echo "Error: Input file '${gwas_list}' not found."
    exit 1
fi

# 2. Process each line
while read -r file; do
  # Skip empty lines (if any)
  [[ -z "$file" ]] && continue

  # 2A. Ensure the file actually exists
  if [[ ! -f "$file" ]]; then
    echo "Error: File '$file' not found."
    exit 1
  fi
		    
  # 2B. Check if the file ends in .gz
  if [[ "$file" == *.gz ]]; then
    # Attempt to decompress
    if ! gzip -d "$file"; then
       echo "Error decompressing '$file'."
       exit 1
    fi
    echo "Decompressed: $file"
  else
    # It's not gzipped, so just "process" it (example echo)
    echo "File '$file' is not compressed; proceeding with processing..."
    # put your code/commands for uncompressed files here
  fi
done < "$gwas_list"
##### 2025
# Maxico City Preopective Study (MCPS)
# /broad/hptmp/mesbah/dataset/ch_gwas/mcps/gwama.DNMT3A.GCST90435342.tsv
# /broad/hptmp/mesbah/dataset/ch_gwas/mcps/gwama.ASXL1.GCST90435344.tsv
# /broad/hptmp/mesbah/dataset/ch_gwas/mcps/gwama.PPM1D.GCST90435345.tsv
# /broad/hptmp/mesbah/dataset/ch_gwas/mcps/gwama.TET2.GCST90435343.tsv
# /broad/hptmp/mesbah/dataset/ch_gwas/mcps/gwama.CH.GCST90435341.tsv
# /broad/hptmp/mesbah/dataset/ch_gwas/mcps/gwama.SF3B1_SRSF2.GCST90435352.tsv
# /broad/hptmp/mesbah/dataset/ch_gwas/mcps/gwama.TP53.GCST90435346.tsv

# CHIP|DNMT3A|TET2: while read ANC; do while read pheno; do rm /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.ukbb200k_ukb250k_AoU250k_TOPMed72k_MGBB53k_BioVU54k_MCPS136k.txt; for files in $(echo -e "/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/chr1_22.MultiANC.ukb200k_N193342.has${pheno}.regenie.noNA.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/chr1_22.MultiANC.ukb250k_N243350.has${pheno}.regenie.noNA.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/aou/chr1_22.${ANC}.AoU_v71.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/topmed/GWAMA.${ANC}.TOPMed74k.chr1_22.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mgbb/has${pheno}.chr1_22.${ANC}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/biovu/BioVU.saige_has${pheno}_results_merged_subset.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/mcps/gwama.${pheno}.GCST*.tsv"); do ls ${files} | awk '{print $NF}' >> /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.ukbb200k_ukb250k_AoU250k_TOPMed72k_MGBB53k_BioVU54k_MCPS136k.txt; done; done < <(echo -e "CH\nDNMT3A\nTET2"); done < <(echo -e "MultiANC")

# ASXL1: 
# while read ANC; do while read pheno; do rm /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.ukbb200k_ukb250k_AoU250k_TOPMed72k_MGBB53k_MCPS136k.txt; for files in $(echo -e "/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/chr1_22.MultiANC.ukb200k_N193342.has${pheno}.regenie.noNA.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/chr1_22.MultiANC.ukb250k_N243350.has${pheno}.regenie.noNA.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/aou/chr1_22.${ANC}.AoU_v71.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/topmed/GWAMA.${ANC}.TOPMed74k.chr1_22.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mgbb/has${pheno}.chr1_22.${ANC}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/mcps/gwama.${pheno}.GCST*.tsv"); do ls ${files} | awk '{print $NF}' >> /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.ukbb200k_ukb250k_AoU250k_TOPMed72k_MGBB53k_MCPS136k.txt; done; done < <(echo -e "ASXL1"); done < <(echo -e "MultiANC")

#SF: while read ANC; do while read pheno; do rm /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.ukbb200k_ukb250k_AoU250k_TOPMed72k_MGBB53k_MCPS136k.txt; for files in $(echo -e "/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/chr1_22.MultiANC.ukb200k_N193342.has${pheno}.regenie.noNA.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/chr1_22.MultiANC.ukb250k_N243350.has${pheno}.regenie.noNA.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/aou/chr1_22.${ANC}.AoU_v71.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/topmed/GWAMA.${ANC}.TOPMed74k.chr1_22.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mgbb/has${pheno}.chr1_22.${ANC}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/mcps/gwama.SF3B1_SRSF2.GCST90435352.tsv"); do ls $files | awk '{print $NF}' >> /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.ukbb200k_ukb250k_AoU250k_TOPMed72k_MGBB53k_MCPS136k.txt; done; done < <(echo -e "SF"); done < <(echo -e "MultiANC")

#DDR:TP53: while read ANC; do while read pheno; do rm /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.ukbb200k_ukb250k_AoU250k_TOPMed72k_MGBB53k_tp53MCPS136k.txt;for files in $(echo -e "/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/chr1_22.MultiANC.ukb200k_N193342.has${pheno}.regenie.noNA.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/chr1_22.MultiANC.ukb250k_N243350.has${pheno}.regenie.noNA.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/aou/chr1_22.${ANC}.AoU_v71.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/topmed/GWAMA.${ANC}.TOPMed74k.chr1_22.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mgbb/has${pheno}.chr1_22.${ANC}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/mcps/gwama.TP53.GCST*.tsv"); do ls $files | awk '{print $NF}' >> /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.ukbb200k_ukb250k_AoU250k_TOPMed72k_MGBB53k_tp53MCPS136k.txt; done; done < <(echo -e "DDR"); done < <(echo -e "MultiANC")
#DDR PPM1D: while read ANC; do while read pheno; do rm /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.ukbb200k_ukb250k_AoU250k_TOPMed72k_MGBB53k_ppm1dMCPS136k.txt; for files in $(echo -e "/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/chr1_22.MultiANC.ukb200k_N193342.has${pheno}.regenie.noNA.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/chr1_22.MultiANC.ukb250k_N243350.has${pheno}.regenie.noNA.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/aou/chr1_22.${ANC}.AoU_v71.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/topmed/GWAMA.${ANC}.TOPMed74k.chr1_22.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mgbb/has${pheno}.chr1_22.${ANC}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/mcps/gwama.PPM1D.GCST*.tsv"); do ls $files | awk '{print $NF}' >> /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.ukbb200k_ukb250k_AoU250k_TOPMed72k_MGBB53k_ppm1dMCPS136k.txt; done; done < <(echo -e "DDR"); done < <(echo -e "MultiANC")
## ANC-specific
# AFR/AMR: while read pheno; do while read ANC; do rm /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.ukbb450k_AoU250k_TOPMed72k_MGBB53k_MCPS136k.txt; for files in $(echo -e "/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/has${pheno}.UKB450k.chr1_22.${ANC}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/aou/chr1_22.${ANC}.AoU_v71.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/topmed/GWAMA.${ANC}.TOPMed74k.chr1_22.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mgbb/has${pheno}.chr1_22.${ANC}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/mcps/gwama.${pheno}.GCST*.tsv"); do ls $files | awk '{print $NF}' >> /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.ukbb450k_AoU250k_TOPMed72k_MGBB53k_MCPS136k.txt; done; done < <(echo -e "AFR\nAMR"); done < <(echo -e "CH\nDNMT3A\nTET2\nASXL1")

# SF: while read pheno; do while read ANC; do rm /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.ukbb450k_AoU250k_TOPMed72k_MGBB53k_MCPS136k.txt; for files in $(echo -e "/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/has${pheno}.UKB450k.chr1_22.${ANC}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/aou/chr1_22.${ANC}.AoU_v71.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/topmed/GWAMA.${ANC}.TOPMed74k.chr1_22.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mgbb/has${pheno}.chr1_22.${ANC}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/mcps/gwama.SF3B1_SRSF2.GCST90435352.tsv"); do ls $files | awk '{print $NF}' >> /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.ukbb450k_AoU250k_TOPMed72k_MGBB53k_MCPS136k.txt;done; done < <(echo -e "AFR\nAMR"); done < <(echo -e "SF")

# DDR: while read pheno; do while read ANC; do rm /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.ukbb450k_AoU250k_TOPMed72k_MGBB53k_tp53MCPS136k.txt; for files in $(echo -e "/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/has${pheno}.UKB450k.chr1_22.${ANC}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/aou/chr1_22.${ANC}.AoU_v71.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/topmed/GWAMA.${ANC}.TOPMed74k.chr1_22.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mgbb/has${pheno}.chr1_22.${ANC}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/mcps/gwama.TP53.GCST*.tsv"); do ls ${files} | awk '{print $NF}' >> /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.ukbb450k_AoU250k_TOPMed72k_MGBB53k_tp53MCPS136k.txt;done; done < <(echo -e "AFR\nAMR"); done < <(echo -e "DDR")

# DDR: while read pheno; do while read ANC; do rm /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.ukbb450k_AoU250k_TOPMed72k_MGBB53k_ppm1dMCPS136k.txt; for files in $(echo -e "/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/has${pheno}.UKB450k.chr1_22.${ANC}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/aou/chr1_22.${ANC}.AoU_v71.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/topmed/GWAMA.${ANC}.TOPMed74k.chr1_22.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mgbb/has${pheno}.chr1_22.${ANC}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/mcps/gwama.PPM1D.GCST*.tsv"); do ls ${files} | awk '{print $NF}' >> /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.ukbb450k_AoU250k_TOPMed72k_MGBB53k_ppm1dMCPS136k.txt;done; done < <(echo -e "AFR\nAMR"); done < <(echo -e "DDR")

# list : ls /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/*MCPS136k.txt | awk '{print $NF}' > /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/input_files/meta5_vs_MCPSlist_of_21summary.list

###################### MultiANC Meta-GWAS ####################
## MultiANC v1: 
# CH | DNMT3A | TET2 | CHvaf10 | ASXL1 | SF | DDR : UKB450k + AoU + TOPMed + MGBB
# while read ANC; do while read pheno; do echo -e "/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/chr1_22.has${pheno}.${ANC}.ukbb450k.ukb200k_N193342.ukb250k_N243350.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/aou/chr1_22.${ANC}.AoU_v71.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/topmed/GWAMA.${ANC}.TOPMed74k.chr1_22.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mgbb/has${pheno}.chr1_22.${ANC}.regenie.tsv" > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.ukbb450k_AoU250k_TOPMed72k_MGBB53k.txt; qsub -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/tmpdir -pe smp 4 -binding linear:4 -l h_vmem=30G -l h_rt=10:00:00 -N ${ANC}_${pheno}.ukb_aou_top_mgb /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/01.GWAMA.Meta_Analysis.all.sh /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.ukbb450k_AoU250k_TOPMed72k_MGBB53k.txt /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/chr1_22.has${pheno}.${ANC}.UKBB450k_AoU250k_TOPMed72k_MGBB53k; done < <(echo -e "CH\nCHvaf10\nDNMT3A\nTET2\nASXL1\nSF\nDDR"); done < <(echo -e "MultiANC")

# CH | DNMT3A | TET2: UKB450k + AoU + TOPMed + MGBB + BioVU
# while read ANC; do while read pheno; do echo -e "/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/chr1_22.has${pheno}.${ANC}.ukbb450k.ukb200k_N193342.ukb250k_N243350.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/aou/chr1_22.${ANC}.AoU_v71.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/topmed/GWAMA.${ANC}.TOPMed74k.chr1_22.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mgbb/has${pheno}.chr1_22.${ANC}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/biovu/BioVU.saige_has${pheno}_results_merged_subset.tsv" > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.ukbb450k_AoU250k_TOPMed72k_MGBB53k_BioVU54k.txt; qsub -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/tmpdir -pe smp 4 -binding linear:4 -l h_vmem=30G -l h_rt=10:00:00 -N ${ANC}_${pheno}.ukb_aou_top_mgb_biovu /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/01.GWAMA.Meta_Analysis.all.sh /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.ukbb450k_AoU250k_TOPMed72k_MGBB53k_BioVU54k.txt /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/chr1_22.has${pheno}.${ANC}.UKBB450k_AoU250k_TOPMed72k_MGBB53k_BioVU54k; done < <(echo -e "CH\nDNMT3A\nTET2"); done < <(echo -e "MultiANC")

#################################################################
## Input for GWAMA
# while read pheno; do while read ANC; do echo -e "/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/has${pheno}.UKB450k.chr1_22.${ANC}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/aou/chr1_22.${ANC}.AoU_v71.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/topmed/GWAMA.${ANC}.TOPMed74k.chr1_22.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mgbb/has${pheno}.chr1_22.${ANC}.regenie.tsv" > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.ukbb450k_AoU250k_TOPMed72k_MGBB53k.txt; done < <(echo -e "AFR\nAMR\nEUR\nFemale\nMale"); done < <(echo -e "CH\nCHvaf10\nDNMT3A\nTET2\nASXL1\nSF\nDDR")

## EUR with bioVU
# while read pheno; do while read ANC; do echo -e "/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/has${pheno}.UKB450k.chr1_22.${ANC}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/aou/chr1_22.${ANC}.AoU_v71.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/topmed/GWAMA.${ANC}.TOPMed74k.chr1_22.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mgbb/has${pheno}.chr1_22.${ANC}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/biovu/BioVU.saige_has${pheno}_results_merged_subset.tsv" > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.ukbb450k_AoU250k_TOPMed72k_MGBB53k_BioVU54k.txt; done < <(echo -e "EUR"); done < <(echo -e "CH\nDNMT3A\nTET2")
## MultiANC UKB 200k 250k
# while read ANC; do while read pheno; do echo -e "/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/chr1_22.MultiANC.ukb200k_N193342.has${pheno}.regenie.noNA.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/chr1_22.MultiANC.ukb250k_N243350.has${pheno}.regenie.noNA.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/aou/chr1_22.${ANC}.AoU_v71.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/topmed/GWAMA.${ANC}.TOPMed74k.chr1_22.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mgbb/has${pheno}.chr1_22.${ANC}.regenie.tsv" > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.ukbb200k_ukb250k_AoU250k_TOPMed72k_MGBB53k.txt; done < <(echo -e "CH\nCHvaf10\nDNMT3A\nTET2\nASXL1\nSF\nDDR"); done < <(echo -e "MultiANC")
#  MultiANC UKB 200k 250k with BioVU
# while read ANC; do while read pheno; do echo -e "/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/chr1_22.MultiANC.ukb200k_N193342.has${pheno}.regenie.noNA.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/chr1_22.MultiANC.ukb250k_N243350.has${pheno}.regenie.noNA.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/aou/chr1_22.${ANC}.AoU_v71.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/topmed/GWAMA.${ANC}.TOPMed74k.chr1_22.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mgbb/has${pheno}.chr1_22.${ANC}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/biovu/BioVU.saige_has${pheno}_results_merged_subset.tsv" > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.ukbb200k_ukb250k_AoU250k_TOPMed72k_MGBB53k_BioVU54k.txt; done < <(echo -e "CH\nDNMT3A\nTET2"); done < <(echo -e "MultiANC")
## No UKBB
# while read pheno; do while read ANC; do echo -e "/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/aou/chr1_22.${ANC}.AoU_v71.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/topmed/GWAMA.${ANC}.TOPMed74k.chr1_22.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mgbb/has${pheno}.chr1_22.${ANC}.regenie.tsv" > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.AoU250k_TOPMed72k_MGBB53k.noukbb.txt; done < <(echo -e "AFR\nAMR\nEUR\nFemale\nMale\nMultiANC"); done < <(echo -e "CH\nCHvaf10\nDNMT3A\nTET2\nASXL1\nSF\nDDR")

## EUR with bioVU
# while read pheno; do while read ANC; do echo -e "/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/aou/chr1_22.${ANC}.AoU_v71.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/topmed/GWAMA.${ANC}.TOPMed74k.chr1_22.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mgbb/has${pheno}.chr1_22.${ANC}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/biovu/BioVU.saige_has${pheno}_results_merged_subset.tsv" > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.AoU250k_TOPMed72k_MGBB53k_BioVU54k.noukbb.txt; done < <(echo -e "EUR"); done < <(echo -e "CH\nDNMT3A\nTET2")

# while read pheno; do while read ANC; do echo -e "/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/aou/chr1_22.${ANC}.AoU_v71.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/topmed/GWAMA.${ANC}.TOPMed74k.chr1_22.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mgbb/has${pheno}.chr1_22.${ANC}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/biovu/BioVU.saige_has${pheno}_results_merged_subset.tsv" > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.AoU250k_TOPMed72k_MGBB53k_BioVU54k.noukbb.txt; done < <(echo -e "MultiANC"); done < <(echo -e "CH\nDNMT3A\nTET2")

## ls /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/*.txt | awk 'NR>17{print $NF}' > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/list_of_96_gwma_inputfiles.list

###################################################################
###################################################################


#################################################################

############### Task: Worked
## 1. Run Meta-Analysis
###########################

######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)
##################################
# zcat /medpop/esp/mesbah/GWAS_CHIP/inputUKB/eurUKB/chr1_22.has${traits}.21Aug2021_ukbEUR.tsv.gz | awk '(NR==1){print "SNPID\tN\tREF\tALT\tAAF\tBETA\tSE\tP"}(NR>1 && $2~/^[0-9]+/){print $1"\t"($14+$18)"\t"$4"\t"$5"\t"$13"\t"$22"\t"$23"\t"$12}' | awk '!seen[$1]++'| gzip -c > /medpop/esp/mesbah/GWAS_CHIP/inputUKB/eurUKB/chr1_22.ukb200k.hg37_eur_${traits}.tsv.gz
##
## 
# gwas_list=${1}

# output_prefix=${2}

# outDir=${3} # summary tsv file for plotting

## Run GWAMA
# GWAMA=/medpop/esp2/mesbah/tools/GWAMA_v2.2.2/GWAMA 
${GWAMA} \
	-i ${gwas_list} \
	-qt \
	--name_marker SNPID \
	--name_n N \
	--name_ea ALT \
	--name_nea REF \
	--name_eaf AAF \
	--name_beta BETA \
	--name_se SE \
	--indel_alleles \
	-o ${output_prefix}

## compress
# gzip ${output_prefix}.out
gzip ${output_prefix}.log.out
# Recompress the files after processing
# while read -r lines; do
  #      if ! gzip -f "${lines}"; then
   #             echo "Error recompressing: ${lines}"
    #            exit 1
     #   fi
      #  echo "Recompressed: ${lines}"
# done < "${gwas_list}"
# Decompress the GWAMA output file
if ! gzip -f "${output_prefix}.out"; then
    echo "Error compressing: ${output_prefix}.out"
    exit 1
fi
echo "Compressed output file: ${output_prefix}.out.gz"
###########


######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########

