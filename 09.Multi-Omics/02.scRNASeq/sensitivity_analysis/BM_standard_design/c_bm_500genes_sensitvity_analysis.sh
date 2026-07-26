#!/bin/bash
#SBATCH --job-name=bm_gene_sub
#SBATCH --mem=64G
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=1
#SBATCH --time=6:00:00
#SBATCH --output=logs/bm_500genesub_%A_%a.log
#SBATCH --array=1-11

module load conda
source /apps/software/Miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate scdrs

h5ad_file="/data/cteu/2_users/zyu_lab/linke/chip_gwas_rev/Data/sc_annot_files/BM_standard_design.h5ad"
cov_file="/data/cteu/2_users/zyu_lab/linke/chip_gwas_rev/Data/sc_annot_files_qc/BM_standard_design_covs.tsv"
gs_file="/data/cteu/2_users/zyu_lab/linke/chip_gwas_rev/Data/scDRS/munged/sensitivity_analysis/gene_subsample_batches_genesymbol/batch$(printf '%02d' ${SLURM_ARRAY_TASK_ID}).gs"
out_folder="/data/cteu/2_users/zyu_lab/linke/chip_gwas_rev/Results/sensitivity_analysis/bm_500genes_sub/batch${SLURM_ARRAY_TASK_ID}"
mkdir -p "$out_folder" "${out_folder}/downstream"

echo "=== Gene subsample batch ${SLURM_ARRAY_TASK_ID} ==="
echo "gs: $gs_file"

# scdrs compute-score \
#     --h5ad-file    "$h5ad_file" \
#     --h5ad-species human \
#     --cov-file     "$cov_file" \
#     --gs-file      "$gs_file" \
#     --gs-species   human \
#     --out-folder   "$out_folder" \
#     --flag-filter-data True \
#     --flag-raw-count   True \
#     --n-ctrl 1000 \
#     --flag-return-ctrl-raw-score  False \
#     --flag-return-ctrl-norm-score True \
#     2>&1 | tee "${out_folder}/log.txt"

score_count=$(ls "${out_folder}"/*.full_score.gz 2>/dev/null | wc -l)

scdrs perform-downstream \
    --h5ad-file      "$h5ad_file" \
    --score-file     "${out_folder}/@.full_score.gz" \
    --group-analysis "anno" \
    --out-folder     "${out_folder}/downstream" \
    2>&1 | tee "${out_folder}/downstream/log.txt"

echo "=== Done batch ${SLURM_ARRAY_TASK_ID} ==="