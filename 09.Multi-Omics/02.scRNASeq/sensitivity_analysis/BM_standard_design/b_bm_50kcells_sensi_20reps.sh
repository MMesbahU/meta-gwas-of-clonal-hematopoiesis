#!/bin/bash
#SBATCH --job-name=bm_cell_sub
#SBATCH --mem=64G
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --output=/data/cteu/2_users/zyu_lab/linke/chip_gwas_rev/logs/bm_50kcells_sensi_%A_%a.log
#SBATCH --array=20

module load conda
source /apps/software/Miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate scdrs

manifest="/data/cteu/2_users/zyu_lab/linke/chip_gwas_rev/Data/sc_annot_files_subsets/helper_files/bm_cell_subsample_manifest.txt"
gs_file="/data/cteu/2_users/zyu_lab/linke/chip_gwas_rev/Data/scDRS/munged/chip_dnmt3a_tet2_munged_n1000.gs"

h5ad_file=$(awk -F'\t' -v line="${SLURM_ARRAY_TASK_ID}" 'NR==line {print $1}' "$manifest")
cov_file=$(awk  -F'\t' -v line="${SLURM_ARRAY_TASK_ID}" 'NR==line {print $2}' "$manifest")
cell_type=$(awk -F'\t' -v line="${SLURM_ARRAY_TASK_ID}" 'NR==line {print $3}' "$manifest")

rep=$((SLURM_ARRAY_TASK_ID - 1))
out_folder="/data/cteu/2_users/zyu_lab/linke/chip_gwas_rev/Results/sensitivity_analysis/bm_50kcells_sub/rep${rep}"
mkdir -p "$out_folder" "${out_folder}/downstream"

echo "Rep $rep | h5ad: $h5ad_file | cov: $cov_file"

scdrs compute-score \
    --h5ad-file    "$h5ad_file" \
    --h5ad-species human \
    --cov-file     "$cov_file" \
    --gs-file      "$gs_file" \
    --gs-species   human \
    --out-folder   "$out_folder" \
    --flag-filter-data True \
    --flag-raw-count   True \
    --n-ctrl 1000 \
    --flag-return-ctrl-raw-score  False \
    --flag-return-ctrl-norm-score True \
    2>&1 | tee "${out_folder}/log.txt"

score_count=$(ls "${out_folder}"/*.full_score.gz 2>/dev/null | wc -l)
[ "$score_count" -eq 0 ] && echo "ERROR: no score files" && exit 1

scdrs perform-downstream \
    --h5ad-file      "$h5ad_file" \
    --score-file     "${out_folder}/@.full_score.gz" \
    --group-analysis "$cell_type" \
    --out-folder     "${out_folder}/downstream" \
    2>&1 | tee "${out_folder}/downstream/log.txt"

echo "Done rep $rep"