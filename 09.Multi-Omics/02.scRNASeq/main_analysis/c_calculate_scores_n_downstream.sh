#!/bin/bash
#SBATCH --job-name=chip_gwas_scDRS
#SBATCH --chdir=/data/cteu/2_users/zyu_lab/linke/chip_gwas_rev/
#SBATCH --mem=64G
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00
#SBATCH --output=/data/cteu/2_users/zyu_lab/linke/chip_gwas_rev/logs/chip_gwas_scDRS_%A_%a.log
#SBATCH --error=/data/cteu/2_users/zyu_lab/linke/chip_gwas_rev/logs/chip_gwas_scDRS_%A_%a.log
#SBATCH --array=2

module load conda
source /apps/software/Miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate scdrs

# Paths
manifest="/data/cteu/2_users/zyu_lab/linke/chip_gwas_rev/Data/sc_annot_files/file_list.txt"
gs_file="/data/cteu/2_users/zyu_lab/linke/chip_gwas_rev/Data/scDRS/munged/chip_dnmt3a_tet2_munged_n1000.gs"
results_base="/data/cteu/2_users/zyu_lab/linke/chip_gwas_rev/Results"

# Parse manifest line for this task
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$manifest")
h5ad_file=$(echo "$line" | awk -F'\t' '{print $1}')
cov_file=$(echo "$line" | awk -F'\t' '{print $2}')

# Derive output folder from h5ad filename
sample_name=$(basename "$h5ad_file" .h5ad)
out_folder="${results_base}/${sample_name}"
mkdir -p "$out_folder"

# Sanity checks before running
if [ ! -f "$h5ad_file" ]; then
    echo "ERROR: h5ad file not found: $h5ad_file"
    exit 1
fi
if [ ! -f "$cov_file" ]; then
    echo "ERROR: cov file not found: $cov_file"
    exit 1
fi

echo "=== Task $SLURM_ARRAY_TASK_ID ==="
echo "h5ad:    $h5ad_file"
echo "species: human"
echo "cov:     $cov_file"
echo "out:     $out_folder"

scdrs compute-score \
    --h5ad-file    "$h5ad_file" \
    --h5ad-species "human" \
    --cov-file     "$cov_file" \
    --gs-file      "$gs_file" \
    --gs-species   human \
    --out-folder   "$out_folder" \
    --flag-filter-data True \
    --flag-raw-count True \
    --n-ctrl 1000 \
    --flag-return-ctrl-raw-score False \
    --flag-return-ctrl-norm-score True \
    2>&1 | tee "${out_folder}/log.txt"

echo "=== Done: $sample_name ==="


