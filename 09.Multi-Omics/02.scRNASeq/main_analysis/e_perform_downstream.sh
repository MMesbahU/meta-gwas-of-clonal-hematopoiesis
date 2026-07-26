#!/bin/bash
#SBATCH --job-name=chip_gwas_scDRS_downstream
#SBATCH --chdir=/data/cteu/2_users/zyu_lab/linke/chip_gwas_rev/
#SBATCH --mem=64G
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=4
#SBATCH --time=3:00:00
#SBATCH --output=/data/cteu/2_users/zyu_lab/linke/chip_gwas_rev/logs/downstream_%A_%a.log
#SBATCH --error=/data/cteu/2_users/zyu_lab/linke/chip_gwas_rev/logs/downstream_%A_%a.log
#SBATCH --array=3

module load conda
source /apps/software/Miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate scdrs

# Paths
manifest="/data/cteu/2_users/zyu_lab/linke/chip_gwas_rev/Data/sc_annot_files/file_list.txt"
results_base="/data/cteu/2_users/zyu_lab/linke/chip_gwas_rev/Results"

# Parse manifest — 3 columns now
h5ad_file=$(awk  -F'\t' -v line="${SLURM_ARRAY_TASK_ID}" 'NR==line {print $1}' "$manifest")
cov_file=$(awk   -F'\t' -v line="${SLURM_ARRAY_TASK_ID}" 'NR==line {print $2}' "$manifest")
cell_type=$(awk  -F'\t' -v line="${SLURM_ARRAY_TASK_ID}" 'NR==line {print $3}' "$manifest")

# Derive paths
sample_name=$(basename "$h5ad_file" .h5ad)
out_folder="${results_base}/${sample_name}"
downstream_folder="${out_folder}/downstream"
mkdir -p "$downstream_folder"

echo "=== Task $SLURM_ARRAY_TASK_ID: downstream for $sample_name ==="
echo "h5ad:       $h5ad_file"
echo "species:    human"
echo "cell_type:  $cell_type"
echo "score dir:  $out_folder"
echo "output dir: $downstream_folder"

# Sanity checks
if [ ! -f "$h5ad_file" ]; then
    echo "ERROR: h5ad not found: $h5ad_file"
    exit 1
fi

if [ -z "$cell_type" ]; then
    echo "ERROR: cell_type column is empty for line $SLURM_ARRAY_TASK_ID"
    echo "Check column 4 of your manifest"
    exit 1
fi

# Check that compute-score finished — need at least one score file
score_count=$(ls "${out_folder}"/*.full_score.gz 2>/dev/null | wc -l)
if [ "$score_count" -eq 0 ]; then
    echo "ERROR: no .full_score.gz files found in $out_folder"
    echo "Did compute-score finish successfully for $sample_name?"
    exit 1
fi
echo "Found $score_count score file(s) — proceeding"

# List the score files being used
ls "${out_folder}"/*.full_score.gz

# Run downstream analysis
scdrs perform-downstream \
    --h5ad-file      "$h5ad_file" \
    --score-file     "${out_folder}/@.full_score.gz" \
    --group-analysis "$cell_type" \
    --flag_raw_count True \
    --out-folder     "$downstream_folder" \
    2>&1 | tee "${downstream_folder}/log.txt"

# Check if downstream produced output
out_count=$(ls "${downstream_folder}"/*.csv 2>/dev/null | wc -l)
if [ "$out_count" -eq 0 ]; then
    echo "WARNING: no .csv output files found in $downstream_folder"
    echo "Check ${downstream_folder}/log.txt for errors"
    exit 1
fi

echo "=== Done: $sample_name — $out_count output files written ==="