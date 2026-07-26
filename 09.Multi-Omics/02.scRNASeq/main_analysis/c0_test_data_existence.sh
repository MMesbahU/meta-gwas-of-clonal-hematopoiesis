manifest="/data/cteu/2_users/zyu_lab/linke/chip_gwas_rev/Data/sc_annot_files/file_list.txt"
n_lines=$(wc -l < "$manifest")

echo "Testing all $n_lines manifest lines..."
echo "=================================="

for i in $(seq 1 $n_lines); do
    h5ad_file=$(awk -F'\t' -v line="$i" 'NR==line {print $1}' "$manifest")
    cov_file=$(awk -F'\t' -v line="$i" 'NR==line {print $2}' "$manifest")

    echo "Line $i:"
    echo "  h5ad:    $h5ad_file"
    echo "  cov:     $cov_file"

    # Check files exist
    [ ! -f "$h5ad_file" ] && echo "  ERROR: h5ad MISSING" || echo "  OK: h5ad exists"
    [ ! -f "$cov_file"  ] && echo "  ERROR: cov  MISSING" || echo "  OK: cov exists"
    echo "----------------------------------"
done