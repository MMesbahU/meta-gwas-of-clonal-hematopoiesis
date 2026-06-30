#!/usr/bin/env python

import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
# import matplotlib.lines as mlines
# Read input and output paths from arguments
input_file = sys.argv[1]
output_plot_mean_3sd = sys.argv[2]
output_plot_median_3mad = sys.argv[3]
output_plot_iqr = sys.argv[4]
title_var = sys.argv[5]
EAF_threshold = float(sys.argv[6])
# X_percent = float(sys.argv[7])

# Load GWAS summary data
df = pd.read_csv(input_file, sep='\t', compression='gzip', usecols=['SNPID', 'EAF', 'Z', 'P', 'N'])

# Split the SNPID column into CHR and BP
df[['CHR', 'BP', 'REF', 'ALT']] = df['SNPID'].str.split(':', expand=True)
df['CHR'] = df['CHR'].str.replace('chr', '').astype(int)
df['BP'] = df['BP'].astype(int)

# Set the threshold for significance
threshold = 5e-8

# Keep EAF EAF_threshold

eaf_low_bound = EAF_threshold
eaf_high_bound = float( 1 - EAF_threshold)
df = df[(df['EAF'] >= eaf_low_bound) & (df['EAF'] <= eaf_high_bound)]

# Sort data by chromosome and position
df = df.sort_values(by=['CHR', 'BP']).reset_index(drop=True)
df['ind'] = range(len(df))  # For plotting

# Set color scheme for chromosomes
chromosome_colors = {i: 'skyblue' if i % 2 == 0 else 'grey' for i in sorted(df['CHR'].unique())}

# Define chromosome ticks for x-axis
ticks = [df[df['CHR'] == chrom]['ind'].mean() for chrom in sorted(df['CHR'].unique())]
labels = list(range(1, 23))

# Calculate statistics for filtering
mean_n = df['N'].mean()
std_n = df['N'].std()
median_n = df['N'].median()
mad_n = stats.median_abs_deviation(df['N'])
min_n = df['N'].min()
max_n = df['N'].max()
q1, q3 = df['N'].quantile([0.25, 0.75])
iqr_n = q3 - q1

# Define bounds for each filtering criterion
mean_3sd_lower, mean_3sd_upper = mean_n - 3 * std_n, mean_n + 3 * std_n
median_3mad_lower, median_3mad_upper = median_n - 3 * mad_n, median_n + 3 * mad_n
iqr_lower, iqr_upper = q1 - 1.5 * iqr_n, q3 + 1.5 * iqr_n

# Print N summary statistics
print(f"Mean of N: {mean_n}")
print(f"Median of N: {median_n}")
print(f"Min of N: {min_n}")
print(f"Max of N: {max_n}")
print(f"Standard Deviation of N: {std_n}")
print(f"Mean ± 3 SD: [{mean_3sd_lower}, {mean_3sd_upper}]")
print(f"Median Absolute Deviation of N: {mad_n}")
print(f"Median ± 3 MAD: [{median_3mad_lower}, {median_3mad_upper}]")
print(f"IQR of N: {iqr_n}")
print(f"IQR Range (Q1 - 1.5*IQR, Q3 + 1.5*IQR): [{iqr_lower}, {iqr_upper}]")

# Filter datasets based on each criterion
df_mean_3sd = df[(df['N'] >= mean_3sd_lower)]
df_median_3mad = df[(df['N'] >= median_3mad_lower) ]
df_iqr = df[(df['N'] >= iqr_lower)]

# keep SNPs with x% N e.g. 80%
# df = df[(df['N'] >= (max_n * X_percent) )]

# Function to create a Manhattan plot and save it
def plot_manhattan(df, output_path, title):
    fig, ax = plt.subplots(figsize=(12, 8.5))

    # Plot each chromosome with specific color
    for chrom in sorted(df['CHR'].unique()):
        df_chrom = df[df['CHR'] == chrom]
        chrom_color = chromosome_colors[chrom]

        # Plot non-significant points
        below_thresh = df_chrom[df_chrom['P'] >= threshold]
        ax.scatter(below_thresh['ind'], -np.log10(below_thresh['P']), c=chrom_color, marker='.', alpha=0.3)

        # Plot significant points based on EAF and Z score direction
        above_thresh_pos_dark = df_chrom[(df_chrom['P'] < threshold) & (df_chrom['Z'] > 0) & (df_chrom['EAF'] >= 0.01) & (df_chrom['EAF'] <= 0.99)]
        above_thresh_pos_light = df_chrom[(df_chrom['P'] < threshold) & (df_chrom['Z'] > 0) & ((df_chrom['EAF'] < 0.01) | (df_chrom['EAF'] > 0.99))]
        above_thresh_neg_dark = df_chrom[(df_chrom['P'] < threshold) & (df_chrom['Z'] < 0) & (df_chrom['EAF'] >= 0.01) & (df_chrom['EAF'] <= 0.99)]
        above_thresh_neg_light = df_chrom[(df_chrom['P'] < threshold) & (df_chrom['Z'] < 0) & ((df_chrom['EAF'] < 0.01) | (df_chrom['EAF'] > 0.99))]

        # Use 'tab:blue' for EAF between 0.01 and 0.99, and 'tab:cyan' for EAF outside this range
        ax.scatter(above_thresh_pos_dark['ind'], -np.log10(above_thresh_pos_dark['P']), c='tab:blue', marker='^', alpha=0.9)
        ax.scatter(above_thresh_pos_light['ind'], -np.log10(above_thresh_pos_light['P']), c='tab:cyan', marker='^', alpha=0.5)
        ax.scatter(above_thresh_neg_dark['ind'], -np.log10(above_thresh_neg_dark['P']), c='tab:blue', marker='v', alpha=0.9)
        ax.scatter(above_thresh_neg_light['ind'], -np.log10(above_thresh_neg_light['P']), c='tab:cyan', marker='v', alpha=0.5)

    # Threshold line and formatting
    ax.axhline(y=-np.log10(threshold), color='gray', linestyle='--', label=r'$P = 5 \times 10^{-8}$')
    ax.set_ylabel('-log10(P-value)')
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, fontsize=6, fontweight='bold')
    ax.set_title(title)
    ax.legend(loc='upper right')
        
    # Legend for significant SNPs
    # Legend for significant SNPs with unfilled shapes and color-coded text labels

    # Unfilled '^' marker for +ve Z with gray outline
    ax.scatter([], [], edgecolor='gray', facecolor='none', marker='^', label='+ve Z-score')

    # Unfilled 'v' marker for -ve Z with gray outline
    ax.scatter([], [], edgecolor='gray', facecolor='none', marker='v', label='-ve Z-score')

    # EAF categories 
    ax.scatter([], [], c='tab:blue', marker='s', alpha=0.9, label='EAF>=1%')
    ax.scatter([], [], c='tab:cyan', marker='s', alpha=0.5, label='0.1%<EAF<1%')
	       
    ax.legend(loc='upper right')

###################
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)

    plt.savefig(output_path, dpi=300)

    plt.close(fig)

###################################################
# Generate and save plots for each filtered dataset
plot_manhattan(df_mean_3sd, output_plot_mean_3sd, title_var)

plot_manhattan(df_median_3mad, output_plot_median_3mad, title_var)

plot_manhattan(df_iqr, output_plot_iqr, title_var)

