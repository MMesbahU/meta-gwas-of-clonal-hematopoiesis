# Load required packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)  # For arranging multiple plots into one figure
library(RColorBrewer)
library(ggpubr)
library(cowplot)
theme_set(theme_cowplot())
### Plot a: CHIP Prevalence by Age ###
load("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/rerun/fig1.data_TOPMed2019_UKB200k_UKB250k_MGBB53k.rda")

# Define the flag variables
flag_vars <- c("hasCHIP", "hasExpCHIP", "hasDNMT3A", 
               "hasTET2", "hasASXL1", "hasSF", "hasDDR")

# Filter data to include only ages between 20 and 89 and reshape to long format
long_df <- person_df_AncPred_wgs_N245388 %>% 
  filter(Age_biosample_collection >= 20 & Age_biosample_collection < 90) %>%
  pivot_longer(
    cols = all_of(flag_vars),
    names_to = "flag",
    values_to = "var"
  )

# Define custom colors for each flag
# custom_colors <- c(
# "hasCHIP"    = "black",   # Overall CHIP
# "hasExpCHIP" = "red",     # Expanded CHIP
# "hasDNMT3A"  = "blue",    # DNMT3A
# "hasTET2"    = "green",   # TET2
# "hasASXL1"   = "orange",  # ASXL1
# "hasSF"      = "purple",  # Splicing Factors
# "hasDDR"     = "brown"    # DNA damage
# )

custom_colors <- c(
  "hasCHIP"    = "black",     # Overall CHIP remains black
  "hasExpCHIP" = "#E41A1C",   # Red for Expanded CHIP
  "hasDNMT3A"  = "#377EB8",   # Blue for DNMT3A
  "hasTET2"    = "#4DAF4A",   # Green for TET2
  "hasASXL1"   = "#FF7F00",   # Bright Orange for ASXL1
  "hasSF"      = "#F781BF",   # Magenta/Pink for Splicing Factors
  "hasDDR"     = "#00BFC4"    # Turquoise for DNA damage
)


# Define custom legend labels for each flag
custom_labels <- c(
  "hasCHIP"    = "Overall CHIP",
  "hasExpCHIP" = "Expanded CHIP",
  "hasDNMT3A"  = "DNMT3A",
  "hasTET2"    = "TET2",
  "hasASXL1"   = "ASXL1",
  "hasSF"      = "Splicing Factors",
  "hasDDR"     = "DNA damage"
)




P1 <- ggplot(data = topmed2019.ukbb.mgbb, aes(x = Age, y = CHIP_status, group = Cohort)) +
  geom_smooth(aes(colour = Cohort), method = "glm", 
              method.args = list(family = "binomial"), se = TRUE) +
  xlab("Age") +
  ylab("Prevalence") +
  ggtitle("a") +
  scale_y_continuous(breaks = seq(0, 0.6, 0.1)) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  # scale_colour_manual(
  #   name = "CHIP Category",
  #   values = custom_colors,
  #   labels = custom_labels
  # ) +
  theme(
    legend.position = "right", 
    plot.title = element_text(size = 20, face = "bold")
  )

P1

#---------------- 2 -------------------
### Plot b: Proportion of Individuals by Number of CHIP Mutations ###

# Create table for mutation counts and compute proportions
count.aou250k <- as.data.frame(round(prop.table(table(table(noDup.aou.ch_var_withHemeCA$SampleID))) * 100, 1))
count.aou250k$Cohort <- "AoU"


p2 <- ggplot(data = chip_per_sample_2019, aes(x = reorder(Var1, -Freq), 
                                              y = Freq, fill=Cohort)) +
  xlab("Number of CHIP mutations") +
  ylab("Proportion of Individuals (%)") +
  geom_bar(stat = "identity",  width = 0.8, position = position_dodge()) +
  geom_text(aes(label = Freq), vjust = -0.5, color = "black",
            position = position_dodge(0.9), size = 3) +
  theme(legend.title = element_blank(),
        legend.position = "none", 
        plot.title = element_text(size = 20, face = "bold")) +
  ggtitle("b")

p2
#--------------------------------------

#----------------------- 3 -------------------
### Plot c: Proportion of Individuals by Top 20 CHIP Genes ###

# Create table for CHIP genes and compute proportions
CHIP_Gene_Table.AoU250k <- as.data.frame(round(prop.table(table(noDup.aou.ch_var_withHemeCA$Gene)) * 100, 1),
                                         stringsAsFactors = FALSE)
CHIP_Gene_Table.AoU250k$Cohort <- "AoU"
CHIP_Gene_Table.AoU250k <- CHIP_Gene_Table.AoU250k[order(CHIP_Gene_Table.AoU250k$Freq, decreasing = TRUE),]

CHIP_Gene_Table.2019 <- CHIP_Gene_Table.2019[order(CHIP_Gene_Table.2019$Freq, decreasing = TRUE), ]

# --- Define the gene color mapping (same as in Plot d) ---
gene_colors <- c(
  "DNMT3A" = "#377EB8",  # Blue for DNMT3A
  "TET2"   = "#4DAF4A",  # Green for TET2
  "ASXL1"  = "#FF7F00",  # Bright Orange for ASXL1
  # Splicing Factors get the same vivid magenta/pink:
  "SF3B1"  = "#F781BF",
  "SRSF2"  = "#F781BF",
  "U2AF1"  = "#F781BF",
  "ZRSR2"  = "#F781BF",
  # DNA damage genes get turquoise:
  "PPM1D"  = "#00BFC4",
  "TP53"   = "#00BFC4"
)

# --- Define the color mapping for Plot c ---
# Extract the top 20 genes from the CHIP gene table (assumed to be in column 'Var1')
top20_genes <- head(unique(CHIP_Gene_Table.2019$Var1), 20)
unique_genes_c <- sort(unique(top20_genes))
mapped_genes_c <- intersect(unique_genes_c, names(gene_colors))
unmapped_genes_c <- setdiff(unique_genes_c, names(gene_colors))

if (length(unmapped_genes_c) > 0) {
  n_needed <- length(unmapped_genes_c)
  # Use up to 8 colors from the "Set2" palette; if more are needed, expand the palette
  if(n_needed <= 8){
    fallback_colors <- brewer.pal(n = n_needed, name = "Set2")
  } else {
    fallback_colors <- colorRampPalette(brewer.pal(8, "Set2"))(n_needed)
  }
  fallback_mapping <- setNames(fallback_colors, unmapped_genes_c)
} else {
  fallback_mapping <- c()
}

# Combine the pre-defined gene colors with fallback colors
gene_colors_plot_c <- c(gene_colors, fallback_mapping)
# Order the mapping to match the order of genes in the top20 list
gene_colors_plot_c <- gene_colors_plot_c[unique_genes_c]

# --- Plot c: Proportion of Individuals by Top 20 CHIP Genes with Gene Colors ---
p3 <- CHIP_Gene_Table.2019 %>% 
  filter(Var1 %in% top20_genes) %>% 
  ggplot(aes(x = reorder(Var1, -Freq), y = Freq, fill = Cohort)) +
  xlab("") +
  ylab("Proportion of Individuals (%)") +
  geom_bar(stat = "identity", width = 0.8, position = position_dodge()) +
  geom_text(aes(label = Freq), vjust = -0.5, color = "black",
            position = position_dodge(0.9), size = 3) +
  scale_fill_manual(values = gene_colors_plot_c) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = "none", 
        plot.title = element_text(size = 20, face = "bold")) +
  ggtitle("c")

# Display Plot c
p3

#-------------- 4
### Plot d: Variant Allele Fraction by Top 20 Genes ###

AoU_top20_genes <- noDup.aou.ch_var_withHemeCA %>% 
  filter(Gene %in% head(CHIP_Gene_Table.AoU250k, n = 20)$Var1)

# Ensure the Gene factor is ordered by the table order
AoU_top20_genes$Gene <- factor(AoU_top20_genes$Gene, 
                               levels = head(CHIP_Gene_Table.AoU250k, n = 20)$Var1)


# Define the color mapping for genes consistent with Plot a
gene_colors <- c(
  "DNMT3A" = "#377EB8",  # Blue for DNMT3A
  "TET2"   = "#4DAF4A",  # Green for TET2
  "ASXL1"  = "#FF7F00",  # Bright Orange for ASXL1
  # Splicing Factors (all get the same vivid magenta/pink)
  "SF3B1"  = "#F781BF",
  "SRSF2"  = "#F781BF",
  "U2AF1"  = "#F781BF",
  "ZRSR2"  = "#F781BF",
  # DNA damage genes (all get the same turquoise)
  "PPM1D"  = "#00BFC4",
  "TP53"   = "#00BFC4"
)

# Get the unique genes from your top 20 gene dataset
unique_genes <- sort(unique(AoU_top20_genes$Gene))
# Determine which genes are already in our mapping
mapped_genes <- intersect(unique_genes, names(gene_colors))
# Identify genes not yet assigned a color
unmapped_genes <- setdiff(unique_genes, names(gene_colors))

# For unmapped genes, assign distinct colors from the "Set2" palette
if (length(unmapped_genes) > 0) {
  # RColorBrewer's Set2 has up to 8 colors. If more are needed, expand the palette.
  n_needed <- length(unmapped_genes)
  if(n_needed <= 8){
    fallback_colors <- brewer.pal(n = n_needed, name = "Set2")
  } else {
    fallback_colors <- colorRampPalette(brewer.pal(8, "Set2"))(n_needed)
  }
  fallback_mapping <- setNames(fallback_colors, unmapped_genes)
} else {
  fallback_mapping <- c()
}

# Combine the defined colors with the fallback mapping
gene_colors_plot <- c(gene_colors, fallback_mapping)
# Optionally, reorder to match the order of genes in your dataset
gene_colors_plot <- gene_colors_plot[unique_genes]

# Create Plot d with the updated color scheme
p4 <- ggplot(data = AoU_top20_genes, aes(x = Gene, y = VAF, fill = Gene)) + 
  xlab("") +
  ylab("Variant Allele Fraction") + 
  geom_boxplot() + 
  scale_y_log10(breaks = c(0.02, 0.1, 0.2, 0.3, 0.5, 1)) +
  scale_fill_manual(values = gene_colors_plot) +
  stat_summary(fun = "median", geom = "point", color = "white", 
               position = position_dodge(0.9)) +
  ggtitle("d") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.title = element_blank(), 
        legend.position = "none",
        plot.title = element_text(size = 20, face = "bold"))

# Display Plot d
p4

### Combine the Four Plots into a 2×2 Multipanel Figure ###

# Arrange the plots: first row (P1 and p2) and second row (p3 and p4)
multi_plot <- (P1 | p2) / (p3 | p4)

# Display the combined multipanel plot
multi_plot

### Save the Combined Plot as PNG and PDF ###

# Save as PNG
# ggsave("~/FigS1.AoU.multipanel_plot.png", multi_plot, width = 16, height = 10, dpi = 300)

# Save as PDF
# ggsave("~/FigS1.AoU.multipanel_plot.pdf", multi_plot, width = 16, height = 10)