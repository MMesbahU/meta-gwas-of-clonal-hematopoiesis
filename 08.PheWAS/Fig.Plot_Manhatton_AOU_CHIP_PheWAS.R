#-----------------------
# Manhattan Plot of CHIP PheWAS in AoU
#-----------------------
library(data.table)
library(dplyr)
setwd("~/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/PheWAS/AoU2025/")
#-----------------------
phewas <- fread("AoU_phewas.final_results.MultiANC.min40case_counts.2025.csv", header = TRUE)

# 
pheinfo <- fread("~/Documents/Project/all0fus/Phecodes/data/pheinfo.csv")
pheinfo$Phe_phecode <- paste0("Phe_",pheinfo$phecode)
head(pheinfo$Phe_phecode)

#------------------------
phewas_annot <- merge(phewas[,c(2:10)], pheinfo[,c(2:7)], 
                     by.x="Outcome", 
                     by.y = "Phe_phecode")

write.csv(phewas_annot, "aou_chip.phewas_annot.csv", row.names = F)
#write.csv(phewas_annot[!(phewas_annot$Exposure=="hasIDH1IDH2"),], "aou_chip.phewas_annot.2025.csv", row.names = F)

phewas_annot_noidh <- phewas_annot[!(phewas_annot$Exposure=="hasIDH1IDH2"),]

phewas_annot_noidh$EXPOSURE <- case_when(
  phewas_annot_noidh$Exposure == "hasCH" ~ "Overall CHIP",
  phewas_annot_noidh$Exposure == "hasCHvaf10" ~ "Expanded CHIP",
  phewas_annot_noidh$Exposure == "hasDDR" ~ "DNA damage response genes",
  phewas_annot_noidh$Exposure == "hasSF" ~ "Splicing factors",
  phewas_annot_noidh$Exposure == "hasDNMT3A" ~ "DNMT3A CHIP",
  phewas_annot_noidh$Exposure == "hasTET2" ~ "TET2 CHIP",
  phewas_annot_noidh$Exposure == "hasASXL1" ~ "ASXL1 CHIP",
  phewas_annot_noidh$Exposure == "hasPPM1D" ~ "PPM1D CHIP",
  phewas_annot_noidh$Exposure == "hasTP53" ~ "TP53 CHIP",
  phewas_annot_noidh$Exposure == "hasSF3B1" ~ "SF3B1 CHIP",
  phewas_annot_noidh$Exposure == "hasJAK2" ~ "JAK2 CHIP",
  TRUE ~ "NA" # Default value
)

table(phewas_annot_noidh$EXPOSURE, exclude = NULL)

#-----------------------------------
# Manhattan Plot
#-----------------------------------
library(ggplot2)
library(dplyr)
library(ggrepel)

Exposures <- unique(phewas_annot_noidh$EXPOSURE)
#-----------------------------------
# Main Plot
library(ggplot2)
library(dplyr)
library(ggrepel)
n_outcome <- length(unique(phewas_annot_noidh$description))
# 1774
n_exposure <- length(unique(phewas_annot_noidh$Exposure))
# 11
p_threshold <- as.numeric(formatC(x = 0.05/n_outcome, digits = 1, format = "E"))
p_bonferoni <- as.numeric(formatC(x = 0.05/(n_outcome*n_exposure), digits = 1, format = "E"))

#### Adjusted P
phewas_annot_noidh$FDR <- p.adjust(phewas_annot_noidh$P,method = "fdr")

res_1$FDR <- p.adjust(res_1$P,method = "fdr")
res_1 <- subset(phewas_annot_noidh, phewas_annot_noidh$Exposure=="hasCH")


#---
Exposures <- unique(phewas_annot_noidh$EXPOSURE)

for(i in Exposures){
  # 1. Subset data for the current Exposure
  df <- phewas_annot_noidh[phewas_annot_noidh$EXPOSURE == i, c(5,6,10:14)]
  
  # 2. Summarize group info
  #    We'll add a 'gap' in groupStart so that each group has spacing on the x-axis.
  df_group <- df %>%
    group_by(group, groupnum) %>%
    summarise(n_group = n(), .groups = "drop") %>%
    arrange(groupnum) %>%
    # Create an incrementing gap for each group so there's extra space between groups
    mutate(
      # row_number() gives the group index in ascending order
      group_index = row_number(),
      # Add some constant spacing (e.g., 2 units) * group_index
      # Use cumsum(n_group) to track the end of the last group + spacing
      groupStart = cumsum(lag(n_group, default = 0)) + 2 * (group_index - 1),
      mid = groupStart + n_group / 2
    )
  
  # 3. Merge group info and compute x, shape_var, etc.
  df_plot <- df %>%
    left_join(df_group %>% select(group, groupStart), by = "group") %>%
    group_by(group) %>%
    mutate(
      x = groupStart + row_number(),
      shape_var = ifelse(Z > 0, 24, 25),
      negLogP = -log10(P)
    ) %>%
    ungroup()
  
  # 4. Plot
  sig_line <- -log10(p_threshold)
  sig_line_bonferoni <- -log10(p_bonferoni)
  p <- ggplot(df_plot, aes(x = x, y = negLogP)) +
    geom_point(aes(shape = shape_var, color = color)) +
    scale_shape_identity() +
    scale_color_identity() +
    geom_hline(yintercept = sig_line, 
               linetype = "dashed", color = "blue") +
    geom_hline(yintercept = sig_line_bonferoni, 
               linetype = "dashed", color = "red") +
    geom_text_repel(
      data = subset(df_plot, P < p_threshold),
      aes(label = description),
      box.padding = 0.4,
      max.overlaps = Inf,
      force = 1,
      segment.alpha = 0.5,
      size = 2.5
    ) +
    scale_x_continuous(
      breaks = df_group$mid,
      labels = df_group$group
    ) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    ) +
    labs(
      x = "Phecode Groups",
      y = expression(-log[10](P)),
      title = paste("PheWAS Manhattan Plot for", i)
    )
  
  # 5. Save the plot as a high-resolution PDF
  ggsave(
    filename = paste0("AoU.PheWAS_Manhattan_Plot_", i, ".pdf"),
    plot = p,
    device = "pdf",
    width = 12,       # adjust as needed
    height = 8,       # adjust as needed
    dpi = 600         # high resolution
  )
}



write.csv(phewas_annot_noidh, "AoU.phewas_annot_noidh.chip_gwas.2025.csv", row.names = FALSE)
#-------------------------------


