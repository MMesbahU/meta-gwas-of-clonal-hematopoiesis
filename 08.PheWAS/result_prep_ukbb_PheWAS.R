rm(list=ls())
#------------------------------------
# UKBB Updated CHIP PheWAS by Maryam
#------------------------------------
library(data.table)
library(dplyr)
setwd("~/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/PheWAS/UKBB_MGBB_MVP_BioVU/ukb_phewas")
#-----------------------
# 
pheinfo <- fread("~/Documents/Project/all0fus/Phecodes/data/pheinfo.csv")
# pheinfo$Phe_phecode <- paste0("Phe_",pheinfo$phecode)
# head(pheinfo$Phe_phecode)
pheinfo$desc_string <- gsub("[^[:alnum:]]+", "_", pheinfo$description)

#----------
# phewas result
phewas_circ <- fread("NewCHIP_Phewas_CoxPH.PheCodePhenos.CirculatoryRespiratory.txt", header = TRUE)
table(phewas_circ$y %in% pheinfo$desc_string)

phewas_derm <- fread("NewCHIP_Phewas_CoxPH.PheCodePhenos.DermDigestGU.txt", header = TRUE)
table(phewas_derm$y %in% pheinfo$desc_string)

phewas_heam <- fread("NewCHIP_Phewas_CoxPH.PheCodePhenos.HematologicNeoplasmInfectious.txt", header = TRUE)
table(phewas_heam$y %in% pheinfo$desc_string)

phewas_injur <- fread("NewCHIP_Phewas_CoxPH.PheCodePhenos.InjuriesPoisoningMSK.txt", header = TRUE)
table(phewas_injur$y %in% pheinfo$desc_string)

phewas_mental <- fread("NewCHIP_Phewas_CoxPH.PheCodePhenos.MentalNeuroSensorySymptoms.txt", header = TRUE)
table(phewas_mental$y %in% pheinfo$desc_string)

phewas_preg <- fread("NewCHIP_Phewas_CoxPH.PheCodePhenos.PregnancyCongenitalEndocrine.txt", header = TRUE)
table(phewas_preg$y %in% pheinfo$desc_string)

res.ukbb_phewas <- as.data.frame(rbind(phewas_circ, 
                                       phewas_derm,
                                       phewas_heam, 
                                       phewas_injur,
                                       phewas_mental,
                                       phewas_preg))
rm(phewas_circ, 
   phewas_derm,
   phewas_heam, 
   phewas_injur,
   phewas_mental,
   phewas_preg)

names(res.ukbb_phewas)
names(res.ukbb_phewas) <- c("Beta", "HR", "SE", "Z", "P", "Outcome", "Exposure",names(res.ukbb_phewas)[8:11])
names(res.ukbb_phewas)
## annotate
annot.res.ukbb_phewas <- merge(res.ukbb_phewas, 
                               pheinfo[,c(2:7)], 
                               by.x="Outcome", 
                               by.y="desc_string")

names(annot.res.ukbb_phewas)

annot.res.ukbb_phewas$EXPOSURE <- case_when(
  annot.res.ukbb_phewas$Exposure == "NEWCALLS_hasCH" ~ "Overall CHIP",
  annot.res.ukbb_phewas$Exposure == "NEWCALLS_hasCHvaf10" ~ "Expanded CHIP",
  annot.res.ukbb_phewas$Exposure == "NEWCALLS_hasDDR" ~ "DNA damage response genes",
  annot.res.ukbb_phewas$Exposure == "NEWCALLS_hasSF" ~ "Splicing factors",
  annot.res.ukbb_phewas$Exposure == "NEWCALLS_hasDNMT3A" ~ "DNMT3A CHIP",
  annot.res.ukbb_phewas$Exposure == "NEWCALLS_hasTET2" ~ "TET2 CHIP",
  annot.res.ukbb_phewas$Exposure == "NEWCALLS_hasASXL1" ~ "ASXL1 CHIP",
  annot.res.ukbb_phewas$Exposure == "NEWCALLS_hasDTA" ~ "DNMT3A/TET2/ASXL1 CHIP",
  annot.res.ukbb_phewas$Exposure == "NEWCALLS_hasCHvaf05" ~ "CHIP (VAF>=5%)",
    TRUE ~ "NA" # Default value
  )

table(annot.res.ukbb_phewas$Exposure, exclude = NULL)

#write.csv(annot.res.ukbb_phewas, "NewCHIP_Phewas_CoxPH.PheCodePhenos.All_result.csv", row.names = F)

# filter data
qc_annot.res.ukbb_phewas <- subset(annot.res.ukbb_phewas, !(annot.res.ukbb_phewas$Exposure %in% c("NEWCALLS_hasDTA", "NEWCALLS_hasCHvaf05")) )

qc_annot.res.ukbb_phewas$FDR <- p.adjust(qc_annot.res.ukbb_phewas$P, method = "fdr")

# write.csv(qc_annot.res.ukbb_phewas, "final_result.NewCHIP_Phewas_CoxPH.PheCodePhenos.fdr.csv", row.names = F)

## 
p_threshold <- 0.05/length(unique(qc_annot.res.ukbb_phewas$Outcome))
# 
# 7 CHIP * 1451 Phecodes = 10157 Tests
p_bonferoni <- 0.05/(length(unique(qc_annot.res.ukbb_phewas$Outcome)) * length(unique(qc_annot.res.ukbb_phewas$EXPOSURE)))

#---------------------
### Plot
Exposures <- unique(qc_annot.res.ukbb_phewas$EXPOSURE)

for(i in Exposures){
  # 1. Subset data for the current Exposure
  df <- qc_annot.res.ukbb_phewas[qc_annot.res.ukbb_phewas$EXPOSURE == i, ]
  
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
    filename = paste0("UKBB.PheWAS_Manhattan_Plot_", i, ".pdf"),
    plot = p,
    device = "pdf",
    width = 12,       # adjust as needed
    height = 8,       # adjust as needed
    dpi = 600         # high resolution
  )
}

