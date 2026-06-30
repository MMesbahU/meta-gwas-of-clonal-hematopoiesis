
# Install/load if not already
# install.packages("phewas")
library(phewas)
library(data.table)
library(dplyr)
library(ggrepel)
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

# write.csv(phewas_annot, "aou_chip.phewas_annot.csv", row.names = F)
#write.csv(phewas_annot[!(phewas_annot$Exposure=="hasIDH1IDH2"),], "aou_chip.phewas_annot.2025.csv", row.names = F)

# phewas_annot_noidh <- phewas_annot[!(phewas_annot$Exposure=="hasIDH1IDH2"),]
# phewas_annot <- phewas_annot_noidh
phewas_annot$EXPOSURE <- case_when(
  phewas_annot$Exposure == "hasCH" ~ "Overall CHIP",
  phewas_annot$Exposure == "hasCHvaf10" ~ "Expanded CHIP",
  phewas_annot$Exposure == "hasDDR" ~ "DNA damage response genes",
  phewas_annot$Exposure == "hasSF" ~ "Splicing factors",
  phewas_annot$Exposure == "hasDNMT3A" ~ "DNMT3A CHIP",
  phewas_annot$Exposure == "hasTET2" ~ "TET2 CHIP",
  phewas_annot$Exposure == "hasASXL1" ~ "ASXL1 CHIP",
  phewas_annot$Exposure == "hasPPM1D" ~ "PPM1D CHIP",
  phewas_annot$Exposure == "hasTP53" ~ "TP53 CHIP",
  phewas_annot$Exposure == "hasSF3B1" ~ "SF3B1 CHIP",
  phewas_annot$Exposure == "hasJAK2" ~ "JAK2 CHIP",
  phewas_annot$Exposure == "hasIDH1IDH2" ~ "IDH1 IDH2 CHIP",
  TRUE ~ "NA" # Default value
)

table(phewas_annot$EXPOSURE, exclude = NULL)

#-----------------------------------
# Manhattan Plot
#-----------------------------------
n_outcome <- length(unique(phewas_annot$description))
# 1774
n_exposure <- length(unique(phewas_annot$Exposure))
# 12
p_threshold <- as.numeric(formatC(x = 0.05/n_outcome, digits = 1, format = "E"))
p_bonferoni <- as.numeric(formatC(x = 0.05/(n_outcome*n_exposure), digits = 1, format = "E"))

#### Adjusted P
phewas_annot$FDR <- p.adjust(phewas_annot$P,method = "fdr")
max(phewas_annot$P[phewas_annot$FDR<=0.05])

#---
Exposures <- unique(phewas_annot$EXPOSURE)


# Make sure you've sourced/defined your phenotypeManhattan() function
for(i in Exposures){
# Subset for a single exposure
df_sub <- phewas_annot[phewas_annot$EXPOSURE == i, ]

# Create the data frame with the correct columns
d <- df_sub %>%
  dplyr::mutate(
    phenotype = as.character(Outcome),  # e.g., "Phe_10", "Phe_110", etc.
    p         = P,
    OR        = exp(Beta)              # if Beta is log(OR) and we want OR-based visualization
  ) %>%
  dplyr::select(phenotype, p, OR, Beta, SE, Z, everything()) 
# ^ includes other columns if you like, but minimally phenotype & p

# Now pass this to phenotypeManhattan()
my_plot <- phenotypeManhattan(
  d,
  suggestive.line = NA,# p_threshold,       # typical "suggestive" threshold
  significant.line=p_threshold, #p_bonferoni,
#  OR.size = TRUE,              # use bubble size by OR
  OR.direction = TRUE,         # use shape for direction (OR>1 vs. OR<1)
  annotate.level = p_threshold,        # label points above -log10(1e-4)
  annotate.phenotype.description=TRUE, 
  annotate.size=3,
  #annotate.angle=0,
 # size.x.labels=12,	
  #size.y.labels=12,
  title = paste("PheWAS Manhattan Plot for", i),
# title = paste("PheWAS Manhattan Plot for", i, "(P<=",p_threshold,")"),
  sort.by.category.value=T, 
# point.size=3, 
y.axis.interval=10, 
# switch.axis=F, 
# base.labels=T
)
# Print or save your plot
print(my_plot)

ggsave(
  filename = paste0("../../../NG_v2/Figs/phenotypeManhattan.AoU_glm.PheWAS_Manhattan_Plot_", i, ".pdf"),
  plot   = my_plot,
  width  = 10,
  height = 6,
  dpi    = 300
)
}
