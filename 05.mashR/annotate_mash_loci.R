## Annotate Loci ####
# +-500kb

### Annovar Annotations ####
library(data.table)
library(dplyr)

mash_sorted <- fread("/Users/muddin/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/mash/sorted.chip_dnmt_tet_ltl_mCA_mpn.lfsr05.oct2024.csv", header = TRUE)
str(mash_sorted)

var_annotations <- fread("/Users/muddin/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/mash/annot_hg37.varID.hg19_multianno.txt.gz", header = TRUE)
str(var_annotations)

mas.annot <- merge(mash_sorted, var_annotations, by.x="varID", by.y="Otherinfo6")
head(mas.annot)
# Sort by chr and pos
mas.annot_sorted <- mas.annot %>% arrange(Otherinfo4, Otherinfo5)

# fwrite(mas.annot, "/Users/muddin/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/mash/annoated.chip_dnmt_tet_ltl_mCA_mpn.lfsr05.oct2024.csv", 
#        row.names = F, col.names = T, sep=",")
# 
# fwrite(mas.annot_sorted, "/Users/muddin/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/mash/annoated.sorted.chip_dnmt_tet_ltl_mCA_mpn.lfsr05.oct2024.csv", 
#        row.names = F, col.names = T, sep=",")
# fwrite(mas.annot_sorted, "/Users/muddin/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/mash/annoated.sorted.chip_dnmt_tet_ltl_mCA_mpn.lfsr05.oct2024.tsv", 
#        row.names = F, col.names = T, sep="\t", quote = T)

## #### lfsr<=0.05 | GWAS P <= 5e-8
head(sort(table(mas.annot_sorted$Gene.refGene[mas.annot_sorted$`CHIP LFSR`<0.05 & mas.annot_sorted$`CHIP GWAS P`>5e-8]), decreasing = T))

# new_chip <- mas.annot_sorted[mas.annot_sorted$`CHIP LFSR`<0.05 & mas.annot_sorted$`CHIP GWAS P`>5e-8,]
# new_dnmt <- mas.annot_sorted[mas.annot_sorted$`DNMT3A LFSR`<0.05 & mas.annot_sorted$`DNMT3A GWAS P`>5e-8,]
# new_tet <- mas.annot_sorted[mas.annot_sorted$`TET2 LFSR`<0.05 & mas.annot_sorted$`TET2 GWAS P`>5e-8,]
# new_ltl <- mas.annot_sorted[mas.annot_sorted$`LTL LFSR`<0.05 & mas.annot_sorted$`LTL GWAS P`>5e-8,]
# new_mpn <- mas.annot_sorted[mas.annot_sorted$`MPN LFSR`<0.05 & mas.annot_sorted$`MPN GWAS P`>5e-8,]
# new_mca <- mas.annot_sorted[mas.annot_sorted$`mCA LFSR`<0.05 & mas.annot_sorted$`mCA GWAS P`>5e-8,]

new_chip <- mas.annot_sorted[mas.annot_sorted$`CHIP LFSR`<0.05,]
# add annotations 
# Convert to data.table
dt <- as.data.table(new_chip)

# Sort by chromosome and position
setorder(dt, Chr, Start)

# Initialize locus count
dt[, Locus_count := NA_integer_]
locus_id <- 1
i <- 1

while (i <= nrow(dt)) {
  # Define window for current locus
  chr_i <- dt[i, Chr]
  pos_i <- dt[i, Start]
  window_start <- pos_i - 500000
  window_end <- pos_i + 500000
  
  # Find all variants within the ±500kb window on the same chromosome
  idx <- which(dt$Chr == chr_i & dt$Start >= window_start & dt$Start <= window_end)
  
  # Assign the same locus ID
  dt[idx, Locus_count := locus_id]
  
  # Move to the next unassigned row beyond this locus
  i <- max(idx) + 1
  locus_id <- locus_id + 1
}

# Now annotate New_locus
# For each locus, check if any variant has GWAS P ≤ 5e-8
dt[, New_locus := if (any(`CHIP GWAS P` <= 5e-8, na.rm = TRUE)) "No" else "Yes", by = Locus_count]

new_chip <- dt[,c(111,26,112,1:110)]
names(new_chip) <- c("No.", "Locus","Novel",names(dt[,c(1:110)]))
fwrite(new_chip, "/Users/muddin/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/mash/annoated.sorted.chip.lfsr05.oct2024.csv", row.names = F, col.names = T, sep=",")

# dnmt 
new_dnmt <- mas.annot_sorted[mas.annot_sorted$`DNMT3A LFSR`<0.05,]
dt <- as.data.table(new_dnmt)

# Sort by chromosome and position
setorder(dt, Chr, Start)

# Initialize locus count
dt[, Locus_count := NA_integer_]
locus_id <- 1
i <- 1

while (i <= nrow(dt)) {
  # Define window for current locus
  chr_i <- dt[i, Chr]
  pos_i <- dt[i, Start]
  window_start <- pos_i - 500000
  window_end <- pos_i + 500000
  
  # Find all variants within the ±500kb window on the same chromosome
  idx <- which(dt$Chr == chr_i & dt$Start >= window_start & dt$Start <= window_end)
  
  # Assign the same locus ID
  dt[idx, Locus_count := locus_id]
  
  # Move to the next unassigned row beyond this locus
  i <- max(idx) + 1
  locus_id <- locus_id + 1
}

# Now annotate New_locus
# For each locus, check if any variant has GWAS P ≤ 5e-8
dt[, New_locus := if (any(`DNMT3A GWAS P` <= 5e-8, na.rm = TRUE)) "No" else "Yes", by = Locus_count]

new_dnmt <- dt[,c(111,26,112,1:110)]
names(new_dnmt) <- c("No.", "Locus","Novel",names(dt[,c(1:110)]))
fwrite(new_dnmt, "/Users/muddin/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/mash/annoated.sorted.dnmt3a.lfsr05.oct2024.csv", row.names = F, col.names = T, sep=",")

# TET
new_tet <- mas.annot_sorted[mas.annot_sorted$`TET2 LFSR`<0.05,]
dt <- as.data.table(new_tet)

# Sort by chromosome and position
setorder(dt, Chr, Start)

# Initialize locus count
dt[, Locus_count := NA_integer_]
locus_id <- 1
i <- 1

while (i <= nrow(dt)) {
  # Define window for current locus
  chr_i <- dt[i, Chr]
  pos_i <- dt[i, Start]
  window_start <- pos_i - 500000
  window_end <- pos_i + 500000
  
  # Find all variants within the ±500kb window on the same chromosome
  idx <- which(dt$Chr == chr_i & dt$Start >= window_start & dt$Start <= window_end)
  
  # Assign the same locus ID
  dt[idx, Locus_count := locus_id]
  
  # Move to the next unassigned row beyond this locus
  i <- max(idx) + 1
  locus_id <- locus_id + 1
}

# Now annotate New_locus
# For each locus, check if any variant has GWAS P ≤ 5e-8
dt[, New_locus := if (any(`TET2 GWAS P` <= 5e-8, na.rm = TRUE)) "No" else "Yes", by = Locus_count]

new_tet <- dt[,c(111,26,112,1:110)]
names(new_tet) <- c("No.", "Locus","Novel",names(dt[,c(1:110)]))

fwrite(new_tet, "/Users/muddin/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/mash/annoated.sorted.tet2.lfsr05.oct2024.csv", row.names = F, col.names = T, sep=",")

# LTL ####
new_ltl <- mas.annot_sorted[mas.annot_sorted$`LTL LFSR`<0.05,]
dt <- as.data.table(new_ltl)

# Sort by chromosome and position
setorder(dt, Chr, Start)

# Initialize locus count
dt[, Locus_count := NA_integer_]
locus_id <- 1
i <- 1

while (i <= nrow(dt)) {
  # Define window for current locus
  chr_i <- dt[i, Chr]
  pos_i <- dt[i, Start]
  window_start <- pos_i - 500000
  window_end <- pos_i + 500000
  
  # Find all variants within the ±500kb window on the same chromosome
  idx <- which(dt$Chr == chr_i & dt$Start >= window_start & dt$Start <= window_end)
  
  # Assign the same locus ID
  dt[idx, Locus_count := locus_id]
  
  # Move to the next unassigned row beyond this locus
  i <- max(idx) + 1
  locus_id <- locus_id + 1
}

# Now annotate New_locus
# For each locus, check if any variant has GWAS P ≤ 5e-8
dt[, New_locus := if (any(`LTL GWAS P` <= 5e-8, na.rm = TRUE)) "No" else "Yes", by = Locus_count]

new_ltl <- dt[,c(111,26,112,1:110)]
names(new_ltl) <- c("No.", "Locus","Novel",names(dt[,c(1:110)]))

fwrite(new_ltl, "/Users/muddin/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/mash/annoated.sorted.ltl.lfsr05.oct2024.csv", row.names = F, col.names = T, sep=",")

# MPN ####
new_mpn <- mas.annot_sorted[mas.annot_sorted$`MPN LFSR`<0.05,]
dt <- as.data.table(new_mpn)

# Sort by chromosome and position
setorder(dt, Chr, Start)

# Initialize locus count
dt[, Locus_count := NA_integer_]
locus_id <- 1
i <- 1

while (i <= nrow(dt)) {
  # Define window for current locus
  chr_i <- dt[i, Chr]
  pos_i <- dt[i, Start]
  window_start <- pos_i - 500000
  window_end <- pos_i + 500000
  
  # Find all variants within the ±500kb window on the same chromosome
  idx <- which(dt$Chr == chr_i & dt$Start >= window_start & dt$Start <= window_end)
  
  # Assign the same locus ID
  dt[idx, Locus_count := locus_id]
  
  # Move to the next unassigned row beyond this locus
  i <- max(idx) + 1
  locus_id <- locus_id + 1
}

# Now annotate New_locus
# For each locus, check if any variant has GWAS P ≤ 5e-8
dt[, New_locus := if (any(`MPN GWAS P` <= 5e-8, na.rm = TRUE)) "No" else "Yes", by = Locus_count]

new_mpn <- dt[,c(111,26,112,1:110)]
names(new_mpn) <- c("No.", "Locus","Novel",names(dt[,c(1:110)]))

fwrite(new_mpn, "/Users/muddin/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/mash/annoated.sorted.mpn.lfsr05.oct2024.csv", row.names = F, col.names = T, sep=",")

# mCA ####
new_mca <- mas.annot_sorted[mas.annot_sorted$`mCA LFSR`<0.05,]
dt <- as.data.table(new_mca)

# Sort by chromosome and position
setorder(dt, Chr, Start)

# Initialize locus count
dt[, Locus_count := NA_integer_]
locus_id <- 1
i <- 1

while (i <= nrow(dt)) {
  # Define window for current locus
  chr_i <- dt[i, Chr]
  pos_i <- dt[i, Start]
  window_start <- pos_i - 500000
  window_end <- pos_i + 500000
  
  # Find all variants within the ±500kb window on the same chromosome
  idx <- which(dt$Chr == chr_i & dt$Start >= window_start & dt$Start <= window_end)
  
  # Assign the same locus ID
  dt[idx, Locus_count := locus_id]
  
  # Move to the next unassigned row beyond this locus
  i <- max(idx) + 1
  locus_id <- locus_id + 1
}

# Now annotate New_locus
# For each locus, check if any variant has GWAS P ≤ 5e-8
dt[, New_locus := if (any(`mCA GWAS P` <= 5e-8, na.rm = TRUE)) "No" else "Yes", by = Locus_count]

new_mca <- dt[,c(111,26,112,1:110)]
names(new_mca) <- c("No.", "Locus","Novel",names(dt[,c(1:110)]))
fwrite(new_mca, "/Users/muddin/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/mash/annoated.sorted.mca.lfsr05.oct2024.csv", row.names = F, col.names = T, sep=",")

###
# fwrite(mas.annot_sorted[mas.annot_sorted$`CHIP LFSR`<0.05,c(1:4,20:110)], "/Users/muddin/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/mash/annoated.sorted.chip_all.lfsr05.oct2024.csv", row.names = F, col.names = T, sep=",")
# fwrite(mas.annot_sorted[mas.annot_sorted$`DNMT3A LFSR`<0.05,c(1,5:7,20:110)], "/Users/muddin/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/mash/annoated.sorted.dnmt3a_all.lfsr05.oct2024.csv", row.names = F, col.names = T, sep=",")
# fwrite(mas.annot_sorted[mas.annot_sorted$`TET2 LFSR`<0.05,c(1,8:10,20:110)], "/Users/muddin/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/mash/annoated.sorted.tet2_all.lfsr05.oct2024.csv", row.names = F, col.names = T, sep=",")
# fwrite(mas.annot_sorted[mas.annot_sorted$`LTL LFSR`<0.05,c(1,11:13,20:110)], "/Users/muddin/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/mash/annoated.sorted.ltl_all.lfsr05.oct2024.csv", row.names = F, col.names = T, sep=",")
# fwrite(mas.annot_sorted[mas.annot_sorted$`mCA LFSR`<0.05,c(1,14:16,20:110)], "/Users/muddin/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/mash/annoated.sorted.mca_all.lfsr05.oct2024.csv", row.names = F, col.names = T, sep=",")
# fwrite(mas.annot_sorted[mas.annot_sorted$`MPN LFSR`<0.05,c(1,17:19,20:110)], "/Users/muddin/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/mash/annoated.sorted.mpn.lfsr05.oct2024.csv", row.names = F, col.names = T, sep=",")


## POPS Annotations ####
library(data.table)
library(dplyr)
library(biomaRt)
require(openxlsx)

## 

# Select the human genes dataset
ensembl <- useEnsembl(biomart = "ensembl", 
                      dataset = "hsapiens_gene_ensembl", 
                      # mirror = "useast", # version or GRCh arguments cannot be used together with the mirror argument
                      GRCh=37)
##
setwd("~/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/pops/") 
ch <- fread("CH_pops_gene.preds.txt", header = T, na.strings = "NA")
dnmt <- fread("DNMT3A_pops_gene.preds.txt", header = T, na.strings = "NA")
tet <- fread("TET2_pops_gene.preds.txt", header = T, na.strings = "NA")

######
# Get the gene intervals
gene_intervals <- getBM(
  attributes = c("ensembl_gene_id", 
                 "chromosome_name", 
                 "start_position", 
                 "end_position"),
  filters = "ensembl_gene_id",
  values = ch$ENSGID,
  mart = ensembl
)
#########
table(gene_intervals$ensembl_gene_id %in% ch$ENSGID)

ch <- merge(ch, gene_intervals, 
            by.x="ENSGID",
            by.y="ensembl_gene_id")

ch$region <- paste0(ch$chromosome_name,":",
                    ch$start_position,"-",
                    ch$end_position)

ch <- ch %>% arrange(chromosome_name, 
                     start_position, 
                     end_position)
head(ch)

gene_intervals <- getBM(
  attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position"),
  filters = "ensembl_gene_id",
  values = dnmt$ENSGID,
  mart = ensembl
)

dnmt <- merge(dnmt, gene_intervals, 
              by.x="ENSGID",
              by.y="ensembl_gene_id")

dnmt$region <- paste0(dnmt$chromosome_name,":",
                      dnmt$start_position,"-",
                      dnmt$end_position)

dnmt <- dnmt %>% arrange(chromosome_name, 
                     start_position, 
                     end_position)
head(dnmt)


gene_intervals <- getBM(
  attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position"),
  filters = "ensembl_gene_id",
  values = tet$ENSGID,
  mart = ensembl
)
tet <- merge(tet, gene_intervals, 
             by.x="ENSGID",
             by.y="ensembl_gene_id")

tet$region <- paste0(tet$chromosome_name,":",
                     tet$start_position,"-",
                     tet$end_position)

tet <- tet %>% arrange(chromosome_name, 
                         start_position, 
                         end_position)
head(tet)

## mash data ####
mash_chip.sorted <- fread("/Users/muddin/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/mash/annoated.sorted.chip_all.lfsr05.oct2024.csv", header=T)
head(mash_chip.sorted)

# Define ±500kb window
buffer_size <- 500000  # 500kb

# CHIP
# Identify genes within ±500kb of GWAS loci
selected_chip_genes <-  ch %>%
  inner_join(mash_chip.sorted, 
             by = c("chromosome_name" = "Otherinfo4"),
             relationship = "many-to-many") %>%  # Match chromosomes
  filter(start_position >= Otherinfo5 - buffer_size & end_position <= Otherinfo5 + buffer_size)  %>% # -+500kb window
  distinct(GENE, .keep_all = TRUE)

head(selected_chip_genes)

selected_chip_genes <- selected_chip_genes %>% arrange(chromosome_name, start_position, Otherinfo5)
head(selected_chip_genes)

fwrite(selected_chip_genes, "~/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/mash/chip_mash_selected_pops_genes.csv", row.names = F, col.names=T, quote = T, sep=","); rm(selected_chip_genes, mash_chip.sorted, ch)

# dnmt
mash_dnmt.sorted <- fread("/Users/muddin/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/mash/annoated.sorted.dnmt3a_all.lfsr05.oct2024.csv", header=T)
head(mash_dnmt.sorted)

# Define ±500kb window
buffer_size <- 500000  # 500kb

# CHIP
# Identify genes within ±500kb of GWAS loci
selected_dnmt_genes <-  dnmt %>%
  inner_join(mash_dnmt.sorted, 
             by = c("chromosome_name" = "Otherinfo4"),
             relationship = "many-to-many") %>%  # Match chromosomes
  filter(start_position >= Otherinfo5 - buffer_size & end_position <= Otherinfo5 + buffer_size)  %>% # -+500kb window
  distinct(GENE, .keep_all = TRUE)

head(selected_dnmt_genes)

selected_dnmt_genes <- selected_dnmt_genes %>% arrange(chromosome_name, start_position, Otherinfo5)
head(selected_dnmt_genes)

fwrite(selected_dnmt_genes, "~/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/mash/dnmt3a_mash_selected_pops_genes.csv", row.names = F, col.names=T, quote = T, sep=","); rm(selected_dnmt_genes, mash_dnmt.sorted, dnmt)

# tet
mash_tet.sorted <- fread("/Users/muddin/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/mash/annoated.sorted.tet2_all.lfsr05.oct2024.csv", header=T)
head(mash_tet.sorted)

# Define ±500kb window
buffer_size <- 500000  # 500kb

# CHIP
# Identify genes within ±500kb of GWAS loci
selected_tet_genes <-  tet %>%
  inner_join(mash_tet.sorted, 
             by = c("chromosome_name" = "Otherinfo4"),
             relationship = "many-to-many") %>%  # Match chromosomes
  filter(start_position >= Otherinfo5 - buffer_size & end_position <= Otherinfo5 + buffer_size)  %>% # -+500kb window
  distinct(GENE, .keep_all = TRUE)

head(selected_tet_genes)

selected_tet_genes <- selected_tet_genes %>% arrange(chromosome_name, start_position, Otherinfo5)
head(selected_tet_genes)

fwrite(selected_tet_genes, "~/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/mash/tet2_mash_selected_pops_genes.csv", row.names = F, col.names=T, quote = T, sep=","); rm(selected_tet_genes, tet, mash_tet.sorted)

## Select top 3 genes per locus #### 
library(dplyr)
mash_tet <- fread("~/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/mash/tet2_mash_selected_pops_genes.csv", header = T)

tet_top3 <- mash_tet %>%
  group_by(Chr, Start, End) %>%
  slice_max(order_by = PoPS_Score, n = 3, with_ties = FALSE) %>%
  ungroup()

fwrite(tet_top3, "~/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/mash/tet2_mash_top3pops_genes.csv", row.names = F, col.names=T, quote = T, sep=","); rm(tet_top3)

mash_dnmt <- fread("~/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/mash/dnmt3a_mash_selected_pops_genes.csv", header = T)

dnmt_top3 <- mash_dnmt %>%
  group_by(Chr, Start, End) %>%
  slice_max(order_by = PoPS_Score, n = 3, with_ties = FALSE) %>%
  ungroup()

fwrite(dnmt_top3, "~/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/mash/dnmt3a_mash_top3pops_genes.csv", row.names = F, col.names=T, quote = T, sep=","); rm(dnmt_top3)

mash_chip <- fread("~/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/mash/chip_mash_selected_pops_genes.csv", header = T)

chip_top3 <- mash_chip %>%
  group_by(Chr, Start, End) %>%
  slice_max(order_by = PoPS_Score, n = 3, with_ties = FALSE) %>%
  ungroup()

fwrite(chip_top3, "~/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/mash/chip_mash_top3pops_genes.csv", row.names = F, col.names=T, quote = T, sep=","); rm(chip_top3)

