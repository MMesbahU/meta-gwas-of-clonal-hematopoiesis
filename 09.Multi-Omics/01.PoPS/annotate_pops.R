## Annotate PoPS #### 
########## Jan 24, 2024
## POPs score
# BiocManager::install("biomaRt")
library(biomaRt)
library(data.table)
require(openxlsx)
setwd("~/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/pops/") 
ch <- fread("CH_pops_gene.preds.txt", header = T, na.strings = "NA")
dnmt <- fread("DNMT3A_pops_gene.preds.txt", header = T, na.strings = "NA")
tet <- fread("TET2_pops_gene.preds.txt", header = T, na.strings = "NA")

######
# Select the human genes dataset
ensembl <- useEnsembl(biomart = "ensembl", 
                      dataset = "hsapiens_gene_ensembl", 
                      mirror = "useast", GRCh=37)
######
# Get the gene intervals
gene_intervals <- getBM(
  attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position"),
  filters = "ensembl_gene_id",
  values = ch$ENSGID,
  mart = ensembl
)
#########
table(gene_intervals$ensembl_gene_id %in% ch$ENSGID)

ch <- merge(ch, gene_intervals, 
            by.x="ENSGID",
            by.y="ensembl_gene_id")

gene_intervals <- getBM(
  attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position"),
  filters = "ensembl_gene_id",
  values = dnmt$ENSGID,
  mart = ensembl
)

dnmt <- merge(dnmt, gene_intervals, 
              by.x="ENSGID",
              by.y="ensembl_gene_id")

tet <- merge(tet, gene_intervals, 
             by.x="ENSGID",
             by.y="ensembl_gene_id")

########
pops_cdt <- list("Table S1 CHIP"=ch, 
                 "Table S2 DNMT3A"=dnmt,
                 "Table S3 TET2"=tet)

write.xlsx(pops_cdt, 
           file = "Table_S1_3.chr_pos.POPs.metaCHIP.N900k.Jan2024.xlsx")
#-------------------------------------------------------
#Read files
library(readxl)
library(dplyr)
library(tidyr)
chip <- read_excel("Table_S1_3.chr_pos.POPs.metaCHIP.N900k.Jan2024.xlsx", 
                   sheet = "Table S1 CHIP")
chip$region <- paste0("chr", chip$chromosome_name,":",
                      chip$start_position,"-",
                      chip$end_position)

chip_d <- chip %>% 
  filter(! (chromosome_name %in% c("HSCHR4_6_CTG12", "KI270713.1") ) ) %>% 
  mutate(chrom = as.numeric(chromosome_name)) %>%
  arrange(chrom, start_position, end_position)



dnmt <- read_excel("Table_S1_3.chr_pos.POPs.metaCHIP.N900k.Jan2024.xlsx", 
                   sheet = "Table S2 DNMT3A")
dnmt$region <- paste0("chr", dnmt$chromosome_name,":",
                      dnmt$start_position,"-",
                      dnmt$end_position)
dnmt_d <- dnmt %>% 
  filter(! (chromosome_name %in% c("HSCHR4_6_CTG12", "KI270713.1") ) ) %>% 
  mutate(chrom = as.numeric(chromosome_name)) %>%
  arrange(chrom, start_position, end_position)

tet2 <- read_excel("Table_S1_3.chr_pos.POPs.metaCHIP.N900k.Jan2024.xlsx", 
                   sheet = "Table S3 TET2")
tet2$region <- paste0("chr", tet2$chromosome_name,":",
                      tet2$start_position,"-",
                      tet2$end_position)
tet2_d <- tet2 %>% 
  filter(! (chromosome_name %in% c("HSCHR4_6_CTG12", "KI270713.1") ) ) %>% 
  mutate(chrom = as.numeric(chromosome_name)) %>%
  arrange(chrom, start_position, end_position)

##
sorted.pops_cdt <- list("Table S1 CHIP"=chip_d, 
                        "Table S2 DNMT3A"=dnmt_d,
                        "Table S3 TET2"=tet2_d)
write.xlsx(sorted.pops_cdt, 
           file = "Table_S1_3.POPs.metaCHIP.N900k.Jan2025.sorted.xlsx")
#-------------------------------------------------------
#### GWAS loci +-500kb
chip_loci <- fread("~/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/Summary/Final2024/summary/annovar/eaf001/CHIP/eaf001p5e7.annovar.pval_0.0001.GWAMA.chr1_22.hasCH.MultiANC.ukbb200k_ukb250k_AoU250k_TOPMed72k_MGBB53k_BioVU54k.csv")
summary(chip_loci$AAF)
summary(chip_loci$N_Studies)
summary(chip_loci$P)

chip_loci <- chip_loci %>% filter(N_Studies>=2 & P<=5e-8) 
chip_d <- chip_d %>% mutate(chromosome_name = paste0("chr", chromosome_name))
head(chip_d)

# Define ±500kb window
buffer_size <- 500000  # 500kb

# CHIP
# Identify genes within ±500kb of GWAS loci
selected_chip_genes <- chip_d %>%
  inner_join(chip_loci, by = c("chromosome_name" = "Otherinfo4")) %>%  # Match chromosomes
  filter(start_position <= Otherinfo5 + buffer_size & end_position >= Otherinfo5 - buffer_size)  %>% # ±500kb window
  distinct(GENE, .keep_all = TRUE)

# dnmt3a
dnmt3a_loci <- fread("~/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/Summary/Final2024/summary/annovar/eaf001/DNMT3A/eaf001p5e7.annovar.pval_0.0001.GWAMA.chr1_22.hasDNMT3A.MultiANC.ukbb200k_ukb250k_AoU250k_TOPMed72k_MGBB53k_BioVU54k.csv")
summary(dnmt3a_loci$AAF)
summary(dnmt3a_loci$N_Studies)
summary(dnmt3a_loci$P)

dnmt3a_loci <- dnmt3a_loci %>% filter(N_Studies>=2 & P<=5e-8) 
dnmt_d <- dnmt_d %>% mutate(chromosome_name = paste0("chr", chromosome_name))
head(dnmt_d)
# Define ±500kb window
buffer_size <- 500000  # 500kb
# Identify genes within ±500kb of GWAS loci
selected_dnmt_genes <- dnmt_d %>%
  inner_join(dnmt3a_loci, by = c("chromosome_name" = "Otherinfo4")) %>%  # Match chromosomes
  filter(start_position <= Otherinfo5 + buffer_size & end_position >= Otherinfo5 - buffer_size)  %>% # ±500kb window
  distinct(GENE, .keep_all = TRUE)


# tet2
tet2_loci <- fread("~/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/Summary/Final2024/summary/annovar/eaf001/TET2/eaf001p5e7.annovar.pval_0.0001.GWAMA.chr1_22.hasTET2.MultiANC.ukbb200k_ukb250k_AoU250k_TOPMed72k_MGBB53k_BioVU54k.csv")
summary(tet2_loci$AAF)
summary(tet2_loci$N_Studies)
summary(tet2_loci$P)

tet2_loci <- tet2_loci %>% filter(N_Studies>=2 & P<=5e-8) 
tet2_d <- tet2_d %>% mutate(chromosome_name = paste0("chr", chromosome_name))
head(tet2_d)

# Define ±500kb window
buffer_size <- 500000  # 500kb
# Identify genes within ±500kb of GWAS loci
selected_tet2_genes <- tet2_d %>%
  inner_join(tet2_loci, by = c("chromosome_name" = "Otherinfo4")) %>%  # Match chromosomes
  filter(start_position <= Otherinfo5 + buffer_size & end_position >= Otherinfo5 - buffer_size)  %>% # ±500kb window
  distinct(GENE, .keep_all = TRUE)

## save
##
gwasloci.pops_cdt <- list("Table S1 CHIP"=selected_chip_genes, 
                          "Table S2 DNMT3A"=selected_dnmt_genes,
                          "Table S3 TET2"=selected_tet2_genes)

# write.xlsx(gwasloci.pops_cdt, 
#            file = "Table_S1_3.POPs.metaCHIP.N900k.Jan2025.sorted_ALLgwasVars.xlsx")

write.xlsx(gwasloci.pops_cdt, 
           file = "Table_S1_3.POPs.metaCHIP.N900k.Jan2025.sorted_gwas_loci.xlsx")
######################################
