##### 01 - library #####
library(dplyr)
library(data.table)

gc()
rm(list=ls())


setwd("~/linke/chip_gwas_rev/")



##### 02 - files #####
chip_raw <- fread("Data/MAGMA_files/ukb200k_eur/magma_out.ukb_eur.chr1_22.hasCHIP.genes.out", data.table = F)

dnmt3a_raw <- fread("Data/MAGMA_files/ukb200k_eur/magma_out.ukb_eur.chr1_22.hasDNMT3A.genes.out", data.table = F)
tet2_raw <- fread("Data/MAGMA_files/ukb200k_eur/magma_out.ukb_eur.chr1_22.hasTET2.genes.out", data.table = F)

annot <- fread("Data/MAGMA_files/ukb200k_eur/magma_genes_gprofiler_annot.csv", data.table = F)
annot_short <- annot %>%
    select(GENE=initial_alias, gene_symbol = name) 


## renaming cols
chip_short <- chip_raw %>%
    select(GENE, CHIP=ZSTAT)
dnmt3a_short <- dnmt3a_raw %>%
    select(GENE, DNMT3A=ZSTAT)
tet2_short <- tet2_raw %>%
    select(GENE, TET2=ZSTAT)

full <- left_join(chip_short, dnmt3a_short, by = "GENE") %>%
    left_join(tet2_short, by = "GENE") %>%
    left_join(annot_short, by = "GENE") %>%
    filter(GENE %notin% c("ENSG00000258724","ENSG00000105501")) %>%
    select(-GENE) %>%
    rename(GENE = gene_symbol) %>%
    relocate(GENE, .before=1) %>%
    filter(GENE != "")


full_ensg <- left_join(chip_short, dnmt3a_short, by = "GENE") %>%
    left_join(tet2_short, by = "GENE")

write.table(full, "Data/MAGMA_files/magma_chip_dnmt3a_tet2_merged.tsv",
            sep="\t", row.names = F, quote = F)
write.table(full_ensg, "Data/MAGMA_files/magma_chip_dnmt3a_tet2_merged_ensg.tsv",
            sep="\t", row.names = F, quote = F)





# subsampling to 500 genes ------------------------------------------------
full <- left_join(chip_short, dnmt3a_short, by = "GENE") %>%
    left_join(tet2_short, by = "GENE") %>%
    left_join(annot_short, by = "GENE") %>%
    filter(GENE %notin% c("ENSG00000258724","ENSG00000105501")) %>%
    select(-GENE) %>%
    rename(GENE = gene_symbol) %>%
    relocate(GENE, .before=1) %>%
    filter(GENE != "")


full_ensg <- left_join(chip_short, dnmt3a_short, by = "GENE") %>%
    left_join(tet2_short, by = "GENE")

write.table(full, "Data/MAGMA_files/magma_chip_dnmt3a_tet2_merged.tsv",
            sep="\t", row.names = F, quote = F)
write.table(full_ensg, "Data/MAGMA_files/magma_chip_dnmt3a_tet2_merged_ensg.tsv",
            sep="\t", row.names = F, quote = F)




