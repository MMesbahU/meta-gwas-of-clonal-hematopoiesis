##### 01 - library #####
library(data.table)
library(survival)
library(dplyr)
library(optparse)
library(broom)
library(glue)

gc()
rm(list = ls())

setwd("/medpop/esp2/lli/chip_gwas_rev/")


##### 02 - files #####
results_full <- fread("Results/ukb_phewas/ukbb_phewas_results_full.tsv", data.table = FALSE)
phecodes_list <- fread("Data/ukb_phewas/phecode_list.tsv", data.table = F, header=F)$V1

# i=1
incd_full <- data.frame()
for (i in seq_along(phecodes_list)){
    phecode <- phecodes_list[i]
    file_name <- glue("Results/ukb_phewas/Counts_ukbb_mesbahCHIP_Phewas.PheCode{phecode}.tsv")
    
    # only EUR
    dat <- fread(file_name, data.table = F) %>%
        mutate(phecode_cat = phecode, .before=1)
    
    incd_full <- rbind(incd_full, dat)
}



write.table(incd_full, "Results/ukb_phewas/ukbb_phewas_incd_counts_full.tsv", sep="\t", row.names = F, quote = F)


intersect(names(incd_full), names(results_full))




##### 03 - merge for mesbah #####
results_n_incd_counts <- merge(incd_full, results_full, by = c("phecode_cat", "outcomes", "exposures", "knn"))

results_n_incd_counts_final <- results_n_incd_counts %>%
    select(-n_event, -n)
    

# test <- results_n_incd_counts %>%
#     filter(!is.na(n)) %>%
#     filter((n_incd_cases + n_controls) != n) # 0 rows, perfect


write.table(results_n_incd_counts_final, "Results/ukb_phewas/ukbb_phewas_results_n_counts_full.tsv",
            sep="\t", row.names = F, quote = F)




