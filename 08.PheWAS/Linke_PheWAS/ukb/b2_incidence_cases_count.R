##### 01 - library #####
library(data.table)
library(survival)
library(dplyr)
library(optparse)
library(broom)
library(glue)
library(progress)

gc()
rm(list = ls())

setwd("/medpop/esp2/lli/chip_gwas_rev/")


### optparse
option_list <- list(
    make_option(
        c("--phenCategory"), type = "character",
        default = NULL, help = "start column (characters)", metavar = "character"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
phenCategory <- opt$phenCategory



###### 02 - files ######
# mesbah's CHIP file
load(file='/medpop/esp2/mesbah/projects/Meta_GWAS/n650k/chip_call/ukb450k_mgbb53k.CHIP_vars.for_phewas.rda')
rm(ukb450k_var, mgbb_ch, mgbb_var)

knn_df <- fread("/medpop/esp2/skoyama/passing/ukb_kgpprojection/v03/out/ukb.kgp_projected.tsv.gz", data.table = F)
knn_df <- knn_df %>%
    select(id = eid, knn) %>%
    mutate(id = as.numeric(id)) %>%
    filter(!is.na(id))

covs_df <- fread("Data/ukb_phewas/ukb_phewas_covs_final.tsv", data.table = F)

phecodes = fread(
    glue("/medpop/esp2/mzekavat/UKBB/PhenoFiles/PheCODEs_new/2021-01-08_ukb_phecode_{phenCategory}_March2020fu.csv"), 
    data.table = F)

phecodes <- phecodes %>%
    select(id = f.eid, matches("_(ANY|PREV|DAYS|INCID)$")) %>%
    select(-contains("Tobacco_use_disorder"))

# get a disease list within each phecode file for looping
outcome_cols <- unique(colnames(phecodes)[grep("_ANY", colnames(phecodes))])
outcome_vars <- gsub("_ANY","", outcome_cols)

exposure_vars <- c('CHIP', 'DNMT3A', 'TET2', 'ASXL1', 'JAK2', 'DDR', 'Splice', 
                   'expandedCHIP', 'expandedDNMT3A', 'expandedTET2', 'expandedASXL1', 
                   'expandedJAK2', 'expandedDDR', 'expandedSplice')

ancestries <- names(table(knn_df$knn))


## merging master df
master_df <- left_join(ukb450k_ch %>% rename(id = eid_7089), covs_df, by ="id") %>%
    left_join(phecodes, by = "id") %>%
    left_join(knn_df, by = "id") %>%
    filter(Prev_Maryam_HemeCa_Phenos != 1, cytopenia == 0) %>%
    select(-c(Prev_Maryam_HemeCa_Phenos, cytopenia))

cov_cols <- setdiff(names(covs_df), c("id", "Prev_Maryam_HemeCa_Phenos", "cytopenia"))


##### 03 - counts only (no model fitting) #####
counts <- data.frame()
pb <- progress_bar$new(
    format = "[:bar] :current/:total (:percent) | ETA: :eta | elapsed: :elapsed",
    total = length(outcome_vars) * length(exposure_vars) * length(ancestries),
    clear = FALSE,          # keep the final bar
    width = 60,
    show_after = 0
)

for (i in 1:length(outcome_vars)){
    
    disease <- outcome_vars[i]
    
    prev_col  <- glue("{disease}_PREV")
    incid_col <- glue("{disease}_INCID")
    
    temp <- master_df %>%
        filter(.data[[prev_col]] != 1)
    
    for (j in 1:length(exposure_vars)){
        
        exposure_var <- exposure_vars[j]
        
        for (k in 1:length(ancestries)){
            
            pb$tick()
            
            ancestry <- ancestries[k]
            
            print(glue("Counts | Disease: {disease} | Exposure: {exposure_var} | Ancestry: {ancestry}"))
            
            temp1 <- temp[temp$knn == ancestry,] %>%
                select(contains(disease), !!exposure_var, all_of(cov_cols)) %>%
                tidyr::drop_na(.)
            
            count_result <- data.frame(
                outcomes = disease,
                exposures = exposure_var,
                knn = ancestry,
                n_incd_cases      = NA_integer_,
                n_controls        = NA_integer_,
                n_incd_cases_chip = NA_integer_,
                n_controls_chip   = NA_integer_
            )
            
            if (nrow(temp1) > 0) {
                incid_vec <- temp1[[incid_col]]
                chip_vec  <- temp1[[exposure_var]]
                
                count_result$n_incd_cases      <- sum(incid_vec == 1, na.rm = TRUE)
                count_result$n_controls        <- sum(incid_vec == 0, na.rm = TRUE)
                count_result$n_incd_cases_chip <- sum(incid_vec == 1 & chip_vec == 1, na.rm = TRUE)
                count_result$n_controls_chip   <- sum(incid_vec == 0 & chip_vec == 1, na.rm = TRUE)
            }
            
            counts <- rbind(counts, count_result)
            
            write.table(counts, 
                        glue("Results/ukb_phewas/Counts_ukbb_mesbahCHIP_Phewas.PheCode{phenCategory}.tsv"),
                        row.names = FALSE, quote = FALSE, sep = "\t")
        }
    }
}

print("all done - counts")

write.table(counts, glue("Results/ukb_phewas/Counts_ukbb_mesbahCHIP_Phewas.PheCode{phenCategory}.tsv"), 
            row.names = F, quote = F, sep = "\t")


