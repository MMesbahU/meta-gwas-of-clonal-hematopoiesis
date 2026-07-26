###### 01 - library #####
library(dplyr)
library(data.table)
library(optparse)
library(glue)
library(survival)

gc()
rm(list=ls())

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


##### 02 - files #####
load(file='/medpop/esp2/mesbah/projects/Meta_GWAS/n650k/chip_call/ukb450k_mgbb53k.CHIP_vars.for_phewas.rda')
rm(ukb450k_var, ukb450k_ch, mgbb_var)

mgbb_ch <- mgbb_ch %>%
    rename(id = Biobank_Subject_ID)

covs_df <- fread("Data/mgb_phewas/mgbb_covs_final.tsv", data.table = F)
cytopenia <- fread("Data/mgb_phewas/mgbb_cytopenia_final.tsv", data.table = F)

# phenCategory = 'HematologicNeoplasmInfectious'
phecodes <- fread(
    glue("/medpop/esp2/mzekavat/Partners_BB/phenos/PheCodes/clean_mgbb.ICDdiagnosis.phecodes2020.{phenCategory}.20201105_gender.restricted.csv"), 
        data.table = F)

phecodes <- phecodes %>%
    select(id = Biobank.Subject.ID, matches("_(ANY|PREV|DAYS|INCID)$")) %>%
    select(-contains("Tobacco_use_disorder"))

pheno_list <- unique(colnames(phecodes)[grep("_ANY", colnames(phecodes))])
outcome_vars <- gsub("_ANY","",pheno_list)
    

exposure_vars <- c('CHIP', 'DNMT3A', 'TET2', 'ASXL1', 'JAK2', 'DDR', 'Splice', 
                   'expandedCHIP', 'expandedDNMT3A', 'expandedTET2', 'expandedASXL1', 'expandedJAK2', 
                   'expandedDDR', 'expandedSplice')
mgbb_ch <- mgbb_ch %>%
    select(id, all_of(exposure_vars))

ancestries <- names(table(covs_df$knn))

covs <- paste(setdiff(names(covs_df), c("id", "Prev_Heme_CA", "knn")), collapse ="+")
cov_cols <- setdiff(names(covs_df), c("id", "Prev_Heme_CA", "knn"))



## merging master df
master_df <- left_join(cytopenia, covs_df, by ="id") %>%
    left_join(phecodes, by = "id") %>%
    left_join(mgbb_ch, by = "id") %>%
    filter(any_cytopenia == 0) %>%
    select(-c(any_cytopenia, Prev_Heme_CA))




##### 03 - looping for analysis #####
results <- data.frame()

# i=j=k=1
for (i in 1:length(outcome_vars)){
    
    disease <- outcome_vars[i]
    
    prev_col <- glue("{disease}_PREV")
    incid_col <- glue("{disease}_INCID")
    
    ## filtering out previous cases
    temp <- master_df %>%
        filter(.data[[prev_col]] != 1)
    
    
    for (j in 1:length(exposure_vars)){
        
        exposure_var <- exposure_vars[j]
        # running survival analysis
        fmla <- as.formula(glue("Surv({disease}_DAYS, {disease}_INCID) ~ {exposure_var} + {covs}"))
        
        for (k in 1:length(ancestries)){
            
            ancestry = ancestries[k]
            
            print(glue("Disease: {disease} | Exposure: {exposure_var} | Ancestry: {ancestry}"))
            
            # ancestry stratified
            temp1 <- temp[temp$knn == ancestry,] %>%
                select(contains(disease), !!exposure_var, all_of(cov_cols)) %>%
                tidyr::drop_na(.)
            
            # default NA result (will be overwritten only if coxph succeeds)
            result <- data.frame(
                outcomes = disease,
                exposures = exposure_var,
                knn = ancestry,
                n = NA, n_event = NA,
                hr = NA, se = NA, pval = NA,
                lci = NA, uci = NA,
                n_incd_cases = NA, n_controls = NA,
                n_incd_cases_chip = NA, n_controls_chip = NA
            )
            
            if (nrow(temp1) > 0) {
                
                ## getting incidence cases
                incid_vec <- temp1[[incid_col]]
                chip_vec  <- temp1[[exposure_var]]
                
                n_incd_cases      <- sum(incid_vec == 1, na.rm = TRUE)
                n_controls        <- sum(incid_vec == 0, na.rm = TRUE)
                n_incd_cases_chip <- sum(incid_vec == 1 & chip_vec == 1, na.rm = TRUE)
                n_controls_chip   <- sum(incid_vec == 0 & chip_vec == 1, na.rm = TRUE)
                
                reg <- tryCatch(
                    coxph(fmla, data = temp1),
                    error = function(e) NULL
                )
                
                if (!is.null(reg) && !any(is.na(coef(reg))) && all(is.finite(coef(reg)))) {
                    
                    tidy_res <- broom::tidy(reg, conf.int = TRUE, exponentiate = TRUE) %>%
                        filter(term == exposure_var)
                    
                    result <- tidy_res %>%
                        mutate(
                            outcomes = disease,
                            exposures = exposure_var,
                            knn = ancestry,
                            n = reg$n,
                            n_event = reg$nevent,
                            n_incd_cases      = n_incd_cases,
                            n_controls        = n_controls,
                            n_incd_cases_chip = n_incd_cases_chip,
                            n_controls_chip   = n_controls_chip
                        ) %>%
                        select(
                            outcomes, exposures, knn, n, n_event,
                            hr = estimate, se = std.error, pval = p.value,
                            lci = conf.low, uci = conf.high,
                            n_incd_cases, n_controls, n_incd_cases_chip, n_controls_chip
                        )
                    
                } else {
                    # model failed, but still record the incidence counts we computed
                    result$n_incd_cases      <- n_incd_cases
                    result$n_controls        <- n_controls
                    result$n_incd_cases_chip <- n_incd_cases_chip
                    result$n_controls_chip   <- n_controls_chip
                }
            }
            
            ## append (always, even if NA)
            results <- rbind(results, result)
            
            write.table(results, 
                        glue("Results/mgbb_phewas/HR_mgbb_mesbahCHIP_Phewas_CoxPH.PheCode{phenCategory}.tsv"),
                        row.names = FALSE, quote = FALSE, sep = "\t")
        }
    }
}

print("all done")

write.table(results, glue("Results/mgbb_phewas/HR_mgbb_mesbahCHIP_Phewas_CoxPH.PheCode{phenCategory}.tsv"), 
            row.names = F, quote = F, sep = "\t")




