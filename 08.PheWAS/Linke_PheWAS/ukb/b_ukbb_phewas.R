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

### testing townsend section
knn_townsend <- left_join(knn_df, covs_df, by ="id") %>%
    mutate(townsend_exist_flag = if_else(!is.na(townsend),1,0))
table(knn_townsend$knn, knn_townsend$townsend_exist_flag)
rm(knn_townsend)


# testing section
# phenCategory = 'DermDigestGU'
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

covs <- paste(setdiff(names(covs_df), c("id", "Prev_Maryam_HemeCa_Phenos", "cytopenia")), collapse = "+")
cov_cols <- setdiff(names(covs_df), c("id", "Prev_Maryam_HemeCa_Phenos", "cytopenia"))

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
                lci = NA, uci = NA
            )
            
            if (nrow(temp1) > 0) {
                
                reg <- tryCatch(
                    coxph(fmla, data = temp1),
                    error = function(e) NULL
                )
                
                # only extract results when the model actually fitted
                if (!is.null(reg) && !any(is.na(coef(reg))) && all(is.finite(coef(reg)))) {
                    
                    tidy_res <- broom::tidy(reg, conf.int = TRUE, exponentiate = TRUE) %>%
                        filter(term == exposure_var)
                    
                    result <- tidy_res %>%
                        mutate(
                            outcomes = disease,
                            exposures = exposure_var,
                            knn = ancestry,
                            n = reg$n,
                            n_event = reg$nevent
                        ) %>%
                        select(
                            outcomes, exposures, knn, n, n_event,
                            hr = estimate, se = std.error, pval = p.value,
                            lci = conf.low, uci = conf.high
                        )
                }
            }
            
            ## append (always, even if NA)
            results <- rbind(results, result)
            
            write.table(results, 
                        glue("Results/ukb_phewas/HR_ukbb_mesbahCHIP_Phewas_CoxPH.PheCode{phenCategory}.tsv"),
                        row.names = FALSE, quote = FALSE, sep = "\t")
        }
    }
}

print("all done")

write.table(results, glue("Results/ukb_phewas/HR_ukbb_mesbahCHIP_Phewas_CoxPH.PheCode{phenCategory}.tsv"), 
            row.names = F, quote = F, sep = "\t")

