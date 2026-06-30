###############################################################################
#             0) Clean workspace and load necessary libraries
###############################################################################
rm(list = ls())         # clear everything
gc()                    # force garbage collection

# Libraries
library(data.table)     # for fread, data.table manipulations
library(dplyr) 
library(survival)       # for coxph
library(fst)            # fast I/O to disk
library(future.apply)   # parallel 'future_lapply'
library(progressr)      # nice progress bars for parallel tasks

###############################################################################
##
ARGs <- commandArgs(TRUE)
outdir <- ARGs[1]
Phecode_Group <- ARGs[2]
##
# setwd("/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/PheWAS/")
setwd(outdir)
#########

###############################################################################
#             1) Read and merge your UKB datasets + PheCODE outcomes
#                (Replace these paths and merges with your own)
###############################################################################

cat("1) Loading/merging large data...\n")
#### CHIP data in UKB
# # CHIP VAR
# ukbbch_var <- fread("/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/ch_var_all_ukb200k_n_ukb250k.7Mar2023.csv.gz", header=T)
# head(ukbbch_var)
# sort(table(ukbbch_var$Gene), decreasing = T)
# 
# sort(table(ukbbch_var$Protein_Change[ukbbch_var$Gene%in% c("IDH1", "IDH2")]), decreasing = T)

# CHIP: /medpop/esp2/mesbah/tools/CHIP_metaAnalysis//04.SNP_heritability/QuickPRS
ukb200k_chip <- fread("/medpop/esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/CH_phenoCovar.ukb200k_N193342.28cols.03_05_2024.tsv.gz", header=T)

ukb250k_chip <- fread("/medpop/esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/CH_phenoCovar.ukb250k_N243350.28cols.03_05_2024.tsv.gz", header=T)
table(ukb200k_chip$FID %in% ukb250k_chip$FID)
ukb450k_chip <- as.data.frame(rbind(ukb200k_chip, ukb250k_chip)); rm(ukb200k_chip, ukb250k_chip)

cat("unique sample size", length(unique(ukb450k_chip$FID)))

head(ukb450k_chip)  

table(ukb450k_chip$Batch, exclude = NULL)

## UKB baseline info 
## UKB Baseline
d_base <- fread("/medpop/esp2/mesbah/projects/gxe/UKB_baseline/UKB500k_baseline_info_recoded_participant.tsv.gz")
# "ever_smoked.p20160_i0",         
# "Smoking_status.p20116_i0",
table(d_base$p20116_i0,  d_base$p20160_i0, exclude = NULL)
## Smoking status
ukb450k_chip$smking_status <- factor(ifelse(ukb450k_chip$FID %in% d_base$eid[d_base$p20160_i0=="Yes"], 
                                            1, 
                                            ifelse(ukb450k_chip$FID %in% d_base$eid[d_base$p20160_i0=="No"],
                                                   0, NA)))
table(ukb450k_chip$smking_status, exclude = NULL)

rm(d_base)

##### regroup
table(ukb450k_chip$knn, exclude = NULL)
table(ukb450k_chip$Batch, exclude = NULL)
table(ukb450k_chip$Genetic_Sex, exclude = NULL)
table(ukb450k_chip$GenoBatch, exclude = NULL)
table(ukb450k_chip$smking_status, exclude = NULL)

ukb450k_chip$Genetic_Ancestry_group <- case_when(
    ukb450k_chip$knn == "EUR" ~ "EUR", 
    ukb450k_chip$knn == "AFR" ~ "AFR", 
    ukb450k_chip$knn == "SAS" ~ "SAS", 
TRUE ~ "AMR_or_EAS")

names(ukb450k_chip)
ls()

table(ukb450k_chip$knn, ukb450k_chip$Genetic_Ancestry_group, exclude = NULL)


######
## PRS with EUR Reference Panel
load("/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/prs/prs.eur.ukb450k.rda")

ukb450k_chip_eur <- merge(ukb450k_chip, dat_eur, 
                          by.x="FID", by.y="ID2"); rm(ukb450k_chip, dat_eur)
head(ukb450k_chip_eur)
names(ukb450k_chip_eur)
# multianc gwas prs
# ukb450k_chip_eur_multianc <- ukb450k_chip_eur[,c(1,2,13:21,27:29,35,38,44,50,32,41,47)]
# cat("unique sample size", length(unique(ukb450k_chip_eur_multianc$FID)))
# names(ukb450k_chip_eur_multianc)
### EUR-ancestry gwas PRS
# ukb450k_chip_eur_eurprs <- ukb450k_chip_eur[,c(1,2,13:21,27:29,34,37,43,49,31,40,46)]
# cat("unique sample size", length(unique(ukb450k_chip_eur_eurprs$FID)))
# names(ukb450k_chip_eur_eurprs)
# 
ukb450k_chip_eur.multi_n_eur_prs <- ukb450k_chip_eur[,c(1,2,13:21,27,29,30,35,36,
                                                        38,39,44,45,50,51,32,33,
                                                        41,42,47,48)]
cat("unique sample size", length(unique(ukb450k_chip_eur.multi_n_eur_prs$FID)))

names(ukb450k_chip_eur.multi_n_eur_prs)



# Phecode categories
# phenCategory <- c("CirculatoryRespiratory", 
  #                "HematologicNeoplasmInfectious",
   #               "DermDigestGU",
    #              "MentalNeuroSensorySymptoms", 
     #             "InjuriesPoisoningMSK",
      #            "PregnancyCongenitalEndocrine")


# for(k in 1:length(phenCategory) ){
  gc()
  # cat(k,"\n")
  # cat(phenCategory[k],"\n")

cat(Phecode_Group, "\n")

dat <- fread(paste0("/medpop/esp2/mzekavat/UKBB/PhenoFiles/PheCODEs_new/2021-01-08_ukb_phecode_", Phecode_Group, "_March2020fu.csv"), sep=",")
#  dat <- fread(paste("/medpop/esp2/mzekavat/UKBB/PhenoFiles/PheCODEs_new/2021-01-08_ukb_phecode_", phenCategory[k], "_March2020fu.csv", sep=""), sep=",")
  # dat <- fread(paste("../../../PheWAS/UKBB_MGBB_MVP_BioVU/PheCODEs_new/2021-01-08_ukb_phecode_", phenCategory[k], "_March2020fu.csv.gz", sep=""), sep=",")
  
  dat <- dat %>% select(1, contains("_INCID"), contains("_DAYS"))
  
  # Multi-ancestry and EUR GWAS PRS with EUR ref
  ukb_prs_multianc <- merge(ukb450k_chip_eur.multi_n_eur_prs, 
                            dat,
                            by.x="FID", 
                            by.y="f.eid")
  rm(dat)
  gc()
  
   
  cat("Data dimension:", dim(ukb_prs_multianc), "\n")
  
  
  ###############################################################################
  #             2) Identify your outcomes and exposures
  ###############################################################################
  # E.g. you found outcomes by searching for "_INCID" columns, then removing suffix:
  outcomes <- gsub(
    pattern     = "_INCID", 
    replacement = "", 
    x           = names(ukb_prs_multianc)[grepl("INCID", names(ukb_prs_multianc), ignore.case = TRUE)]
  )
  cat("Number of outcomes found:", length(outcomes), "\n")
  
  # 7 exposures from your code:
  # Define exposures
  # exposures <- c("CH_MultiANC","CHvaf10_MultiANC","DNMT3A_MultiANC",
    #             "TET2_MultiANC","ASXL1_MultiANC",
     #            "DDR_MultiANC","SF_MultiANC")
    
  exposures <- c('CH_EUR', 'CH_MultiANC', 'CHvaf10_EUR', 'CHvaf10_MultiANC', 
                 'DNMT3A_EUR', 'DNMT3A_MultiANC', 'TET2_EUR', 
                 'TET2_MultiANC', 'ASXL1_EUR', 'ASXL1_MultiANC', 
                 'DDR_EUR', 'DDR_MultiANC', 'SF_EUR', 'SF_MultiANC')
  
  
  cat("Number of exposures:", length(exposures), "\n")
  
  cat("Total # of Cox fits to run =", length(outcomes)*length(exposures), "\n")
  
    ###############################################################################
  #             3) Write to .fst (fast on-disk format), then remove from memory
  ###############################################################################
  # This avoids sending a giant data object to each worker.
  
    library(fst)
  
  
    cat("Saving ukb_prs_multianc.fst to disk...\n")
  
    write_fst(ukb_prs_multianc, paste0(Phecode_Group,".ukb_prs_multianc.fst"))
  
    rm(ukb_prs_multianc)
  
    gc()
  
  ###############################################################################
  #             4) Chunk the 224 outcomes so we don't do them all at once
  ###############################################################################
  # We'll process them in sets of e.g. 20 outcomes per chunk.
  
    chunk_size <- 20
  
    outcome_indices <- seq(1, length(outcomes), by = chunk_size)
  
  ###############################################################################
  #             5) Define a function to run Cox models for a single outcome
  #                on your 7 exposures
  ###############################################################################
  run_coxph_for_one_outcome <- function(outcome_prefix, exposures, progressor = NULL) {
    # 1) Read from disk inside the worker
    local_data <- as.data.table(read_fst(paste0(Phecode_Group,".ukb_prs_multianc.fst")) )
    
    # 2) Subset rows, e.g., exclude first 30 days
    outcome_incid <- paste0(outcome_prefix, "_INCID")
    outcome_days  <- paste0(outcome_prefix, "_DAYS")
    local_data <- local_data[get(outcome_days) > 30]
    
    # 3) Drop unused factor levels
    local_data <- droplevels(local_data)
    
    # 4) Check if any factor has <2 levels
   # factor_vars <- c("knn", "Genetic_Sex", "smking_status", "Batch", "GenoBatch")
     factor_vars <- c("Genetic_Ancestry_group", "Genetic_Sex") 
    for (fv in factor_vars) {
      if (fv %in% names(local_data)) {
        if (length(unique(local_data[[fv]])) < 2) {
          # Skip this outcome to avoid the contrasts error
          return(NULL)
        }
      }
    }
    
    # 5) Count events
    disease_cases <- sum(local_data[[outcome_incid]] == 1, na.rm = TRUE)
    total_n       <- nrow(local_data)
    if (disease_cases == 0 || total_n == 0) {
      return(NULL)
    }
    
    # 6) Loop over exposures, fit cox, store results
    results_for_this_outcome <- list()
    for (expo in exposures) {
      if (!is.null(progressor)) {
        progressor(sprintf("Outcome=%s, Exposure=%s", outcome_prefix, expo))
      }
      
      fml <- as.formula(
        paste0(
          "Surv(", outcome_days, ", ", outcome_incid, ") ~ scale(", expo, 
          ") + Genetic_Ancestry_group + Genetic_Sex + Age_at_recruitment + smking_status + ",
          "Batch + GenoBatch + sqrtAge_at_recruitment + PC1 + PC2 + PC3 + PC4 + PC5"
        )
      )
      
      model <- coxph(fml, data = local_data)
      coef_data <- summary(model)$coefficients[1, c(1, 3, 4, 5)]  # Beta, SE, Z, P
      
      results_for_this_outcome[[length(results_for_this_outcome) + 1]] <- data.table(
        Outcome       = outcome_prefix,
        Exposure      = expo,
        Beta          = coef_data[1],
        SE            = coef_data[2],
        Z             = coef_data[3],
        P             = coef_data[4],
        Disease_Cases = disease_cases,
        N             = total_n
      )
    }
    
    return(rbindlist(results_for_this_outcome))
  }
  
  ###############################################################################
  #             6) Main loop over outcome CHUNKS (parallel within each chunk)
  ###############################################################################
  # We'll accumulate results in a list. Each chunk returns a data.table of results.
  all_chunks_results <- list()
  chunk_counter <- 1
  
  # Set up parallel plan (keep it modest if memory is a concern)
  library(future)
  plan(multisession, workers = 5)
  
  cat("\nStarting chunked analysis...\n")
  
  for (start_idx in outcome_indices) {
    end_idx <- min(start_idx + chunk_size - 1, length(outcomes))
    chunk_outcomes <- outcomes[start_idx:end_idx]
    
    cat(sprintf("Processing chunk #%d: outcomes [%d..%d] (size %d)\n",
                chunk_counter, start_idx, end_idx, length(chunk_outcomes)))
    
    # Show progress bar for this chunk
    with_progress({
      p <- progressor(steps = length(chunk_outcomes))
      
      # 
      chunk_results_list <- future_lapply(
        X   = chunk_outcomes,
        FUN = function(one_outcome) {
          run_coxph_for_one_outcome(
            outcome_prefix = one_outcome,
            exposures      = exposures,
            progressor     = p
          )
        },
        future.seed = TRUE  # <--- Important for reproducible, parallel-safe RNG
      )
      
      # Combine chunk results
      chunk_results <- rbindlist(chunk_results_list, fill=TRUE)
      all_chunks_results[[chunk_counter]] <- chunk_results
    })
    
    chunk_counter <- chunk_counter + 1
  }
  
  # Combine *all* chunks
  final_results <- rbindlist(all_chunks_results, fill = TRUE)
  
  
  ###############################################################################
  #             7) Save final results
  ###############################################################################
  cat("\nWriting final combined results...\n")
  
  # output_dir <- paste0(phenCategory[k],"_UKB450k_chip_prs")
  
  # dir.create(output_dir, showWarnings = FALSE)
  
  out_file <- paste0(outdir,"/",Phecode_Group,".UKB_Cox_Results.csv")
  
  fwrite(final_results, out_file)
  
  cat("All done! Results in ", out_file,"\n") 
############################### END #########################  
  
  
