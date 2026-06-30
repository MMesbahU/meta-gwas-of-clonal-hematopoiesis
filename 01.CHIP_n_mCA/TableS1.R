#---------------------------
#---------------
# Total N ####
#---------------
library(data.table)
library(dplyr)

# UKBB All ####
# ukb
UKB_anc <- fread("/Volumes/medpop_esp2/skoyama/passing/ukb_kgpprojection/v03/out/ukb.kgp_projected.tsv.gz", 
                 header=T)
ukb450k.chip_table1 <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/n650k/chip_call/ukb450k.CHIP_phewas.csv.gz")
ukb450k.var_table1 <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/n650k/chip_call/ukb450k.variants_phewas.csv.gz")
load("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/ukb500k_covariates.rda");rm(d_base,d_pcs)
ukb450k.chip_table1 <- merge(ukb450k.chip_table1, d_ukb_all, 
                             by.x="eid_7089", by.y = "eid"); rm(d_ukb_all)

ukb450k.chip_table1$eid <- as.character(ukb450k.chip_table1$eid_7089)

summary(ukb450k.chip_table1$AGE_assessment)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 37.00   50.00   58.00   56.55   63.00   73.00 
round(summary(ukb450k.chip_table1$AGE_assessment))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 37      50      58      57      63      73 
 # round(sd(ukb450k.chip_table1$AGE_assessment))
# 8
(table(ukb450k.chip_table1$Genetic_Sex, exclude = NULL))
# Female   Male   <NA> 
#   246161 207673    493
round(prop.table(table(ukb450k.chip_table1$Genetic_Sex, exclude = NULL))*100,1)
# Female   Male   <NA> 
#   54.2   45.7    0.1 
# non-eur N=187889
# (187889/881904) *100
# [1] 21.30493
# EUR N=694015
#
ukb450k.chip_table1_n454327 <- merge(ukb450k.chip_table1, UKB_anc, 
                                     by = "eid")# ; rm(d_ukb_all)
nrow(ukb450k.chip_table1_n454327)
# 453930
# 454327
length(unique(ukb450k.chip_table1_n454327$eid))
# 453928

table(ukb450k.chip_table1_n454327$knn[!duplicated(ukb450k.chip_table1_n454327$eid)], exclude=NULL)
# AFR    AMR    EAS    EUR    SAS 
# 8796   2244   2428 430564   9896
table(ukb450k.chip_table1_n454327$Genetic_Sex[!duplicated(ukb450k.chip_table1_n454327$eid)], exclude=NULL)
# Female   Male   <NA> 
#   246161 207673     94

#### CHIP data in UKB
# CHIP VAR
ukbbch_var <- fread("/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/ch_var_all_ukb200k_n_ukb250k.7Mar2023.csv.gz", header=T)
head(ukbbch_var)
sort(table(ukbbch_var$Gene), decreasing = T)

sort(table(ukbbch_var$Protein_Change[ukbbch_var$Gene%in% c("IDH1", "IDH2")]), decreasing = T)

# CHIP: /medpop/esp2/mesbah/tools/CHIP_metaAnalysis//04.SNP_heritability/QuickPRS
# GWAS samples
# 192960+242846
# 435806
# number of individuals used in analysis = 192960
# * case-control counts for each trait:
#   - 'hasCH': 11384 cases and 181576 controls
# - 'hasCHvaf05': 10597 cases and 182363 controls
# - 'hasCHvaf10': 6349 cases and 186611 controls
# - 'hasDNMT3A': 6880 cases and 181576 controls
# - 'hasTET2': 1987 cases and 181576 controls
# - 'hasASXL1': 941 cases and 181576 controls
# - 'hasSF': 310 cases and 181576 controls
# - 'hasDDR': 432 cases and 181576 controls
# * number of individuals used in analysis = 242846
# * case-control counts for each trait:
#   - 'hasCH': 17079 cases and 225767 controls
# - 'hasCHvaf05': 13851 cases and 228995 controls
# - 'hasCHvaf10': 6615 cases and 236231 controls
# - 'hasDNMT3A': 9043 cases and 225767 controls
# - 'hasTET2': 3501 cases and 225767 controls
# - 'hasASXL1': 1710 cases and 225767 controls
# - 'hasSF': 407 cases and 225767 controls
# - 'hasDDR': 1138 cases and 225767 controls
ukb200k_chip <- fread("/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/CH_phenoCovar.ukb200k_N193342.28cols.03_05_2024.tsv.gz", header=T)

ukb250k_chip <- fread("/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/CH_phenoCovar.ukb250k_N243350.28cols.03_05_2024.tsv.gz", header=T)

ukb450k_chip <- as.data.frame(rbind(ukb200k_chip, ukb250k_chip)); rm(ukb200k_chip, ukb250k_chip)

head(ukb450k_chip)  

table(ukb450k_chip$Batch, exclude = NULL)
# ukb200k ukb250k 
# 193342  243350
## UKB baseline info 
## UKB Baseline
d_base <- fread("/Volumes/medpop_esp2/mesbah/projects/gxe/UKB_baseline/UKB500k_baseline_info_recoded_participant.tsv.gz")
# "ever_smoked.p20160_i0",         
# "Smoking_status.p20116_i0",
table(d_base$p20116_i0,  d_base$p20160_i0, exclude = NULL)
## Smoking status
ukb450k_chip$smking_status <- factor(ifelse(ukb450k_chip$FID %in% d_base$eid[d_base$p20160_i0=="Yes"], 
                                            "ever_smoker", 
                                            ifelse(ukb450k_chip$FID %in% d_base$eid[d_base$p20160_i0=="No"],
                                                   "never_smoker", "unknown")), 
                                     levels = c("never_smoker", "unknown", "ever_smoker"))
table(ukb450k_chip$smking_status, exclude = NULL)

ukb450k_chip$GenANC <- factor(ukb450k_chip$knn, 
                              levels = c("EUR", "AFR", "AMR", "EAS", "SAS"))
table(ukb450k_chip$GenANC)

summary(glm(hasDNMT3A ~ Age_at_recruitment + 
              GenANC + Genetic_Sex + 
      smking_status + Batch + GenoBatch + sqrtAge_at_recruitment + 
      PC1 + PC2 + PC3 + PC4 + PC5, data = ukb450k_chip, family = "binomial"))

 #### GLM 
ukb450k_chip %>% 
  group_by(Genetic_Sex , hasDNMT3A) %>% 
  summarise(n=n(), 
            average=mean(Age_at_recruitment), 
            SD=sd(Age_at_recruitment), 
            Median=median(Age_at_recruitment))

# Define your flag variables
flag_vars <- c("hasCH", "hasCHvaf10", "hasDNMT3A", 
               "hasTET2", "hasASXL1", "hasSF", "hasDDR")
# Create an empty list to store results for each variable
results_list <- list()
# results_list_noAge2 <- list()
# Loop through each flag variable
for (var in flag_vars) {
  # Get the case/control count as a table
  cnt <- table(ukb450k_chip[[var]], useNA = "ifany")
  
  # Fit the GLM model for this outcome variable.
  # Note: We use get(var) to refer to the column by name.
  model <- glm(get(var) ~ Age_at_recruitment + 
                 GenANC + Genetic_Sex + 
                 smking_status + Batch + 
                 GenoBatch + sqrtAge_at_recruitment + 
                 PC1 + PC2 + PC3 + PC4 + PC5,  
               data=ukb450k_chip, family = "binomial")
  
  # Extract selected rows and columns of the model coefficients.
  # Here we take rows: 2,3, 4,5,6,7,9 and columns: 1 to 4.
  coef_mat <- summary(model)$coefficients[c(2:7,9), 1:4]
  
  # Store both the count table and the coefficient matrix in the results list,
  # using the variable name as the list element name.
  results_list[[var]] <- list("case_control_count" = cnt,
                              "glm_coefficients" = coef_mat)
  # results_list_noAge2[[var]] <- list("case_control_count" = cnt,
  #                             "glm_coefficients" = coef_mat)
  # Optionally, print the results to the console:
  cat("Table for", var, ":\n")
  print(cnt)
  cat("\nGLM for", var, ":\n")
  print(coef_mat)
  cat("\n----------------------------\n")
}

results_list

## Mapping #####
# Define the mapping from flag variable names to descriptive outcome labels
outcome_mapping <- c(
  hasCH    = "Overall CHIP",
  hasCHvaf10 = "Expanded CHIP",
  hasDNMT3A  = "DNMT3A",
  hasTET2    = "TET2",
  hasASXL1   = "ASXL1",
  hasSF      = "Splicing Factors",
  hasDDR     = "DNA damage"
)

# Define a mapping for cleaned exposure names
exposure_mapping <- c(
  "Age_at_recruitment" = "Age",
  "Genetic_SexMale"      = "Male Sex",
  "GenANCAFR"               = "AFR",
  "GenANCAMR"               = "AMR",
  "GenANCEAS"               = "EAS",
  "GenANCSAS"               = "SAS",
  "smking_statusever_smoker"  = "Ever Smoker"
)

# Define the flag variables in the desired order
flag_vars <- c("hasCH", "hasCHvaf10", "hasDNMT3A", 
               "hasTET2", "hasASXL1", "hasSF", "hasDDR")

# Create an empty list to store data frames for each outcome
final_list <- list()
results_list <- results_list

for(var in flag_vars) {
  # Extract the GLM coefficient matrix for this outcome and convert it to a data frame
  coef_mat <- results_list[[var]]$glm_coefficients
  df <- as.data.frame(coef_mat)
  
  # Rename columns to Beta, SE, Z, and P
  colnames(df) <- c("Beta", "SE", "Z", "P")
  
  # Create an Exposure column using the row names of the coefficient matrix and clean them
  df$Exposure <- rownames(df)
  df$Exposure <- sapply(df$Exposure, function(x) {
    if(x %in% names(exposure_mapping)) exposure_mapping[x] else x
  })
  
  # Add the Outcome column using the outcome mapping
  df$Outcome <- outcome_mapping[[var]]
  
  # Extract the case/control counts from the table.
  # We assume that the table names "0" and "1" correspond to Controls and Cases, respectively.
  cnt <- results_list[[var]]$case_control_count
  controls <- if("0" %in% names(cnt)) as.numeric(cnt[["0"]]) else NA
  cases    <- if("1" %in% names(cnt)) as.numeric(cnt[["1"]]) else NA
  
  # Add the counts to every row of this outcome's data frame
  df$Cases <- cases
  df$Controls <- controls
  df$N <- cases + controls
  # Select and reorder columns as desired
  df <- df[, c("Outcome", "Exposure", "Beta", "SE", "Z", "P", "N", "Cases", "Controls")]
  
  # Store this data frame in the list
  final_list[[var]] <- df
}

# Combine the data frames for all outcomes into one data frame
final_results <- do.call(rbind, final_list)
rownames(final_results) <- NULL
# Print the final data frame
print(final_results)
# Optionally, save the final results to a file (e.g., as an RDA file)
# save(final_results, file = "final_glm_results.rda")
write.csv(final_results, "~/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/NG_v2/Figs/UKBB_CHIP.risk_factors.results_flag_vars.csv", row.names = FALSE)

#-----------------------------------------
# AoU All ####
# Done online
# * number of individuals used in analysis = 225864
# * case-control counts for each trait:
#   - 'hasCH': 8784 cases and 217080 controls
# - 'hasCHvaf10': 6752 cases and 217080 controls
# - 'hasDTA': 6317 cases and 217080 controls
# - 'hasDDR': 596 cases and 217080 controls
# - 'hasSF': 407 cases and 217080 controls
# - 'hasDNMT3A': 4371 cases and 217080 controls
# - 'hasTET2': 1480 cases and 217080 controls
# - 'hasASXL1': 647 cases and 217080 controls
# - 'hasPPM1D': 399 cases and 217080 controls
# - 'hasTP53': 198 cases and 217080 controls
# - 'hasSF3B1': 202 cases and 217080 controls
#-----------------------------------------

# TOPMED All ####
# GWAS samples
TOPMed_GWAS_pheno.MultiANC.
topmed_gwas <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/pheno/TOPMed_GWAS_pheno.MultiANC.04_23_2024.tsv.gz", header=T)
str(topmed_gwas)

## Data Bick et al. 2020
topmed.chip.2019.table1 <- fread("~/Documents/Project/TopMed/TOPMED_CHIP_100k_7-1-19/annot.TOPMed_100k_CHIP_calls_7_1_19_fordistribution.csv")
table(topmed.chip.2019.table1$Ever_smoker_baseline, exclude = NULL)
# 0     1  <NA> 
#   16202 25926 44988 
qcd.topmed.chip.2019.table1 <- topmed.chip.2019.table1 %>% filter(sample.id %in% topmed_filtered_ch.v2$NWD_ID)
nrow(qcd.topmed.chip.2019.table1)
table(qcd.topmed.chip.2019.table1$Ever_smoker_baseline, exclude = NULL)
#     0     1  <NA> 
#   15111 25401 24455
table(qcd.topmed.chip.2019.table1$Ever_smoker_baseline, 
      qcd.topmed.chip.2019.table1$Current_smoker_baseline, 
      exclude = NULL)
# 0     1  <NA>
# 0    15111     0     0
# 1    15230 10081    90
# <NA>   651    70 23734
qcd.topmed.chip.2019.table1 <- qcd.topmed.chip.2019.table1 %>%
  mutate(
    smoking_status = case_when(
      (Ever_smoker_baseline %in% 1) | (Current_smoker_baseline %in% 1) ~ "ever_smoker",
      
      # never_smoker if Current_smoker_baseline == 0 AND
      #              Ever_smoker_baseline is either 0 OR NA
      (Current_smoker_baseline %in% 0) &
        (Ever_smoker_baseline %in% 0 | is.na(Ever_smoker_baseline))  ~ "never_smoker",
      
      # Everything else is unknown
      TRUE ~ "unknown"
    ),
    smoking_status = factor(smoking_status,
                            levels = c("never_smoker","unknown","ever_smoker")
    )
  )

table(qcd.topmed.chip.2019.table1$smoking_status, exclude = NULL)
# never_smoker      unknown  ever_smoker 
# 15762        23734        25471                                                      
#
topmed_gwas$smoking_status <- factor(ifelse(topmed_gwas$NWD_ID %in% qcd.topmed.chip.2019.table1$sample.id[qcd.topmed.chip.2019.table1$smoking_status=="never_smoker"],
                                     "never_smoker",
                                     ifelse(topmed_gwas$NWD_ID %in% qcd.topmed.chip.2019.table1$sample.id[qcd.topmed.chip.2019.table1$smoking_status=="ever_smoker"], 
                                            "ever_smoker", 
                                            "unknown")),
                                     levels = c("never_smoker",
                                                "unknown",
                                                "ever_smoker"))

table(topmed_gwas$smoking_status, exclude = NULL)
# never_smoker      unknown  ever_smoker 
# 15762        33389        25471 
topmed_gwas$GenANC <- factor(topmed_gwas$knn, 
                             levels = c("EUR", "AFR", "AMR", "EAS", "SAS"))
### TOPMed Study
# N>3000
n3k_studies <- c(names(table(topmed_gwas$STUDY.y)[(table(topmed_gwas$STUDY.y)>=3000)]))

topmed_gwas$STUDY <- ifelse( !(topmed_gwas$STUDY.y %in% n3k_studies), 
                             "not_3k",  
                             topmed_gwas$STUDY.y)


table(topmed_gwas$STUDY, exclude = NULL)
### Sequencing Centers
table(topmed_gwas$SEQ_CTR, exclude = NULL)
# Baylor    Broad Illumina Macrogen     NYGC       UW    WashU 
# 23984    31370      714      854     1726    11360     4614 
topmed_gwas$Sequencing_Center <- ifelse(!(topmed_gwas$SEQ_CTR %in% 
                                            c("Baylor", "Broad","UW") ), 
                                        "Others", topmed_gwas$SEQ_CTR)


table(topmed_gwas$Sequencing_Center, exclude = NULL)
# Baylor  Broad Others     UW 
# 23984  31370   7908  11360
## TOPMed Phase
table(topmed_gwas$PHASE, exclude = NULL)

topmed_gwas$TOPMed_Phase <- factor(ifelse(!(topmed_gwas$PHASE %in% c(1,2,3) ), 
                                   "Others", topmed_gwas$PHASE))

table(topmed_gwas$TOPMed_Phase, exclude = NULL)

# GLM
summary(glm(hasCHvaf02 ~ AgeAtBloodDraw + Sqrd_AgeAtBloodDraw + 
              Genetic_Sex + GenANC + smoking_status + 
              STUDY + Sequencing_Center + TOPMed_Phase + 
              PC1 + PC2 + PC3 + PC4 + PC5, 
            data=topmed_gwas, family = "binomial"))  
#### 
names(topmed_gwas)
# Define your flag variables
flag_vars <- c("hasCHvaf02", "hasCHvaf10", "hasDNMT3A", 
               "hasTET2", "hasASXL1", "hasSF", "hasDDR")
# Create an empty list to store results for each variable
results_list <- list()
# results_list_noAge2 <- list()
# Loop through each flag variable
for (var in flag_vars) {
  # Get the case/control count as a table
  cnt <- table(topmed_gwas[[var]], useNA = "ifany")
  # Fit the GLM model for this outcome variable.
  # Note: We use get(var) to refer to the column by name.
  model <- glm(get(var) ~ AgeAtBloodDraw + Genetic_Sex + 
                 GenANC + smoking_status + 
                 STUDY + Sequencing_Center + TOPMed_Phase + 
                 Sqrd_AgeAtBloodDraw +
                 PC1 + PC2 + PC3 + PC4 + PC5,
               data=topmed_gwas, family = "binomial")
  
  # Extract selected rows and columns of the model coefficients.
  # Here we take rows: 2,3, 4,5,6,7,9 and columns: 1 to 4.
  coef_mat <- summary(model)$coefficients[c(2:7,9), 1:4]
  
  # Store both the count table and the coefficient matrix in the results list,
  # using the variable name as the list element name.
  results_list[[var]] <- list("case_control_count" = cnt,
                              "glm_coefficients" = coef_mat)
    # Optionally, print the results to the console:
  cat("Table for", var, ":\n")
  print(cnt)
  cat("\nGLM for", var, ":\n")
  print(coef_mat)
  cat("\n----------------------------\n")
}
results_list

## Mapping #####
# Define the mapping from flag variable names to descriptive outcome labels
outcome_mapping <- c(
  hasCHvaf02    = "Overall CHIP",
  hasCHvaf10 = "Expanded CHIP",
  hasDNMT3A  = "DNMT3A",
  hasTET2    = "TET2",
  hasASXL1   = "ASXL1",
  hasSF      = "Splicing Factors",
  hasDDR     = "DNA damage"
)

# Define a mapping for cleaned exposure names
exposure_mapping <- c(
  "AgeAtBloodDraw" = "Age",
  "Genetic_SexMale"      = "Male Sex",
  "GenANCAFR"               = "AFR",
  "GenANCAMR"               = "AMR",
  "GenANCEAS"               = "EAS",
  "GenANCSAS"               = "SAS",
  "smoking_statusever_smoker"  = "Ever Smoker"
)

# Define the flag variables in the desired order
flag_vars <- c("hasCHvaf02", "hasCHvaf10", "hasDNMT3A", 
               "hasTET2", "hasASXL1", "hasSF", "hasDDR")

# Create an empty list to store data frames for each outcome
final_list <- list()
for(var in flag_vars) {
  
  # Extract the GLM coefficient matrix for this outcome and convert it to a data frame
  coef_mat <- results_list[[var]]$glm_coefficients
  df <- as.data.frame(coef_mat)
  
  # Rename columns to Beta, SE, Z, and P
  colnames(df) <- c("Beta", "SE", "Z", "P")
  
  # Create an Exposure column using the row names of the coefficient matrix and clean them
  df$Exposure <- rownames(df)
  df$Exposure <- sapply(df$Exposure, function(x) {
    if(x %in% names(exposure_mapping)) exposure_mapping[x] else x
  })
  
  # Add the Outcome column using the outcome mapping
  df$Outcome <- outcome_mapping[[var]]
  
  # Extract the case/control counts from the table.
  # We assume that the table names "0" and "1" correspond to Controls and Cases, respectively.
  cnt <- results_list[[var]]$case_control_count
  controls <- if("0" %in% names(cnt)) as.numeric(cnt[["0"]]) else NA
  cases    <- if("1" %in% names(cnt)) as.numeric(cnt[["1"]]) else NA
  
  # Add the counts to every row of this outcome's data frame
  df$Cases <- cases
  df$Controls <- controls
  df$N <- cases + controls
  # Select and reorder columns as desired
  df <- df[, c("Outcome", "Exposure", "Beta", "SE", "Z", "P", "N", "Cases", "Controls")]
  
  # Store this data frame in the list
  final_list[[var]] <- df
}

# Combine the data frames for all outcomes into one data frame
final_results <- do.call(rbind, final_list)
rownames(final_results) <- NULL
# Print the final data frame
print(final_results)

# Optionally, save the final results to a file (e.g., as an RDA file)
# save(final_results, file = "final_glm_results.rda")
write.csv(final_results, "~/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/NG_v2/Figs/TOPMed_CHIP.risk_factors.results_flag_vars.csv", 
          row.names = FALSE)

#---------------------------------------------------------------------
# MGBB All ####
# "/medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.Feb2024.tsv.gz"
#  --covarColList PC{1:10},Age_Genotyping,sqrAge_Genotyping \
# --catCovarList Sex,ancestry_pred,Batch_CHIP_call 
## CHIP: U2AF1 included
# 12 new cases

#-------------------- MGBB CHIP ####
MGBB_Phenos <- fread("~/Documents/Project/Partners_CHIP/Kavvya/Data_2022/MGBB_Phenos_2CODE_2022-06-21.csv")

mgbb_var <- fread("~/Documents/Project/CHIP_annotation/2022_CHIP_Call/MGB_40k/MGBB53k_CHIP_return/MGBB_53k_WES_CHIP_variants_20Mar2023.csv", header=T)
mgbb_chip <- fread("~/Documents/Project/CHIP_annotation/2022_CHIP_Call/MGB_40k/MGBB53k_CHIP_return/MGBB_53k_WES_CHIP_calls_20Mar2023.csv", header=T)

mgbb_chip_demo <- merge(mgbb_chip, MGBB_Phenos, by="Biobank_Subject_ID")
nrow(mgbb_chip_demo)

length(unique(mgbb_chip_demo$Biobank_Subject_ID))
# 53306
round(summary(mgbb_chip_demo$Age_Genotyping))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2      40      56      54      67     102 
round(sd(mgbb_chip_demo$Age_Genotyping))
# 17
table(mgbb_chip_demo$Sex)
# Female   Male 
# 29662  23644
round(prop.table(table(mgbb_chip_demo$Sex))*100,1)
# Female   Male 
# 55.6   44.4

### MGB KNN Skyoma
# load(file = "~/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/genetic_ancestry.mgbb_topmed_ukbb.sk2024.rda")
MGB_anc <- fread("/Volumes/medpop_esp2/projects/MGB_Biobank/genotype/53K_GSA/release/GSA_53K.hg38.1KGPCA.tsv.gz",
                 header=T)

table(MGB_anc$knn, exclude = NULL)
# AFR   AMR   EAS   EUR   SAS 
# 2892  3851  1064 44892   607 
MGB_anc$SampleID <- as.integer(stringr::str_split_fixed(string = MGB_anc$IID, pattern = "[-]", n=2)[,2])
length(unique(MGB_anc$SampleID))
# 53306
# table(mgbb_chip$Biobank_Subject_ID %in% MGB_anc$SampleID, exclude = NULL)
# # FALSE  TRUE 
# # 352 52993 
mgbb53k.chip_predANC <- merge(mgbb_chip_demo, MGB_anc[,c(14,3)], 
                              by.x="Biobank_Subject_ID",
                              by.y="SampleID")


# summary for Suppl table 1
round(summary(mgbb53k.chip_predANC$Age_Genotyping))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2      40      56      54      67     102 
round(sd(mgbb53k.chip_predANC$Age_Genotyping))
# 17
table(mgbb53k.chip_predANC$knn, exclude=NULL)
# AFR   AMR   EAS   EUR   SAS
# 2869  3802  1059 44653   601 
round(prop.table(table(mgbb53k.chip_predANC$knn, exclude=NULL))*100,1)
# AFR  AMR  EAS  EUR  SAS 
# 5.4  7.2  2.0 84.3  1.1
table(mgbb53k.chip_predANC$Sex, exclude=NULL)
# Female   Male 
# 29477  23507
round(prop.table(table(mgbb53k.chip_predANC$Sex, exclude=NULL)) *100,1)
# Female   Male 
# 55.6   44.4
##
### MGBB Smoking ####
qcd.mgbb_survey_dat <- fread("~/Documents/Project/Partners_CHIP/Health Information Survey/20220510_Survey_Data_1706_clean_ver_1_0.csv", header=T) 
length(unique(qcd.mgbb_survey_dat$subject_id))
# 99193
table(qcd.mgbb_survey_dat$subject_id %in% mgbb_chip_demo$Biobank_Subject_ID)
table(qcd.mgbb_survey_dat$smoke, exclude = NULL)
#   1     2     3  <NA> 
#   36843 21886  3072 37392
table(qcd.mgbb_survey_dat$smoke[qcd.mgbb_survey_dat$subject_id %in% mgbb_chip_demo$Biobank_Subject_ID], exclude = NULL)
# 1     2     3  <NA> 
#   15304  9240  1437 15172 
# "Have you smoked at least 100 cigarettes in your lifetime?"
# [1] ""No""
# [2] ""Yes, smoked in past, but quit""
# [3] ""Yes, currently smoke"""
mgbb53k.chip_predANC$smk_100cig_ever <- factor(ifelse(mgbb53k.chip_predANC$Biobank_Subject_ID %in%qcd.mgbb_survey_dat$subject_id[!is.na(qcd.mgbb_survey_dat$smoke) & qcd.mgbb_survey_dat$smoke==1],
                                               "No", 
                                               ifelse(mgbb53k.chip_predANC$Biobank_Subject_ID %in%qcd.mgbb_survey_dat$subject_id[!is.na(qcd.mgbb_survey_dat$smoke) & (qcd.mgbb_survey_dat$smoke==2 | qcd.mgbb_survey_dat$smoke==3)], 
                                                      "Yes", "Unknown")),
                                               levels = c("No","Unknown", "Yes"))
table(mgbb53k.chip_predANC$smk_100cig_ever, exclude=NULL)
# No Unknown     Yes 
# 15208   27150   10626
# GLM
summary(mgbb53k.chip_predANC %>% glm(hasCHIP ~ Age_Genotyping + Sex + knn + 
                                       smk_100cig_ever + (Age_Genotyping)^2 + factor(Batch_CHIP_call), 
                                     data = ., family = "binomial"))

table(mgbb53k.chip_predANC$smk_100cig_ever[mgbb53k.chip_predANC$Biobank_Subject_ID %in% mgbb_gwas_samp_ch$Biobank_Subject_ID])
# smoke_age_1	"During which age ranges did you smoke? (please check all ages that apply) (choice=Less than 15 years old)"                                                                                                                                                                                                                                                                     
# <15 years old
# 15-19
# 20-29
# 30-39
# 40-49
# 50-59
# 60-69
# 70-79
# 80-89
# 90 and over
# smoke_no_1	"Number of cigarettes smoked per day when less than 15 years old"                                                                                                                                                                                                                                                                                                               	"[1] ""None""
# [2] "1-4"
# [3] "5-14"
# [4] "15-24"
# [4] "25-35"
# [5] "36-44"
# [6] "45+"

# Physical data
phy <- fread("/Volumes/mesbah/dataset/mgbb/mgbb_phenos/Phy.txt.gz", header = T, sep="|")
sort(table(phy$Concept_Name))

smking_concepts <- unique(phy$Concept_Name)[grepl(pattern = "Smoking Tobacco Use",x = unique(phy$Concept_Name), ignore.case = T)]

smk_phy <- phy %>% filter(Concept_Name %in% smking_concepts); rm(phy)


smking_concepts_code <- smk_phy[,c(4,6)] %>% filter(Concept_Name %in% smking_concepts & !duplicated(Code) )
never_smoker <- c("SH-SMKTBCO-5", "SH-SMKTBCO-7")

ever_smoker <- c("SH-SMKTBCO-1", "SH-SMKTBCO-2", 
                 "SH-SMKTBCO-4", "SH-SMKTBCO-9", 
                 "SH-SMKTBCO-10")
not_known <- c("SH-SMKTBCO-8", "SH-SMKTBCO-6", "SH-SMKTBCO-3")

table(unique(smk_phy$Subject_Id[smk_phy$Code %in% never_smoker]) %in% 
        unique(smk_phy$Subject_Id[smk_phy$Code %in% ever_smoker]) )
# FALSE  TRUE 
# 33223  5103
all_mgbb_samples <- unique(smk_phy$Subject_Id)

table(unique(smk_phy$Subject_Id[smk_phy$Code %in% never_smoker]) %in% 
        unique(smk_phy$Subject_Id[smk_phy$Code %in% not_known]) )
# FALSE  TRUE 
# 36783  1543
head(sort((table(smk_phy$Subject_Id[(smk_phy$Subject_Id[smk_phy$Code %in% never_smoker] %in% smk_phy$Subject_Id[smk_phy$Code %in% ever_smoker])])), decreasing = T))
smk_phy[smk_phy$Subject_Id=="10052116",]
# Convert the date column to Date class (assuming "m/d/yyyy" format)
smk_phy <- smk_phy %>%
  mutate(DATE = as.Date(Date, format = "%m/%d/%Y")) %>%
  arrange(DATE)

# # Convert character to Date
# smk_phy$DATE <- as.Date(smk_phy$Date, format = "%m/%d/%Y")
# 
# # Sort by Date and store in a new object (or overwrite)
# smk_phy_sorted <- smk_phy[order(smk_phy$DATE), ]


### get list of samples
list_never_smoker <- unique(smk_phy$Subject_Id[smk_phy$Code %in% never_smoker])
list_ever_smoker <- unique(smk_phy$Subject_Id[smk_phy$Code %in% ever_smoker])
list_notKnown <- unique(smk_phy$Subject_Id[smk_phy$Code %in% not_known])
#
mgbb53k.chip_predANC$Ever_Smoker <- factor(ifelse(mgbb53k.chip_predANC$Biobank_Subject_ID %in% list_never_smoker & 
                                             !(mgbb53k.chip_predANC$Biobank_Subject_ID %in% list_ever_smoker), 
                                           "Never Smoker", 
                                           ifelse(mgbb53k.chip_predANC$Biobank_Subject_ID %in% list_ever_smoker, 
                                                  "Ever Smoker", "Not Known")), 
                                           levels = c("Never Smoker","Not Known", "Ever Smoker"))
table(mgbb53k.chip_predANC$Ever_Smoker, exclude = NULL)
# fwrite(mgbb53k.chip_predANC, "~/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MGBB_Dem_smoking_knn.Feb2025.csv.gz", row.names = F, col.names = T)

############
#------------------
##
# dem <- fread("/Volumes/mesbah/dataset/mgbb/mgbb_phenos/Dem.txt.gz", header = T, sep="|")

summary(mgbb53k.chip_predANC %>% glm(hasCHIP ~ Age_Genotyping + Sex + knn + 
                                       smk_100cig_ever + (Age_Genotyping)^2 + factor(Batch_CHIP_call), 
                                     data = ., family = "binomial"))

summary(mgbb53k.chip_predANC %>% glm(hasCHIP ~ Age_Genotyping + Sex + knn + 
                                       Ever_Smoker + (Age_Genotyping)^2 + factor(Batch_CHIP_call), 
                                     data = ., family = "binomial"))

names(mgbb53k.chip_predANC)
head(mgbb53k.chip_predANC[,c(1:16,87,88,89)])

# 
# MGBB GWAS samples ####
# mgbb_gwas_samp_ch <- fread("/Users/muddin/Documents/Project/CHIP_annotation/2022_CHIP_Call/MGB_40k/MGBB53k_CHIP_return/GWAS_pheno/mgbb53k.imp_new.noRel_sk.Feb2024.tsv.gz", header=TRUE)
mgbb_gwas_samp_ch <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.Feb2024.tsv.gz", header=TRUE)

names(mgbb_gwas_samp_ch)

mgbb_gwas_samp_ch <- merge(mgbb_gwas_samp_ch, 
                           mgbb53k.chip_predANC[,c(1:16,87,88,89)], 
                           by="Biobank_Subject_ID")

# fwrite(mgbb_gwas_samp_ch, "~/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/mgbb_gwas_samples.ch_risk.csv.gz", row.names = F, col.names = TRUE)
names(mgbb_gwas_samp_ch)
# Risk factors for prevalent CHIP 
table(mgbb_gwas_samp_ch$ancestry_pred, mgbb_gwas_samp_ch$knn)

mgbb_gwas_samp_ch$ancestry_pred <- factor(mgbb_gwas_samp_ch$ancestry_pred, 
                                          levels = c("EUR", "AFR", "AMR", "EAS", "SAS"))

table(mgbb_gwas_samp_ch$ancestry_pred)
str(mgbb_gwas_samp_ch$ancestry_pred)

flag_vars <- c("hasCHvaf02", "hasCHvaf10", "hasDNMT3A", 
               "hasTET2", "hasASXL1", "hasSF", "hasDDR")

# Loop through each variable and print its table
for (var in flag_vars) {
  cat("Table for", var, ":\n")
  print( table(mgbb_gwas_samp_ch[[var]], useNA = "ifany"))
  cat("\n")
  
  cat("GLM for", var, ":\n")
  print( 
    summary(glm(get(var) ~ Age_Genotyping.x + 
                  Sex.x + ancestry_pred +
                Ever_Smoker + sqrAge_Genotyping +
                  PC1 + PC2 + PC3 + PC4+ PC5+
                  Batch_CHIP_call.x, 
                data=mgbb_gwas_samp_ch, family="binomial"))
  )
  cat("\n")
}

##
for (var in flag_vars) {
  cat("Table for", var, ":\n")
  print( table(mgbb_gwas_samp_ch[[var]], useNA = "ifany"))
  cat("\n")
  
  cat("GLM for", var, ":\n")
  print( 
    summary(glm(get(var) ~ Age_Genotyping.x + 
                  Sex.x + ancestry_pred + 
                  Ever_Smoker + # sqrAge_Genotyping +
                  PC1 + PC2 + PC3 + PC4+ PC5+
                  Batch_CHIP_call.x, 
                data=mgbb_gwas_samp_ch, family="binomial"))$coefficients[c(2:7,9),1:4]
  )
  cat("\n")
}

##--------------------------------
# Define your flag variables
flag_vars <- c("hasCHvaf02", "hasCHvaf10", "hasDNMT3A", 
                           "hasTET2", "hasASXL1", "hasSF", "hasDDR")


# Create an empty list to store results for each variable
results_list <- list()
# results_list_noAge2 <- list()
# Loop through each flag variable
for (var in flag_vars) {
  
  # Get the case/control count as a table
  cnt <- table(mgbb_gwas_samp_ch[[var]], useNA = "ifany")
  
  # Fit the GLM model for this outcome variable.
  # Note: We use get(var) to refer to the column by name.
  model <- glm(get(var) ~ Age_Genotyping.x + 
                 Sex.x + ancestry_pred + 
                 Ever_Smoker +  sqrAge_Genotyping +
                 PC1 + PC2 + PC3 + PC4+ PC5+
                 Batch_CHIP_call.x, 
               data=mgbb_gwas_samp_ch, family = "binomial")
  
  # Extract selected rows and columns of the model coefficients.
  # Here we take rows: 2,3, 4,5,6,7,9 and columns: 1 to 4.
  coef_mat <- summary(model)$coefficients[c(2:7,9), 1:4]
  
  # Store both the count table and the coefficient matrix in the results list,
  # using the variable name as the list element name.
  results_list[[var]] <- list("case_control_count" = cnt,
                              "glm_coefficients" = coef_mat)
  # results_list_noAge2[[var]] <- list("case_control_count" = cnt,
  #                             "glm_coefficients" = coef_mat)
  # Optionally, print the results to the console:
  cat("Table for", var, ":\n")
  print(cnt)
  cat("\nGLM for", var, ":\n")
  print(coef_mat)
  cat("\n----------------------------\n")
}

results_list

# Save the results list to an RDA file.
# This file will include the table counts (with row names preserved) and
# the GLM coefficient matrices.
# save(results_list, file = "~/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/NG_v2/Figs/MGBB_CHIP.risk_factors.age2adj.results_flag_vars.rda")
# save(results_list_noAge2, file = "~/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/NG_v2/Figs/MGBB_CHIP.risk_factors.noAge2adj.results_flag_vars.rda")
load("~/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/NG_v2/Figs/MGBB_CHIP.risk_factors.age2adj.results_flag_vars.rda")
load("~/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/NG_v2/Figs/MGBB_CHIP.risk_factors.noAge2adj.results_flag_vars.rda")

## Mapping #####
# Define the mapping from flag variable names to descriptive outcome labels
outcome_mapping <- c(
  hasCHvaf02    = "Overall CHIP",
  hasCHvaf10 = "Expanded CHIP",
  hasDNMT3A  = "DNMT3A",
  hasTET2    = "TET2",
  hasASXL1   = "ASXL1",
  hasSF      = "Splicing Factors",
  hasDDR     = "DNA damage"
)

# Define a mapping for cleaned exposure names
exposure_mapping <- c(
  "Age_Genotyping.x" = "Age",
  "Sex.xMale"      = "Male Sex",
  "ancestry_predAFR"               = "AFR",
  "ancestry_predAMR"               = "AMR",
  "ancestry_predEAS"               = "EAS",
  "ancestry_predSAS"               = "SAS",
  "Ever_SmokerEver Smoker"  = "Ever Smoker"
)

# Define the flag variables in the desired order
flag_vars <- c("hasCHvaf02", "hasCHvaf10", "hasDNMT3A", 
               "hasTET2", "hasASXL1", "hasSF", "hasDDR")

# Create an empty list to store data frames for each outcome
final_list <- list()
results_list <- results_list
# results_list <- results_list_noAge2

for(var in flag_vars) {
  
  # Extract the GLM coefficient matrix for this outcome and convert it to a data frame
  coef_mat <- results_list[[var]]$glm_coefficients
  df <- as.data.frame(coef_mat)
  
  # Rename columns to Beta, SE, Z, and P
  colnames(df) <- c("Beta", "SE", "Z", "P")
  
  # Create an Exposure column using the row names of the coefficient matrix and clean them
  df$Exposure <- rownames(df)
  df$Exposure <- sapply(df$Exposure, function(x) {
    if(x %in% names(exposure_mapping)) exposure_mapping[x] else x
  })
  
  # Add the Outcome column using the outcome mapping
  df$Outcome <- outcome_mapping[[var]]
  
  # Extract the case/control counts from the table.
  # We assume that the table names "0" and "1" correspond to Controls and Cases, respectively.
  cnt <- results_list[[var]]$case_control_count
  controls <- if("0" %in% names(cnt)) as.numeric(cnt[["0"]]) else NA
  cases    <- if("1" %in% names(cnt)) as.numeric(cnt[["1"]]) else NA
  
  # Add the counts to every row of this outcome's data frame
  df$Cases <- cases
  df$Controls <- controls
  df$N <- cases + controls
  # Select and reorder columns as desired
  df <- df[, c("Outcome", "Exposure", "Beta", "SE", "Z", "P", "N", "Cases", "Controls")]
  
  # Store this data frame in the list
  final_list[[var]] <- df
}

# Combine the data frames for all outcomes into one data frame
final_results <- do.call(rbind, final_list)
rownames(final_results) <- NULL
final_results$Model <- "withAge2"
final_results_v1 <- final_results

final_results$Model <- "withoutAge2"
final_results_v2 <- final_results
# Print the final data frame
print(final_results_v1)

# comb.final_mgbb <- rbind(final_results_v1, final_results_v2)

# Optionally, save the final results to a file (e.g., as an RDA file)
# save(final_results, file = "final_glm_results.rda")
write.csv(final_results_v1, "~/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/NG_v2/Figs/MGBB_CHIP.risk_factors.age2adj.results_flag_vars.csv", row.names = FALSE)
# write.csv(final_results_v2, "~/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/NG_v2/Figs/MGBB_CHIP.risk_factors.noAge2adj.results_flag_vars.csv", row.names = FALSE)

##

# tp53mgbb_var <- mgbb_var %>%filter(Gene=="TP53")
# summary(tp53mgbb_var)
# tp53mgbb_chip <- mgbb_chip %>%filter(grepl(pattern = "TP53",x = Gene, ignore.case = T))
# summary(tp53mgbb_chip)
# 
# tp53_summary <- mgbb_var %>% 
#   filter(Gene=="TP53") %>% 
#   summarize(
#     nCHIP=n(),
#     N_samp=length(unique(Biobank_Subject_ID)),
#     avgDP=mean(DP),
#     minDP=min(DP),
#     maxDP=max(DP),
#     avgAD=mean(AD_Alt),
#     minAD=min(AD_Alt),
#     maxAD=max(AD_Alt),
#     avgADF=mean(ADF_Alt),
#     minADF=min(ADF_Alt),
#     maxADF=max(ADF_Alt),
#     avgADR=mean(ADR_Alt),
#     minADR=min(ADR_Alt),
#     maxADR=max(ADR_Alt),
#     avgVAF=mean(VAF),
#     minVAF=min(VAF),
#     maxVAF=max(VAF)
# )
#
#----------------------------


#-------------------------------
#---------------
# GWAS N ####
#---------------
# UKBB GWAS ####
# load(file = "~/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/genetic_ancestry.mgbb_topmed_ukbb.sk2024.rda")

#-------------------------------