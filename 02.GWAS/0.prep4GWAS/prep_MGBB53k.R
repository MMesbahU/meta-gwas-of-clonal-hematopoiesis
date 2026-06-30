
############### MGBB CHIP and Covariates 

#### Kaavya pheno file updated on Jun 21 2022 
MGBB_Phenos_2CODE_2022 <- fread("~/Documents/Project/Partners_CHIP/Kavvya/Data_2022/MGBB_Phenos_2CODE_2022-06-21.csv")


##### All WES samples data where CHIP calling is done
# MGB_53445WES <- fread("/Volumes/mesbah/datasets/CHIP/MGBB/MGB_wes_N53445.tsv")
# names(MGB_53445WES) <- "Biobank_Subject_ID"            

#### MGBB Demographic info from Buu
MGBB_Basic_covars <- fread("/Volumes/medpop_esp2/btruong/Projects/sharing/MGBB_basic_covariates.txt", header = T)

# 
mgbb53k.pheno <- merge(MGBB_Phenos_2CODE_2022, MGBB_Basic_covars, 
                         by.x = "Biobank_Subject_ID", by.y =  "IID")
# cor(mgbb53k.pheno$Age_Genotyping, mgbb53k.pheno$age_enroll)

######## MGBB Self Reported Ancestry catagory 
# Coding from Jacqueline
asian <- c("ASIAN","ASIAN@BLACK OR AFRICAN AMERICAN","ASIAN@NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER",
           "ASIAN@OTHER","ASIAN@UNAVAILABLE","ASIAN@WHITE","NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER",
           "NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER@WHITE")

black <- c("BLACK","BLACK OR AFRICAN AMERICAN","BLACK OR AFRICAN AMERICAN@OTHER",
           "BLACK OR AFRICAN AMERICAN@UNAVAILABLE","BLACK OR AFRICAN AMERICAN@WHITE")

mgbb53k.pheno$Ancestry_Self_cat <- ifelse(mgbb53k.pheno$ancestry_self_report %in% asian, "ASIAN",
                                            ifelse(mgbb53k.pheno$ancestry_self_report %in% black, "BLACK",
                                                   ifelse(mgbb53k.pheno$ancestry_self_report=="WHITE", "WHITE", "OTHER")))
prop.table(table(mgbb53k.pheno$Ancestry_Self_cat))

names(mgbb53k.pheno)
############################

## Merge mgbb53k.pheno; n= 53,195
# 53345 with CHIP status (excluding 100 samples used in Panel-of-Normal)
mgbb53k <- merge(mgbb53k.pheno[,c(1:4,6,28:30,82,84,86:97)], 
                 mgbb53k.chip,
                 by="Biobank_Subject_ID")
mgbb53k$sqrAge_Genotyping <- mgbb53k$Age_Genotyping^2
## read genotyped sample file
sam.bgen1 <- fread("/Volumes/medpop_esp2/projects/MGB_Biobank/imputation/53K_GSA/release/bgen/GSA_53K.merged.chr1.sample")

# Biobank ID
sam.bgen1$Biobank_Subject_ID <- as.numeric(stringr::str_split_fixed(string = sam.bgen1$ID_2, pattern = "[-]", n = 2)[,2])

str(sam.bgen1)

mgbb53k.imp <- merge(sam.bgen1[,c(1,2,5)], mgbb53k, by = "Biobank_Subject_ID") 
names(mgbb53k.imp)

mgbb53k.imp <- mgbb53k.imp[, c(2,3,1,4:37)]

names(mgbb53k.imp) <- c("FID", "IID", names(mgbb53k.imp)[3:37])
#############

## 53k updated PCs
pcs.mgb53k <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/PC1_20.qced.mgbb53k.autosomes.eigenvec")

## Updated mgb gwas data
mgbb53k.imp_new <- merge(mgbb53k.imp[,c(2,1,3:13,24:37)], 
                         pcs.mgb53k[,c(2:22)], 
                         by="IID")

mgbb53k.imp_new <- mgbb53k.imp_new[,c(2,1,3:47)]
  # Read MGB CHIP mutation file
mgbb53k.var_table1 <- fread("~/Documents/Project/CHIP_annotation/2022_CHIP_Call/MGB_40k/MGBB_52966_chip_var.b1_2_3.5Aug2022.csv")

########## save files
# N = 52966
# fwrite(mgbb53k.imp_new, "/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.chip_call_all52966.tsv.gz",
#        row.names = F, col.names = T, sep="\t", na = 'NA', quote = F)

# fwrite(mgbb53k.var_table1, "/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.chip_var_N5910_nchip7145.csv.gz",
#        row.names = F, col.names = T, sep=",", na = 'NA')
######### 
## Exclude related samples
## SK rel filter
# mgb.rel_pair <- fread("/Volumes/medpop_esp2/projects/MGB_Biobank/genotype/53K_GSA/release/GSA_53K.related.pair.tsv.gz")
mgb.rel <- fread("/Volumes/medpop_esp2/projects/MGB_Biobank/genotype/53K_GSA/release/GSA_53K.related.tsv", 
                 header = F)

mgbb53k.imp_new.noRel <- subset(mgbb53k.imp_new, !(mgbb53k.imp_new$IID %in% mgb.rel$V1) )

  # multi-ancestry
fwrite(mgbb53k.imp_new.noRel, "/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.21Jul2022.tsv.gz", 
       row.names = F, col.names = T, sep = "\t", 
       quote = F, na = "NA")
  # Europen only
fwrite(mgbb53k.imp_new.noRel[mgbb53k.imp_new.noRel$Ancestry_Self_cat=="WHITE", ], "/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k_WHITE.imp_new.noRel_sk.21Jul2022.tsv.gz", 
       row.names = F, col.names = T, sep = "\t", 
       quote = F, na = "NA")

### 
prop.table(table(mgbb53k.imp_new.noRel$Ancestry_Self_cat, exclude = NULL))*100
# ASIAN     BLACK     OTHER     WHITE 
# 3.065302  5.057151  7.736805 84.140741 
prop.table(table(mgbb53k.imp_new.noRel$Sex, exclude = NULL))*100
# Female     Male 
# 55.36428 44.63572 

### Save all files
# save.image(file = "~/Documents/Project/CHIP_GWAS/rerun/ukb450.mgb53k.topmed127k.rda")

