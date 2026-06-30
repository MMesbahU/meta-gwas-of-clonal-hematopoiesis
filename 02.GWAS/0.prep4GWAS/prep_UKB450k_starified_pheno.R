

#################### Update Pheno with Genetic Ancestry #####
### Load rda file created in script "prep_pheno4GWAS.ukb450k.ch_n_vaf.R"

load("/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/CH_phenoCovar.ukb200k_ukb250k.8Mar2023.rda")
####
library(data.table)

library("glue")

library("tidyverse")

date_label = format(Sys.Date(), "%m_%d_%Y")

# SK's UKB genetic Ancestry
UKB_anc <- fread("/Volumes/medpop_esp2/skoyama/passing/ukb_kgpprojection/v03/out/ukb.kgp_projected.tsv.gz", 
                 header=T)

UKB_anc$FID <- as.numeric(UKB_anc$eid)
head(UKB_anc)
summary(UKB_anc)
str(UKB_anc)
names(UKB_anc)
## UKB200k
d_ukb200k.nomiss_genAnc <- merge (d_ukb200k.nomiss[, c(1,2,31:41,15,16,17,19:29)], 
                                  UKB_anc[, c(19,12)], 
                                  by="FID")
str(d_ukb200k.nomiss_genAnc)
## UKB250k
d_ukb250k.nomiss_genAnc <- merge (d_ukb250k.nomiss[,c(1,2,31:41,15,16,17,19:29)], 
                                  UKB_anc[!duplicated(UKB_anc$eid), c(19,12)], 
                                  by="FID")
str(d_ukb250k.nomiss_genAnc)



######################## Plot
d_ukb250k.nomiss_genAnc %>%
  group_by(knn) %>%
  #filter(n() >= 1e3) %>%
  ggplot(data = ., aes(x = Age_at_recruitment, y = hasCH, colour = knn)) +
  # geom_point(alpha = .4) + 
  geom_smooth(aes(y=hasCH, colour=knn), 
              method ="glm", 
              method.args = list(family = "binomial")
              , se = T) +
  ggtitle("CHIP prevalence by Ancestry") +
  theme_bw() +
  scale_y_continuous(labels = scales::percent) +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13)
  ) + 
  labs(x = "age at blood draw", y = "CHIP prevalence")


######################
### Save file
outDir <- "/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/"

### Combined ukb 200k + 250k 
d_ukb200k_n_250k.nomiss_genAnc <- as.data.frame(rbind(d_ukb200k.nomiss_genAnc, 
                                                      d_ukb250k.nomiss_genAnc))
nrow(d_ukb200k_n_250k.nomiss_genAnc)
head(d_ukb200k_n_250k.nomiss_genAnc)
table(d_ukb200k_n_250k.nomiss_genAnc$knn)
  #
cohortNam <- "ukb200k_n_250k"

col_nums <- ncol( get( paste0("d_",cohortNam,".nomiss_genAnc")) )

N_pop=nrow( get( paste0("d_",cohortNam,".nomiss_genAnc")))

fwrite( get( paste0("d_",cohortNam,".nomiss_genAnc")), 
       glue("{outDir}CH_phenoCovar.{COHORT}_N{Sample_Size}.{NCOLs}cols.{date}.tsv.gz", 
            date = date_label, 
            Sample_Size=N_pop, 
            NCOLs=col_nums,
            OUTDir=outDir,
            COHORT=cohortNam),
       row.names = F, col.names = T, sep="\t", na = 'NA', quote = F)

###############
## UKB 200k
col_nums <- ncol(d_ukb200k.nomiss_genAnc)
cohortNam <- "ukb200k"
N_pop=nrow(d_ukb200k.nomiss_genAnc)
fwrite(d_ukb200k.nomiss_genAnc, 
       glue("{outDir}CH_phenoCovar.{COHORT}_N{Sample_Size}.{NCOLs}cols.{date}.tsv.gz", 
            date = date_label, 
            Sample_Size=N_pop, 
            NCOLs=col_nums,
            OUTDir=outDir,
            COHORT=cohortNam),
       row.names = F, col.names = T, sep="\t", na = 'NA', quote = F)

## UKB 250k
col_nums <- ncol(d_ukb250k.nomiss_genAnc)
cohortNam <- "ukb250k"
N_pop=nrow(d_ukb250k.nomiss_genAnc)

fwrite(d_ukb250k.nomiss_genAnc, 
       glue("{outDir}CH_phenoCovar.{COHORT}_N{Sample_Size}.{NCOLs}cols.{date}.tsv.gz", 
            date = date_label, 
            Sample_Size=N_pop, 
            NCOLs=col_nums,
            OUTDir=outDir,
            COHORT=cohortNam),
       row.names = F, col.names = T, sep="\t", na = 'NA', quote = F)


#####################
## stratified by Sex
## stratified by Sex
library(data.table)
library("glue")
library("tidyverse")
date_label = format(Sys.Date(), "%m_%d_%Y")
## ukb 200k
outDir <- "/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/"
## UKB 200k or 250k 
# or "ukb200k_n_250k"
cohortNam <- "ukb200k_n_250k"

SEX=c("Male", "Female")

for(i in 1:length(SEX)){
  
  cat(cohortNam, ":", SEX[i], "\n")
  
  SEX_lab <- SEX[i]
  
  col_nums <- ncol(get(paste0("d_",cohortNam,".nomiss_genAnc")) %>% 
                     filter(Genetic_Sex==SEX[i]))
  
  N_pop <- nrow(get(paste0("d_",cohortNam,".nomiss_genAnc")) %>% 
                  filter(Genetic_Sex==SEX[i]))
  
  fwrite(get(paste0("d_",cohortNam,".nomiss_genAnc")) %>% 
           filter(Genetic_Sex==SEX[i]), 
         glue("{outDir}CH_phenoCovar.{COHORT}_{sex_lab}_N{Sample_Size}.{NCOLs}cols.{date}.tsv.gz", 
              date = date_label, 
              Sample_Size=N_pop, 
              NCOLs=col_nums,
              OUTDir=outDir,
              COHORT=cohortNam,
              sex_lab=SEX_lab),
         row.names = F, col.names = T, sep="\t", na = 'NA', quote = F)
  
}

################# By Genetic Ancestry

## 
library(data.table)
library("glue")
library("tidyverse")
date_label = format(Sys.Date(), "%m_%d_%Y")
## ukb 200k + 250k
outDir <- "/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/"
## 
# or "ukb200k_n_250k"
cohortNam <- "ukb200k_n_250k"
# cohortNam <- "ukb250k"
ANCESTRY <- unique(d_ukb200k_n_250k.nomiss_genAnc$knn)

for(i in 1:length(ANCESTRY)){
  
  cat(cohortNam, ":", ANCESTRY[i], "\n")
  
  Anc_lab <- ANCESTRY[i]
  
  col_nums <- ncol(get(paste0("d_",cohortNam,".nomiss_genAnc")) %>% 
                     filter(knn==ANCESTRY[i]))
  
  N_pop <- nrow(get(paste0("d_",cohortNam,".nomiss_genAnc")) %>% 
                  filter(knn==ANCESTRY[i]))
  
  fwrite(get(paste0("d_",cohortNam,".nomiss_genAnc")) %>% 
           filter(knn==ANCESTRY[i]), 
         glue("{outDir}CH_phenoCovar.{COHORT}_{anc_lab}_N{Sample_Size}.{NCOLs}cols.{date}.tsv.gz", 
              date = date_label, 
              Sample_Size=N_pop, 
              NCOLs=col_nums,
              OUTDir=outDir,
              COHORT=cohortNam,
              anc_lab=Anc_lab),
         row.names = F, col.names = T, sep="\t", na = 'NA', quote = F)
  
}


