#########
library(data.table)
library(dplyr)

##
setwd("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/prs/")
#########

### CHIP PRS #####

## PRS with AFR Reference Panel
list_prs_score.afr <- system("ls /Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/prs/scores/afr/*.profile | awk '{print $NF}'", 
                             intern = TRUE)

names_prs_score.afr <- system("ls /Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/prs/scores/afr/*.profile | awk '{print $NF}' | cut -d '.' -f7,8  | sed 's:has::g' | tr '.' '_' ", 
                              intern = TRUE)
## Read PRS files
d_afr <- fread(list_prs_score.afr[1], header=T, select = c("ID2", "Profile_1"))

str(d_afr)

names(d_afr)

summary(d_afr)

names(d_afr) <- c("ID2", names_prs_score.afr[1])

str(d_afr)

dat_afr <- d_afr; rm(d_afr)

for(i in 2:length(list_prs_score.afr)){
  print(i)
  
  print(list_prs_score.afr[i])
  
  print(names_prs_score.afr[i])
  
  d <- fread(list_prs_score.afr[i], header=T, select = c("ID2", "Profile_1"))
  
  names(d) <- c("ID2", names_prs_score.afr[i])
  
  dat_afr <- merge(dat_afr, d, by="ID2")
  
  rm(d)
  
}

str(dat_afr)

head(dat_afr)

###
fwrite(dat_afr, "/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/prs/prs.afr.ukb450k.csv.gz", 
       row.names = F, col.names = T, sep=",")

save(dat_afr, file = "/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/prs/prs.afr.ukb450k.rda")

###
## PRS with EUR Reference Panel
list_prs_score.gbr <- system("ls /Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/prs/scores/gbr/*.profile | awk '{print $NF}'", 
                             intern = TRUE)

names_prs_score.gbr <- system("ls /Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/prs/scores/gbr/*.profile | awk '{print $NF}' | cut -d '.' -f7,8  | sed 's:has::g' | tr '.' '_' ", 
                              intern = TRUE)
## Read PRS files
d_eur <- fread(list_prs_score.gbr[1], header=T, select = c("ID2", "Profile_1"))

str(d_eur)

names(d_eur)

summary(d_eur)

names(d_eur) <- c("ID2", names_prs_score.gbr[1])

str(d_eur)

dat_eur <- d_eur; rm(d_eur)
gc()
##
for(i in 2:length(list_prs_score.gbr)){
  print(i)
  
  print(list_prs_score.gbr[i])
  
  print(names_prs_score.gbr[i])
  
  d <- fread(list_prs_score.gbr[i], header=T, select = c("ID2", "Profile_1"))
  
  names(d) <- c("ID2", names_prs_score.gbr[i])
  
  dat_eur <- merge(dat_eur, d, by="ID2")
  
  rm(d)
  gc()
}

str(dat_eur)

head(dat_eur)

###
fwrite(dat_eur, "/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/prs/prs.eur.ukb450k.csv.gz", 
       row.names = F, col.names = T, sep=",")

save(dat_eur, file = "/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/prs/prs.eur.ukb450k.rda")
###
##################
