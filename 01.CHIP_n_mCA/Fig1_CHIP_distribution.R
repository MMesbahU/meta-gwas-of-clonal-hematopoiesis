
##################### Figure 1 #######################################
#### CHIP prevalence increased with the age of the donor at the time of blood sampling.
############### Figure 1abcd
library(data.table)
library(ggplot2)
library(ggpubr)
library(cowplot)
theme_set(theme_cowplot())
############################### 
## Note: sort Gene name for fig 1c,1d
############# Table 1
######### TOPMed
### New data 2020 data
# topmed.chip.2020.table1 <- fread("/Users/muddin/Documents/Project/TopMed/TOPMed_CHIP_calls_8-31-2020/TOPMed_CHIPcalls_with_covariates_2020_08_31.tsv")
# topmed.var.2020.table1 <- fread("/Users/muddin/Documents/Project/TopMed/TOPMed_CHIP_calls_8-31-2020/TOPMed_variant_level_CHIPcalls_with_covariates_2020_08_31.tsv")
## Fig1abcd
# total N with CH data: 127946
# Total N with CH and Age data: 
# topmed2020.chip <- topmed.chip.2020.table1[,c(1,12,6,3,10)]  
# names(topmed2020.chip) <- c("SampleID", "Gender", "Age", 
#                             "Ancestry", "CHIP_status")
# topmed2020.chip$Cohort <- "TOPMed"
# 
# topmed2020.var <- topmed.var.2020.table1[, c(1,6,8)]
# names(topmed2020.var) <- c("SampleID", "Gene", "VAF")
# topmed2020.var$Cohort <- "TOPMed"
# 
### New data 2019 data
topmed.chip.2019.table1 <- fread("~/Documents/Project/TopMed/TOPMED_CHIP_100k_7-1-19/annot.TOPMed_100k_CHIP_calls_7_1_19_fordistribution.csv")
topmed.chip.2019.table1$hasCHIP <- ifelse(topmed.chip.2019.table1$nCHIP>0,1,0)
(table(topmed.chip.2019.table1$hasCHIP, exclude = NULL))
topmed.var.2019.table1 <- fread("~/Documents/Project/TopMed/TOPMED_CHIP_100k_7-1-19/TOPMed_100k_CHIP_variants_7_1_19_fordistribution.csv") 

## Fig1abcd
# total N with CH data: 87,116
topmed2019.chip <- topmed.chip.2019.table1[,c(1,34,4,35,7)]  
names(topmed2019.chip) <- c("SampleID", "Gender", "Age", 
                            "Ancestry", "CHIP_status")
topmed2019.chip$Cohort <- "TOPMed"

topmed2019.var <- topmed.var.2019.table1[, c(1,7,10)]
names(topmed2019.var) <- c("SampleID", "Gene", "VAF")
topmed2019.var$Cohort <- "TOPMed"

######### UKB 450k
# N = 454792 ; nrow(ukb450k.chip_table1) 
# ukb450k.chip_table1 <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/chipCall_CV_AB.ukb450k_all.tsv.gz")
# ukb450k.var_table1 <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/chipVar_CV_AB.ukb450k_all.csv.gz")
ukb450k.chip_table1 <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/chip_call/ukb450k.CHIP_phewas.csv.gz")
ukb450k.var_table1 <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/chip_call/ukb450k.variants_phewas.csv.gz")
load("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/ukb500k_covariates.rda");rm(d_base,d_pcs)
ukb450k.chip_table1 <- merge(ukb450k.chip_table1, d_ukb_all, 
                             by.x="eid_7089", by.y = "eid"); rm(d_ukb_all)
## Fig1abcd
# ukb450k.chip <- ukb450k.chip_table1[,c(3,17,16,43,29)]  
ukb450k.chip <- ukb450k.chip_table1[,c(1,2,3,35,5)]
names(ukb450k.chip) <- c("SampleID", "Gender", "Age", 
                         "Ancestry", "CHIP_status")
ukb450k.chip$Cohort <- "UKB"

# ukb450k.var <- ukb450k.var_table1[, c(1,2,54)]
ukb450k.var <- ukb450k.var_table1[, c(1,2,3)]; rm(ukb450k.var_table1)

names(ukb450k.var) <- c("SampleID", "Gene", "VAF")
ukb450k.var$Cohort <- "UKB"
######### MGBB 53k
######## data for Fig1abcd
# N = 52966 ; nrow(mgbb53k.chip_table1) 
mgbb53k.chip_table1 <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.chip_call_all52966.tsv.gz")
mgbb53k.var_table1 <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.chip_var_N5910_nchip7145.csv.gz")
## Fig1abcd
mgbb53k.chip <- mgbb53k.chip_table1[,c(3,5,6,14,16)]  
names(mgbb53k.chip) <- c("SampleID", "Gender", "Age", 
                         "Ancestry", "CHIP_status")
mgbb53k.chip$Cohort <- "MGBB"

mgbb53k.var <- mgbb53k.var_table1[, c(1,14,9)]
names(mgbb53k.var) <- c("SampleID", "Gene", "VAF")
mgbb53k.var$Cohort <- "MGBB"

############# Plot Figure 1
# topmed2020.ukbb.mgbb <- rbind.data.frame(topmed2020.chip, ukb450k.chip,
#                                          mgbb53k.chip)
# nrow(topmed2020.ukbb.mgbb)
# # Total CHIP samples analysed: 635,704
# table(topmed2020.ukbb.mgbb$Cohort)
# # MGBB TOPMed    UKB 
# # 52966 127946 454792
# table(topmed2020.ukbb.mgbb$CHIP_status)
# # 0      1 
# # 595191  40513 
# prop.table(table(topmed2020.ukbb.mgbb$CHIP_status))
# # 0          1 
# # 0.93627065 0.06372935 
# topmed2020.ukbb.mgbb$Cohort <- factor(topmed2020.ukbb.mgbb$Cohort, 
#                                       levels = c("UKB", "TOPMed", "MGBB"))
# 2019
topmed2019.ukbb.mgbb <- rbind.data.frame(topmed2019.chip, ukb450k.chip,
                                         mgbb53k.chip)
nrow(topmed2019.ukbb.mgbb)
# Total CHIP samples analysed: 594,409
# total_samples = 454327 + 87116 + 54583 + 52966
table(topmed2019.ukbb.mgbb$Cohort)
# MGBB TOPMed    UKB 
# 52966 87116 454327
(table(topmed2019.ukbb.mgbb$CHIP_status))
# 0      1 
# 555918  38491 
prop.table(table(topmed2019.ukbb.mgbb$CHIP_status))
# 0         1 
# 0.93524492 0.06475508

topmed2019.ukbb.mgbb$Cohort <- factor(topmed2019.ukbb.mgbb$Cohort, 
                                      levels = c("UKB", "TOPMed", "MGBB"))


## var
# Var.topmed2020.ukbb.mgbb <- rbind.data.frame(topmed2020.var, ukb450k.var,
#                                              mgbb53k.var)
# Var.topmed2020.ukbb.mgbb$Cohort <- factor(Var.topmed2020.ukbb.mgbb$Cohort, 
#                                           levels = c("UKB", "TOPMed", "MGBB"))
# 2019
Var.topmed2019.ukbb.mgbb <- rbind.data.frame(topmed2019.var, ukb450k.var,
                                             mgbb53k.var)
Var.topmed2019.ukbb.mgbb$Cohort <- factor(Var.topmed2019.ukbb.mgbb$Cohort, 
                                          levels = c("UKB", "TOPMed", "MGBB"))

## number of CHIP mutation per sample
count.topmed2019 <- as.data.frame(round(prop.table(table(table(topmed2019.var$SampleID)))
                                        *100,1))
count.topmed2019$Cohort <- "TOPMed"
# count.topmed2020 <- as.data.frame(round(prop.table(table(table(topmed2020.var$SampleID)))
#                                         *100,1))
# count.topmed2020$Cohort <- "TOPMed"

count.ukb450k <- as.data.frame(round(prop.table(table(table(ukb450k.var$SampleID)))
                                     *100,1))
count.ukb450k$Cohort <- "UKB"

count.mgbb53k <- as.data.frame(round(prop.table(table(table(mgbb53k.var$SampleID)))
                                     *100,1))
count.mgbb53k$Cohort <- "MGBB"

# chip_per_sample_2020 <- rbind(count.topmed2020,
#                               count.ukb450k, 
#                               count.mgbb53k)
# chip_per_sample_2020$Cohort <- factor(chip_per_sample_2020$Cohort, 
#                                       levels = c("UKB", "TOPMed", 
#                                                  "MGBB"))
chip_per_sample_2019 <- rbind(count.topmed2019,
                              count.ukb450k, 
                              count.mgbb53k)
chip_per_sample_2019$Cohort <- factor(chip_per_sample_2019$Cohort, 
                                      levels = c("UKB", "TOPMed", 
                                                 "MGBB"))

# top 10 Genes (proportion)
# CHIP_Gene_Table.topmed2020 <- as.data.frame(round(prop.table(table(topmed2020.var$Gene))
#                                                   *100,1),
#                                             stringsAsFactors = F)
# CHIP_Gene_Table.topmed2020$Cohort <- "TOPMed"
CHIP_Gene_Table.topmed2019 <- as.data.frame(round(prop.table(table(topmed2019.var$Gene))
                                                  *100,1),
                                            stringsAsFactors = F)
CHIP_Gene_Table.topmed2019$Cohort <- "TOPMed"

CHIP_Gene_Table.ukbb <- as.data.frame(round(prop.table(table(ukb450k.var$Gene))
                                            *100,1),
                                      stringsAsFactors = F)
CHIP_Gene_Table.ukbb$Cohort <- "UKB"

CHIP_Gene_Table.mgb <- as.data.frame(round(prop.table(table(mgbb53k.var$Gene))
                                           *100,1),
                                     stringsAsFactors = F)
CHIP_Gene_Table.mgb$Cohort <- "MGBB"

# CHIP_Gene_Table.2020 <- rbind.data.frame(CHIP_Gene_Table.topmed2020, CHIP_Gene_Table.ukbb, CHIP_Gene_Table.mgb)
# CHIP_Gene_Table.2020$Cohort <- factor(CHIP_Gene_Table.2020$Cohort, 
#                                       levels = c("UKB", "TOPMed", "MGBB"))

CHIP_Gene_Table.2019 <- rbind.data.frame(CHIP_Gene_Table.topmed2019, CHIP_Gene_Table.ukbb, CHIP_Gene_Table.mgb)
CHIP_Gene_Table.2019$Cohort <- factor(CHIP_Gene_Table.2019$Cohort, 
                                      levels = c("UKB", "TOPMed", "MGBB"))

#  top 10 chip genes in each cohort
# top10_genes.2020 <- unique(c(head(CHIP_Gene_Table.topmed2020[order(CHIP_Gene_Table.topmed2020$Freq, decreasing = T),], n = 10)$Var1,
#                              head(CHIP_Gene_Table.ukbb[order(CHIP_Gene_Table.ukbb$Freq, decreasing = T),], n = 10)$Var1,
#                              head(CHIP_Gene_Table.mgb[order(CHIP_Gene_Table.mgb$Freq, decreasing = T),], n = 10)$Var1))
# 
# CHIP_Gene_Table.2020 <- subset(CHIP_Gene_Table.2020, 
#                                CHIP_Gene_Table.2020$Var1 %in% 
#                                  top10_genes.2020)

top10_genes.2019 <- unique(c(head(CHIP_Gene_Table.topmed2019[order(CHIP_Gene_Table.topmed2019$Freq, decreasing = T),], n = 10)$Var1,
                             head(CHIP_Gene_Table.ukbb[order(CHIP_Gene_Table.ukbb$Freq, decreasing = T),], n = 10)$Var1,
                             head(CHIP_Gene_Table.mgb[order(CHIP_Gene_Table.mgb$Freq, decreasing = T),], n = 10)$Var1))

CHIP_Gene_Table.2019 <- subset(CHIP_Gene_Table.2019, 
                               CHIP_Gene_Table.2019$Var1 %in% 
                                 top10_genes.2019)


# ## VAF distribution in top 10 genes per cohort
# top10_genes.topmed2020.ukbb.mgbb <- subset(Var.topmed2020.ukbb.mgbb,
#                                            Var.topmed2020.ukbb.mgbb$Gene %in% 
#                                              top10_genes.2020)
# top10_genes.topmed2020.ukbb.mgbb$Cohort <- factor(top10_genes.topmed2020.ukbb.mgbb$Cohort, 
#                                                   levels = c("UKB", "TOPMed", "MGBB"))

top10_genes.topmed2019.ukbb.mgbb <- subset(Var.topmed2019.ukbb.mgbb,
                                           Var.topmed2019.ukbb.mgbb$Gene %in% 
                                             top10_genes.2019)
top10_genes.topmed2019.ukbb.mgbb$Cohort <- factor(top10_genes.topmed2019.ukbb.mgbb$Cohort, 
                                                  levels = c("UKB", "TOPMed", "MGBB"))

####################
  # current files
# save.image(file = "/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/rerun/fig1.data_TOPMed2019_UKB200k_UKB250k_MGBB53k.rda")
load("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/rerun/fig1.data_TOPMed2019_UKB200k_UKB250k_MGBB53k.rda")
  # old, ukb450k
## save.image(file = "/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/rerun/fig1.data_TOPMed_UKB450k_MGBB53k.rda")
## load("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/rerun/fig1.data_TOPMed_UKB450k_MGBB53k.rda")

############################################################
# Fig 1a Prevalence
# p1_2020 <- ggplot(data=topmed2020.ukbb.mgbb, aes(x=Age, y=CHIP_status, group=Cohort)) + 
#   geom_smooth(aes(colour=Cohort), method ="glm", 
#               method.args = list(family = "binomial")) + 
#   xlab("Age") + ylab("Prevalence") + ggtitle("a") + 
#   scale_y_continuous(breaks = c(seq(from = 0, to = 0.6, by = 0.1)) ) + scale_x_continuous(breaks = c(seq(0, 100, 10) )) + theme(legend.position = "", plot.title = element_text(size = 20, face = "bold"))
# 
p1_2019 <- ggplot(data=topmed2019.ukbb.mgbb, aes(x=Age, y=CHIP_status, group=Cohort)) + 
  geom_smooth(aes(colour=Cohort), method ="glm", 
              method.args = list(family = "binomial")) + 
  xlab("Age") + ylab("Prevalence") + ggtitle("a") + 
  scale_y_continuous(breaks = c(seq(from = 0, to = 0.6, by = 0.1)) ) + scale_x_continuous(breaks = c(seq(0, 100, 10) )) + theme(legend.position = "", plot.title = element_text(size = 20, face = "bold"))

# Fig 1b
# p2_2020 <- ggplot(data=chip_per_sample_2020, aes(x=reorder(Var1, -Freq),
#                                                  y=Freq, fill=Cohort)) +
#   xlab("Number of CHIP mutations") +
#   ylab(label = "Proportion of Individuals (%)") +
#   geom_bar(stat="identity", width=0.8, position=position_dodge()) +
#   geom_text(aes(label=Freq), vjust=-0.5, color="black",
#             position = position_dodge(0.9), size=3) +
#   theme(legend.title = element_blank(), legend.position = "right", plot.title = element_text(size = 20, face = "bold")) + ggtitle("b") 
# 
p2_2019 <- ggplot(data=chip_per_sample_2019, aes(x=reorder(Var1, -Freq),
                                                 y=Freq, fill=Cohort)) +
  xlab("Number of CHIP mutations") +
  ylab(label = "Proportion of Individuals (%)") +
  geom_bar(stat="identity", width=0.8, position=position_dodge()) +
  geom_text(aes(label=Freq), vjust=-0.5, color="black",
            position = position_dodge(0.9), size=3) +
  theme(legend.title = element_blank(), legend.position = "right", plot.title = element_text(size = 20, face = "bold")) + ggtitle("b") 

# Fig 1c
# p3_2020 <- ggplot(data=CHIP_Gene_Table.2020, aes(x=reorder(Var1, -Freq),
#                                                  y=Freq, fill=Cohort)) +
#   xlab("") +
#   ylab(label = "Proportion of Individuals (%)") +
#   geom_bar(stat="identity", width=0.8, position=position_dodge()) +
#   geom_text(aes(label=Freq), vjust=-0.5, color="black",
#             position = position_dodge(0.9), size=3) +
#   theme(axis.text.x = element_text(angle = 90,
#                                    vjust = 1,  hjust=1),
#         legend.position = "", 
#         plot.title = element_text(size = 20, face = "bold")) +
#   ggtitle("c") 

p3_2019 <- ggplot(data=CHIP_Gene_Table.2019, aes(x=reorder(Var1, -Freq),
                                                 y=Freq, fill=Cohort)) +
  xlab("") +
  ylab(label = "Proportion of Individuals (%)") +
  geom_bar(stat="identity", width=0.8, position=position_dodge()) +
  geom_text(aes(label=Freq), vjust=-0.5, color="black",
            position = position_dodge(0.9), size=3) +
  theme(axis.text.x = element_text(angle = 45,
                                   vjust = 1,  hjust=1),
        legend.position = "", 
        plot.title = element_text(size = 20, face = "bold")) +
  ggtitle("c") 

# Fig 1d
# top10_genes.topmed2020.ukbb.mgbb$Gene <- factor(top10_genes.topmed2020.ukbb.mgbb$Gene, 
#                                                 levels = c("DNMT3A", "TET2", "ASXL1","PPM1D","TP53","SF3B1","SRSF2","JAK2",
#                                                            "ZNF318","ZBTB33","YLPM1", "NF1", "ASXL2") )
# p4_2020 <- ggplot(data=top10_genes.topmed2020.ukbb.mgbb, aes(x=Gene, y=VAF, fill=Cohort)) + 
#   xlab("") + ylab("Variant Allele Fraction") + 
#   geom_violin(trim = FALSE) + scale_y_log10(breaks = c(0.02,0.1,.2,.3, 0.5,1)) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1,  hjust=1),
#         legend.title = element_blank(), 
#         plot.title = element_text(size = 20, face = "bold")) + 
#   ggtitle("d") +
#   stat_summary(fun = "median", geom = "point",
#                color = "white", 
#                position = position_dodge(0.9))
# 
top10_genes.topmed2019.ukbb.mgbb$Gene <- factor(top10_genes.topmed2019.ukbb.mgbb$Gene, 
                                                levels = c("DNMT3A","TET2", "ASXL1", "PPM1D", "TP53", 
                                                           "JAK2", "SF3B1", "SRSF2", "ZNF318", "YLPM1", 
                                                           "GNB1", "GNAS", "ASXL2") )
p4_2019 <- ggplot(data=top10_genes.topmed2019.ukbb.mgbb, aes(x=Gene, y=VAF, fill=Cohort)) + 
  xlab("") + ylab("Variant Allele Fraction") + 
  geom_violin(trim = FALSE) + scale_y_log10(breaks = c(0.02,0.1,.2,.3, 0.5,1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1,  hjust=1),
        legend.title = element_blank(), 
        plot.title = element_text(size = 20, face = "bold")) + 
  ggtitle("d") +
  stat_summary(fun = "median", geom = "point",
               color = "white", 
               position = position_dodge(0.9))


##############################################################
# 
# png("~/Documents/Project/CHIP_GWAS/rerun/Figures/Fig1abcd.chip_distribution.topmed2020_ukb450k_mgbb53k.6Aug2022.png",
#     width=16, height=10, units= "in", res=300, pointsize = 4)
# ggarrange(p1_2020,p2_2020,p3_2020, p4_2020, 
#           ncol = 2, nrow = 2 )
# dev.off()

png("~/Documents/Project/CHIP_GWAS/rerun/Figures/Fig1abcd.chip_distribution.topmed2019_ukb450k_mgbb53k.28Aug2022.png",
    width=16, height=10, units= "in", res=300, pointsize = 4)
ggarrange(p1_2019,p2_2019,p3_2019, p4_2019, 
          ncol = 2, nrow = 2 )
dev.off()

######################################################################
######################################################################
table(ukb450k.chip_table1$p22006=="Caucasian"&
        ukb450k.chip_table1$p21000_i0=="British")
# UKB: Others=55453 + WB=380801 + Black=7094 + Asian=10979
table(topmed.chip.2019.table1$Race)
total.EUR = sum(c(29036,380801,44563,54583))
sort(table(ukb450k.chip_table1$p21000_i0),decreasing = T)
sort(table(mgbb53k.chip_table1$Ancestry_Self_cat),decreasing = T)
total.asian= sum(c(806, 5326+1630+1620+1415+740+207+41, 1577,0) )
# Indian,Any other Asian background, Pakistani, Chinese,  
# White and Asian,Bangladeshi, Asian or Asian British
total.black= sum(c(12149, 4012+2945+114 + 23, 2665, 0) )
# Caribbean,African, Any other Black background, 
# Black or Black British
total.others <- sum(c(43887+56+4+28+1150, 20233+ 14708+11897+ 4067 + 
                        1479 +920 +562 +504 +485 +373 + 186 + 39, 4161, 0))
# Any other white background, Irish,  Other ethnic group,
# Prefer not to answer, Any other mixed background, White and Black Caribbean,
# White, "", White and Black African, Do not know, Mixed
sum(total.asian,total.black,total.others)/sum(total.EUR, total.asian,total.black,total.others)

round(sum(total.asian,total.black,total.others)/sum(total.EUR, total.asian,total.black,total.others) *100,1)
# 21.6

###### CHIP statistics
# Total CHIP mutations
total_chipmutation = nrow(Var.topmed2019.ukbb.mgbb) + 1583
# nrow(topmed2019.var) + nrow(ukb450k.var) + nrow(mgbb53k.var_table1) + 1583
# 44878
# total unique samples
total_samp = length(unique(Var.topmed2019.ukbb.mgbb$SampleID)) + 1583
# length(unique(topmed2019.var$SampleID)) + length(unique(ukb450k.var$SampleID)) + length(unique(mgbb53k.var_table1$sample_id)) + 1583
# 40638  

## Overall CHIP Prevalence
prop.table(table(topmed2019.ukbb.mgbb$CHIP_status))
# 0          1 
# 0.93524492 0.06475508
table(topmed2019.ukbb.mgbb$CHIP_status)

sample_chip_prev <- (38491+1583)/(555918  + 38491 + 54583)
#  0.06174806
# ## ## Average prev
# ()/3
# biovu=1583/54583
# mean()
# sd()
#################
## Case Control in GWAS
###### CHIP case control counts
pheno.mgbb_gwas <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.21Jul2022.tsv")
pheno.ukb200k_gwas <- fread("~/Documents/Project/CHIP_GWAS/Aug11_2021.GWAS/UKB200k_CHIP_PhenoCovar.11Aug2021.tsv.gz")

pheno.ukb250k_gwas <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/chipCallnomiss_CV_AB.ukb250k_all.tsv")

# CHIP Bick et al GWAS:
cat ("TOPMed cases:",3831 ,"Controls:",65404-3831 )
# UKB 450k
multi.chip_ukb_gwas_ids <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/multi_ancestry_summary/ukb450k.multiAncestry.CHIP.ids")
table(pheno.ukb200k_gwas$hasCHIP[pheno.ukb200k_gwas$eid_7089%in% multi.chip_ukb_gwas_ids$V2])
# 0      1 
# 180105  10986 
table(pheno.ukb250k_gwas$hasCHIP[pheno.ukb250k_gwas$eid %in% multi.chip_ukb_gwas_ids$V2])
# 0      1 
# 219112  15640
cat ("UKBB cases:", sum(10986,15640),"Controls:",
     sum(180105,219112), "N=", sum(180105,10986,219112,15640) )

multi.chip_mgbb_gwas_ids <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/multi_ancestry_summary/mgbb53k.multiAncestry.CHIP.ids")
table(pheno.mgbb_gwas$hasCHIP[pheno.mgbb_gwas$IID %in% multi.chip_mgbb_gwas_ids$V2])
# 0     1 
# 44755  5550
cat ("MGBB cases:", 5550,"Controls:",44755, "N=", 44755+5550)

cat ("BioVU cases:",1583 ,"Controls:",54223-1583 )

cat("CHIP Cases:",sum(c(3831, sum(10986,15640), 5550, 1583)))
# CHIP Cases: 37590 # 17044 (old)
cat("CHIP Controls:",sum(c(65404-3831, sum(180105,219112), 44755, 54223-1583)))
# CHIP Controls: 558185 # 306068
# CHIP GWAS N = 37590 + 558185 = 595775

# Ancestry: 
table(pheno.ukb250k_gwas$Ethnicity[pheno.ukb250k_gwas$eid %in% multi.chip_ukb_gwas_ids$V2])
table(pheno.ukb200k_gwas$Ethnic_Background[pheno.ukb200k_gwas$eid_7089 %in% multi.chip_ukb_gwas_ids$V2])
table(pheno.mgbb_gwas$Ancestry_Self_cat[pheno.mgbb_gwas$IID %in% multi.chip_mgbb_gwas_ids$V2])

# Non-EUR = sum(65404-29373, 234752 - 196951,191091-167713, 50305 - 42327, 0)  
nonEUR <- 105188
round(105188/595775 * 100, 1)

## VAF 
round(median(Var.topmed2019.ukbb.mgbb$VAF[Var.topmed2019.ukbb.mgbb$Gene=="DNMT3A" & Var.topmed2019.ukbb.mgbb$Cohort=="UKB"]),3)
round(median(Var.topmed2019.ukbb.mgbb$VAF[Var.topmed2019.ukbb.mgbb$Gene=="DNMT3A" & Var.topmed2019.ukbb.mgbb$Cohort=="TOPMed"]),3)
round(median(Var.topmed2019.ukbb.mgbb$VAF[Var.topmed2019.ukbb.mgbb$Gene=="DNMT3A" & Var.topmed2019.ukbb.mgbb$Cohort=="MGBB"]),3)

round(median(Var.topmed2019.ukbb.mgbb$VAF[Var.topmed2019.ukbb.mgbb$Gene=="TET2" & Var.topmed2019.ukbb.mgbb$Cohort=="UKB"]),3)
round(median(Var.topmed2019.ukbb.mgbb$VAF[Var.topmed2019.ukbb.mgbb$Gene=="TET2" & Var.topmed2019.ukbb.mgbb$Cohort=="TOPMed"]),3)
round(median(Var.topmed2019.ukbb.mgbb$VAF[Var.topmed2019.ukbb.mgbb$Gene=="TET2" & Var.topmed2019.ukbb.mgbb$Cohort=="MGBB"]),3)

round(median(Var.topmed2019.ukbb.mgbb$VAF[Var.topmed2019.ukbb.mgbb$Gene=="JAK2" & Var.topmed2019.ukbb.mgbb$Cohort=="UKB"]),3)
round(median(Var.topmed2019.ukbb.mgbb$VAF[Var.topmed2019.ukbb.mgbb$Gene=="JAK2" & Var.topmed2019.ukbb.mgbb$Cohort=="TOPMed"]),3)
round(median(Var.topmed2019.ukbb.mgbb$VAF[Var.topmed2019.ukbb.mgbb$Gene=="JAK2" & Var.topmed2019.ukbb.mgbb$Cohort=="MGBB"]),3)


###
# DNMT3A: 
cat ("TOPMed cases:",1269 ,"Controls:",63769-1269 )
# UKB450k 
multi.dnmt3a_ukb_gwas_ids <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/multi_ancestry_summary/ukb450k.multiAncestry.DNMT3A.ids")
table(pheno.ukb200k_gwas$hasDNMT3A[pheno.ukb200k_gwas$eid_7089%in% multi.dnmt3a_ukb_gwas_ids$V2])
# 0      1 
# 180105   6748 
table(pheno.ukb250k_gwas$hasDNMT3A[pheno.ukb250k_gwas$eid %in% multi.dnmt3a_ukb_gwas_ids$V2])
# 0      1 
# 219112   7963
# UKB DNMT3A: (6748 + 7963) + (180105 +219112) = 413928  

multi.dnmt3a_mgbb_gwas_ids <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/multi_ancestry_summary/mgbb53k.multiAncestry.DNMT3A.ids")
table(pheno.mgbb_gwas$hasDNMT3A[pheno.mgbb_gwas$IID %in% multi.dnmt3a_mgbb_gwas_ids$V2])
# 0     1 
# 44755  2559
# MGB DNMT3A: 2559 + 44755 = 47314
# BioVU DNMT3A Cases: 607 + 53616 controls = 54223
cat("DNMT3A Cases:",sum(c(1269,(6748 + 7963),2559,607)))
# DNMT3A Cases: 19146
cat("DNMT3A Controls:",sum(c(63769-1269,(180105 +219112),44755,53616)))
# DNMT3A Controls: 560088

# TET2
cat ("TOPMed cases:",462 ,"Controls:",62250-462 )
# 
multi.tet2_ukb_gwas_ids <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/multi_ancestry_summary/ukb450k.multiAncestry.TET2.ids")
table(pheno.ukb200k_gwas$hasTET2[pheno.ukb200k_gwas$eid_7089%in% multi.tet2_ukb_gwas_ids$V2])
# 0      1 
# 180105   1906 
table(pheno.ukb250k_gwas$hasTET2[pheno.ukb250k_gwas$eid %in% multi.tet2_ukb_gwas_ids$V2])
# 0      1 
# 219112   3346

# UKB TET2: (1906 + 3346) + (180105 + 219112) = 404469  

multi.tet2_mgbb_gwas_ids <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/multi_ancestry_summary/mgbb53k.multiAncestry.TET2.ids")
table(pheno.mgbb_gwas$hasTET2[pheno.mgbb_gwas$IID %in% multi.tet2_mgbb_gwas_ids$V2])
# 0     1 
# 44755  1240 
## MGB TET2: 1240 + 44755 = 45995
## BioVU TET2: 339 + 53884 controls = 54223
cat("TET2 Cases:",sum(c(462,(1906 + 3346),1240,339)))
# TET2 Cases: 7293
cat("TET2 Controls:",sum(c(61788,(180105 + 219112) ,44755,53884)))
# TET2 Controls: 559644


## CHIP
round(26626/425843,4)
# pop.chip: 0.0625
round(37590 /(558185 +37590),4)
# 0.0631
### DNMT3A
# round((6748 + 7963)/413928,4)
# pop.dnmt3a: 0.036
round(19146/(560088+19146),4)
# 0.0331
### TET2 Prevalence
round((1906 + 3346)/404469,4)
# pop.tet2 : 0.013
round(7293/(559644+7293),4)
# 0.0129

## EUR GWAS
###########
eur.chip_ukb_gwas_ids <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/eur_only_summary/ukb450k.eur.CHIP.ids")
round(prop.table(table(c(pheno.ukb200k_gwas$hasCHIP[pheno.ukb200k_gwas$IID %in% eur.chip_ukb_gwas_ids$V2],
                       pheno.ukb250k_gwas$hasCHIP[pheno.ukb250k_gwas$IID %in% eur.chip_ukb_gwas_ids$V2])))*100,2)
# 0      1 
# 341652  23012
# 6.31
eur.dnmt3a_ukb_gwas_ids <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/eur_only_summary/ukb450k.eur.DNMT3A.ids")
round(prop.table(table(c(pheno.ukb200k_gwas$hasDNMT3A[pheno.ukb200k_gwas$IID %in% eur.dnmt3a_ukb_gwas_ids$V2],
                         pheno.ukb250k_gwas$hasDNMT3A[pheno.ukb250k_gwas$IID %in% eur.dnmt3a_ukb_gwas_ids$V2])))*100,2)
# 0      1 
# 341652  12737
# 3.59
eur.tet2_ukb_gwas_ids <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/eur_only_summary/ukb450k.eur.TET2.ids")
round(prop.table(table(c(pheno.ukb200k_gwas$hasTET2[pheno.ukb200k_gwas$IID %in% eur.tet2_ukb_gwas_ids$V2],
                         pheno.ukb250k_gwas$hasTET2[pheno.ukb250k_gwas$IID %in% eur.tet2_ukb_gwas_ids$V2])))*100,2)
# 0      1 
# 341652   4537
# 1.31
  # EUR MGBB
eur.chip_mgbb_gwas_ids <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/eur_only_summary/mgbb53k.eur.CHIP.ids")
round(prop.table(table(pheno.mgbb_gwas$hasCHIP[pheno.mgbb_gwas$IID %in% eur.chip_mgbb_gwas_ids$V2]))*100,2)
# 0     1 
# 88.14 11.86
eur.dnmt3a_mgbb_gwas_ids <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/eur_only_summary/mgbb53k.eur.DNMT3A.ids")
round(prop.table(table(pheno.mgbb_gwas$hasDNMT3A[pheno.mgbb_gwas$IID %in% eur.dnmt3a_mgbb_gwas_ids$V2]))*100,2)
# 0     1 
# 94.23  5.77
eur.tet2_mgbb_gwas_ids <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/eur_only_summary/mgbb53k.eur.TET2.ids")
round(prop.table(table(pheno.mgbb_gwas$hasTET2[pheno.mgbb_gwas$IID %in% eur.tet2_mgbb_gwas_ids$V2]))*100,2)
# 0     1 
# 97.02  2.98

###########################

####### Multi-ancestry GWAS
## CHIP
round(26626/425843,4)
# pop.chip: 0.0625
round(37590 /(558185 +37590),4)
# 0.0631
### DNMT3A
round((6748 + 7963)/413928,4)
# pop.dnmt3a: 0.0355
round(19146/(560088+19146),4)
# 0.0331
### TET2 Prevalence
round((1906 + 3346)/404469,4)
# pop.tet2 : 0.013
round(7293/(559644+7293),4)
# 0.0129

########################## EUR GWAS
  #### Overall CHIP EUR
# Pop Prev. (i.e. in UKB):  
round( 23012/(341652 + 23012 ) *100 , 2)
# 6.31%
# Sample Prev. CHIP (TOPMED + UKBB + MGBB + BioVU)
## 3831 * (29036/65404)
round(( 1700 + 23012 + 5018 + 1583) /(29036 + 341652 + 23012 + 37309 + 5018 + 54223) *100,2)
# 6.39%

  #### DNMT3A EUR
# Pop Prev. (i.e. in UKB):  
round( 12737/( 341652 + 12737) *100 , 2)
# 3.59
# Sample Prev. DNMT3A (UKBB + MGBB + BioVU)
round(( 12737 + 2283 + 607) /( 341652+12737 + 37309+2283 + 54223) *100,2)
# 3.49%

  #### TET2 EUR
  # Pop Prev. (i.e. in UKB):  
round( 4537/( 341652  +4537) *100 , 2)
# 1.31
  # Sample Prev. TET2 (UKBB + MGBB + BioVU)
round((4537 + 1145 + 339) /(341652 +  4537 + 37309 +1145 + 54223) *100,2)
# 1.37%
#############
## Table S1 
table(topmed2019.ukbb.mgbb$Cohort)
# MGBB TOPMed    UKB 
# 52966  87116 454327 
cat("UKB: Age ",round(mean(topmed2019.ukbb.mgbb$Age[topmed2019.ukbb.mgbb$Cohort=="UKB"]))
,"(SD=",round(sd(topmed2019.ukbb.mgbb$Age[topmed2019.ukbb.mgbb$Cohort=="UKB"])),")")
# UKB: Age  57 (SD= 8 )
cat("MGBB: Age ",round(mean(topmed2019.ukbb.mgbb$Age[topmed2019.ukbb.mgbb$Cohort=="MGBB"]))
    ,"(SD=",round(sd(topmed2019.ukbb.mgbb$Age[topmed2019.ukbb.mgbb$Cohort=="MGBB"])),")")
# MGBB: Age  54 (SD= 17 )