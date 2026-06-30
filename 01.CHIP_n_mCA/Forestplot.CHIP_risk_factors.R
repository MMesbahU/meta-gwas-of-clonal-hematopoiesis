##############################
#### Forest plot
##############################
library(data.table) # version 1.14.6
library(meta) # version 6.2-1
library(grid) # version 4.2.2
library(scales) # version 1.2.1
setwd("~/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/NG_v2/")
#---------------
## AoU ####
final_results <- fread("Figs/aou250k.risk_factors_for_CHIP.final_glm_results.csv", header=T)
head(final_results)
final_results$V1 <- NULL
final_results$Pval <- formatC(x = final_results$P, digits = 1, format = "E")
final_results$OR <- formatC(x = round(exp(final_results$Beta),2), digits = 2, format = "f")
final_results$lSE <- final_results$Beta - 1.96*final_results$SE
final_results$uSE <- final_results$Beta + 1.96*final_results$SE
final_results$CI95 <- paste0("[", formatC(x = round(exp(final_results$lSE),2), digits = 2, format = "f"),
                             "-",
                             formatC(x = round(exp(final_results$uSE),2), digits = 2, format = "f"),
                             "]")
head(final_results)

##
final_results$Outcome <- ordered(final_results$Outcome, 
                                 levels=c("Overall CHIP",
                                          "Expanded CHIP", 
                                          "DNMT3A","TET2",
                                          "ASXL1","Splicing Factors",
                                          "DNA damage"))


## w/o ancestry
nonAnc.final_results <- final_results %>% 
  filter(!(Exposure %in% c("AFR", "AMR", "SAS", "EAS", "Other")) )
m <- metagen(TE = Beta,
                lower = lSE,
                upper = uSE,
                studlab = Outcome,
                subgroup = Exposure,
                data = nonAnc.final_results,
                sm="OR")

m$N <- nonAnc.final_results$N
m$Cases <- nonAnc.final_results$Cases
m$Controls <- nonAnc.final_results$Controls
## 
pdf("Figs/AoU.CHIP_Risk.forest_plot_CHIP.noANC.pdf", width = 12, height = 8)
forest(x = m, 
       common=F, 
       random=F, 
       hetstat=F, 
       subgroup=k.w>=1, 
       weight.study="same",  
       level=0.95, 
       xlim=c(0.5, 3), 
       smlab="Effect of Exposures\non CHIP\n", 
       smlab.pos=0, 
       colgap=unit(7, "mm"),
       xlab="Odds Ratio (95% CI)", 
       squaresize=0.6, 
       col.subgroup="black", 
       colgap.left=unit(0.1,"cm"),
       colgap.forest.left="3mm", 
       colgap.forest.right="2mm", 
       leftcols=c("studlab"), 
       leftlabs = c("                     "),
       rightcols=c("OR","CI95","Pval", "N", "Cases", "Controls"),
       rightlabs=c("OR","95% CI", "P", "N", "Cases", "Controls"),
       col.inside="black", 
       plotwidth=unit(6.5, "cm"), 
       print.subgroup.name=F)

dev.off()

## Ancestry assoc
##
Anc.final_results <- final_results %>% 
  filter((Exposure %in% c("AFR", "AMR", "SAS", "EAS", "Other")) )

m <- metagen(TE = Beta,
             lower = lSE,
             upper = uSE,
             studlab = Outcome,
             subgroup = Exposure,
             data = Anc.final_results,
             sm="OR")

m$N <- Anc.final_results$N
m$Cases <- Anc.final_results$Cases
m$Controls <- Anc.final_results$Controls


## 
pdf("Figs/AoU.CHIP_Risk.forest_plot_CHIP.ANC.pdf", width = 12, height = 12)
forest(x = m, 
       common=F, 
       random=F, 
       hetstat=F, 
       subgroup=k.w>=1, 
       weight.study="same",  
       level=0.95, 
       xlim=c(0.5, 3), 
       smlab="Effect of Exposures\non CHIP\n", 
       smlab.pos=0, 
       colgap=unit(7, "mm"),
       xlab="Odds Ratio (95% CI)", 
       squaresize=0.6, 
       col.subgroup="black", 
       colgap.left=unit(0.1,"cm"),
       colgap.forest.left="3mm", 
       colgap.forest.right="2mm", 
       leftcols=c("studlab"), 
       leftlabs = c("                     "),
       rightcols=c("OR","CI95","Pval", "N", "Cases", "Controls"),
       rightlabs=c("OR","95% CI", "P", "N", "Cases", "Controls"),
       col.inside="black", 
       plotwidth=unit(6.5, "cm"), 
       print.subgroup.name=F)

dev.off()

#--------------------------------------

#-------- MGBB ####
## MGBB ####
final_results <- fread("Figs/MGBB_CHIP.risk_factors.age2adj.results_flag_vars.csv", header=T)
head(final_results)
final_results$Pval <- formatC(x = final_results$P, digits = 1, format = "E")
final_results$OR <- formatC(x = round(exp(final_results$Beta),2), digits = 2, format = "f")
final_results$lSE <- final_results$Beta - 1.96*final_results$SE
final_results$uSE <- final_results$Beta + 1.96*final_results$SE
final_results$CI95 <- paste0("[", formatC(x = round(exp(final_results$lSE),2), digits = 2, format = "f"),
                             "-",
                             formatC(x = round(exp(final_results$uSE),2), digits = 2, format = "f"),
                             "]")
head(final_results)

##
final_results$Outcome <- ordered(final_results$Outcome, 
                                 levels=c("Overall CHIP",
                                          "Expanded CHIP", 
                                          "DNMT3A","TET2",
                                          "ASXL1","Splicing Factors",
                                          "DNA damage"))


## w/o ancestry
nonAnc.final_results <- final_results %>% 
  filter(!(Exposure %in% c("AFR", "AMR", "SAS", "EAS")) )
m <- metagen(TE = Beta,
             lower = lSE,
             upper = uSE,
             studlab = Outcome,
             subgroup = Exposure,
             data = nonAnc.final_results,
             sm="OR")

m$N <- nonAnc.final_results$N
m$Cases <- nonAnc.final_results$Cases
m$Controls <- nonAnc.final_results$Controls
## 
pdf("Figs/MGBB.CHIP_Risk.forest_plot_CHIP.noANC.pdf", width = 12, height = 8)
forest(x = m, 
       common=F, 
       random=F, 
       hetstat=F, 
       subgroup=k.w>=1, 
       weight.study="same",  
       level=0.95, 
       xlim=c(0.5, 3), 
       smlab="Effect of Exposures\non CHIP\n", 
       smlab.pos=0, 
       colgap=unit(7, "mm"),
       xlab="Odds Ratio (95% CI)", 
       squaresize=0.6, 
       col.subgroup="black", 
       colgap.left=unit(0.1,"cm"),
       colgap.forest.left="3mm", 
       colgap.forest.right="2mm", 
       leftcols=c("studlab"), 
       leftlabs = c("                     "),
       rightcols=c("OR","CI95","Pval", "N", "Cases", "Controls"),
       rightlabs=c("OR","95% CI", "P", "N", "Cases", "Controls"),
       col.inside="black", 
       plotwidth=unit(6.5, "cm"), 
       print.subgroup.name=F)

dev.off()

## Ancestry assoc
##
Anc.final_results <- final_results %>% 
  filter((Exposure %in% c("AFR", "AMR", "SAS", "EAS")) & abs(SE)<=10  )

m <- metagen(TE = Beta,
             lower = lSE,
             upper = uSE,
             studlab = Outcome,
             subgroup = Exposure,
             data = Anc.final_results,
             sm="OR")

m$N <- Anc.final_results$N
m$Cases <- Anc.final_results$Cases
m$Controls <- Anc.final_results$Controls


## 
pdf("Figs/MGBB.CHIP_Risk.forest_plot_CHIP.ANC.pdf", width = 12, height = 12)
forest(x = m, 
       common=F, 
       random=F, 
       hetstat=F, 
       subgroup=k.w>=1, 
       weight.study="same",  
       level=0.95, 
       xlim=c(0.5, 3), 
       smlab="Effect of Exposures\non CHIP\n", 
       smlab.pos=0, 
       colgap=unit(7, "mm"),
       xlab="Odds Ratio (95% CI)", 
       squaresize=0.6, 
       col.subgroup="black", 
       colgap.left=unit(0.1,"cm"),
       colgap.forest.left="3mm", 
       colgap.forest.right="2mm", 
       leftcols=c("studlab"), 
       leftlabs = c("                     "),
       rightcols=c("OR","CI95","Pval", "N", "Cases", "Controls"),
       rightlabs=c("OR","95% CI", "P", "N", "Cases", "Controls"),
       col.inside="black", 
       plotwidth=unit(6.5, "cm"), 
       print.subgroup.name=F)

dev.off()

#------------------- TOPMed
# TOPMed ####
#-------- 
final_results <- fread("Figs/TOPMed_CHIP.risk_factors.results_flag_vars.csv", header=T)
head(final_results)

# filter extremes 
final_results <- final_results %>% filter(abs(SE)<=10)

final_results$Pval <- formatC(x = final_results$P, digits = 1, format = "E")
final_results$OR <- formatC(x = round(exp(final_results$Beta),2), digits = 2, format = "f")
final_results$lSE <- final_results$Beta - 1.96*final_results$SE
final_results$uSE <- final_results$Beta + 1.96*final_results$SE
final_results$CI95 <- paste0("[", formatC(x = round(exp(final_results$lSE),2), digits = 2, format = "f"),
                             "-",
                             formatC(x = round(exp(final_results$uSE),2), digits = 2, format = "f"),
                             "]")
head(final_results)

##
final_results$Outcome <- ordered(final_results$Outcome, 
                                 levels=c("Overall CHIP",
                                          "Expanded CHIP", 
                                          "DNMT3A","TET2",
                                          "ASXL1","Splicing Factors",
                                          "DNA damage"))


## w/o ancestry
nonAnc.final_results <- final_results %>% 
  filter(!(Exposure %in% c("AFR", "AMR", "SAS", "EAS")) )
m <- metagen(TE = Beta,
             lower = lSE,
             upper = uSE,
             studlab = Outcome,
             subgroup = Exposure,
             data = nonAnc.final_results,
             sm="OR")

m$N <- nonAnc.final_results$N
m$Cases <- nonAnc.final_results$Cases
m$Controls <- nonAnc.final_results$Controls
## 
pdf("Figs/TOPMed.CHIP_Risk.forest_plot_CHIP.noANC.pdf", width = 12, height = 8)
forest(x = m, 
       common=F, 
       random=F, 
       hetstat=F, 
       subgroup=k.w>=1, 
       weight.study="same",  
       level=0.95, 
       xlim=c(0.5, 3), 
       smlab="Effect of Exposures\non CHIP\n", 
       smlab.pos=0, 
       colgap=unit(7, "mm"),
       xlab="Odds Ratio (95% CI)", 
       squaresize=0.6, 
       col.subgroup="black", 
       colgap.left=unit(0.1,"cm"),
       colgap.forest.left="3mm", 
       colgap.forest.right="2mm", 
       leftcols=c("studlab"), 
       leftlabs = c("                     "),
       rightcols=c("OR","CI95","Pval", "N", "Cases", "Controls"),
       rightlabs=c("OR","95% CI", "P", "N", "Cases", "Controls"),
       col.inside="black", 
       plotwidth=unit(6.5, "cm"), 
       print.subgroup.name=F)

dev.off()

## Ancestry assoc
##
Anc.final_results <- final_results %>% 
  filter((Exposure %in% c("AFR", "AMR", "SAS", "EAS")) & abs(SE)<=10  )

m <- metagen(TE = Beta,
             lower = lSE,
             upper = uSE,
             studlab = Outcome,
             subgroup = Exposure,
             data = Anc.final_results,
             sm="OR")

m$N <- Anc.final_results$N
m$Cases <- Anc.final_results$Cases
m$Controls <- Anc.final_results$Controls


## 
pdf("Figs/TOPMed.CHIP_Risk.forest_plot_CHIP.ANC.pdf", width = 12, height = 12)
forest(x = m, 
       common=F, 
       random=F, 
       hetstat=F, 
       subgroup=k.w>=1, 
       weight.study="same",  
       level=0.95, 
       xlim=c(0.5, 3), 
       smlab="Effect of Exposures\non CHIP\n", 
       smlab.pos=0, 
       colgap=unit(7, "mm"),
       xlab="Odds Ratio (95% CI)", 
       squaresize=0.6, 
       col.subgroup="black", 
       colgap.left=unit(0.1,"cm"),
       colgap.forest.left="3mm", 
       colgap.forest.right="2mm", 
       leftcols=c("studlab"), 
       leftlabs = c("                     "),
       rightcols=c("OR","CI95","Pval", "N", "Cases", "Controls"),
       rightlabs=c("OR","95% CI", "P", "N", "Cases", "Controls"),
       col.inside="black", 
       plotwidth=unit(6.5, "cm"), 
       print.subgroup.name=F)

dev.off()

#------------------------------- UKBB ####
final_results <- fread("Figs/UKBB_CHIP.risk_factors.results_flag_vars.csv", header=T)
head(final_results)

# filter extremes 
final_results <- final_results %>% filter(abs(SE)<=10)

final_results$Pval <- formatC(x = final_results$P, digits = 1, format = "E")
final_results$OR <- formatC(x = round(exp(final_results$Beta),2), digits = 2, format = "f")
final_results$lSE <- final_results$Beta - 1.96*final_results$SE
final_results$uSE <- final_results$Beta + 1.96*final_results$SE
final_results$CI95 <- paste0("[", formatC(x = round(exp(final_results$lSE),2), digits = 2, format = "f"),
                             "-",
                             formatC(x = round(exp(final_results$uSE),2), digits = 2, format = "f"),
                             "]")
head(final_results)

##
final_results$Outcome <- ordered(final_results$Outcome, 
                                 levels=c("Overall CHIP",
                                          "Expanded CHIP", 
                                          "DNMT3A","TET2",
                                          "ASXL1","Splicing Factors",
                                          "DNA damage"))


## w/o ancestry
nonAnc.final_results <- final_results %>% 
  filter(!(Exposure %in% c("AFR", "AMR", "SAS", "EAS")) )
m <- metagen(TE = Beta,
             lower = lSE,
             upper = uSE,
             studlab = Outcome,
             subgroup = Exposure,
             data = nonAnc.final_results,
             sm="OR")

m$N <- nonAnc.final_results$N
m$Cases <- nonAnc.final_results$Cases
m$Controls <- nonAnc.final_results$Controls
## 
pdf("Figs/UKBB.CHIP_Risk.forest_plot_CHIP.noANC.pdf", width = 12, height = 8)
forest(x = m, 
       common=F, 
       random=F, 
       hetstat=F, 
       subgroup=k.w>=1, 
       weight.study="same",  
       level=0.95, 
       xlim=c(0.5, 3), 
       smlab="Effect of Exposures\non CHIP\n", 
       smlab.pos=0, 
       colgap=unit(7, "mm"),
       xlab="Odds Ratio (95% CI)", 
       squaresize=0.6, 
       col.subgroup="black", 
       colgap.left=unit(0.1,"cm"),
       colgap.forest.left="3mm", 
       colgap.forest.right="2mm", 
       leftcols=c("studlab"), 
       leftlabs = c("                     "),
       rightcols=c("OR","CI95","Pval", "N", "Cases", "Controls"),
       rightlabs=c("OR","95% CI", "P", "N", "Cases", "Controls"),
       col.inside="black", 
       plotwidth=unit(6.5, "cm"), 
       print.subgroup.name=F)

dev.off()

## Ancestry assoc
##
Anc.final_results <- final_results %>% 
  filter((Exposure %in% c("AFR", "AMR", "SAS", "EAS")) & abs(SE)<=10  )

m <- metagen(TE = Beta,
             lower = lSE,
             upper = uSE,
             studlab = Outcome,
             subgroup = Exposure,
             data = Anc.final_results,
             sm="OR")

m$N <- Anc.final_results$N
m$Cases <- Anc.final_results$Cases
m$Controls <- Anc.final_results$Controls


## 
pdf("Figs/UKBB.CHIP_Risk.forest_plot_CHIP.ANC.pdf", width = 12, height = 12)
forest(x = m, 
       common=F, 
       random=F, 
       hetstat=F, 
       subgroup=k.w>=1, 
       weight.study="same",  
       level=0.95, 
       xlim=c(0.5, 3), 
       smlab="Effect of Exposures\non CHIP\n", 
       smlab.pos=0, 
       colgap=unit(7, "mm"),
       xlab="Odds Ratio (95% CI)", 
       squaresize=0.6, 
       col.subgroup="black", 
       colgap.left=unit(0.1,"cm"),
       colgap.forest.left="3mm", 
       colgap.forest.right="2mm", 
       leftcols=c("studlab"), 
       leftlabs = c("                     "),
       rightcols=c("OR","CI95","Pval", "N", "Cases", "Controls"),
       rightlabs=c("OR","95% CI", "P", "N", "Cases", "Controls"),
       col.inside="black", 
       plotwidth=unit(6.5, "cm"), 
       print.subgroup.name=F)

dev.off()

#-------------- All combined ####
d_all <- fread("Figs/Risk_factors_for_CHIP.csv", header = T)
d_all <- d_all %>% filter(abs(SE)<=10 & Exposure!="Other")
d_all$Exposure_v2 <- ifelse(d_all$Exposure=="Smoked >=100 cigarettes",
                            "Ever Smoker",
                            d_all$Exposure)
table(d_all$Exposure_v2)

d_all$Pval <- formatC(x = d_all$P, digits = 1, format = "E")
d_all$OR <- formatC(x = round(exp(d_all$Beta),2), digits = 2, format = "f")
d_all$lSE <- d_all$Beta - 1.96*d_all$SE
d_all$uSE <- d_all$Beta + 1.96*d_all$SE
d_all$CI95 <- paste0("[", formatC(x = round(exp(d_all$lSE),2), digits = 2, format = "f"),
                             "-",
                             formatC(x = round(exp(d_all$uSE),2), digits = 2, format = "f"),
                             "]")
head(d_all)
d_all$Outcome_Exposure <- paste0(d_all$Outcome, " ~ ", d_all$Exposure)
##
d_all$Outcome <- ordered(d_all$Outcome, 
                                 levels=c("Overall CHIP",
                                          "Expanded CHIP", 
                                          "DNMT3A","TET2",
                                          "ASXL1",
                                          "Splicing Factors",
                                          "DNA damage"))

# flag_vars <- unique(d_all$Outcome)
# 
# for(var in flag_vars){
#   ## 
#   cat(var,"\n")
# d <- d_all %>% 
#   filter(abs(SE)<=10 &
#            Outcome==var &
#            Exposure_v2%in%c("Age", 
#                          "Male Sex", 
#                          "Ever Smoker") )
# 
# m <- metagen(TE = Beta,
#                  lower = lSE,
#                  upper = uSE,
#                  #title = "Age",
#                  studlab = Cohort,
#                  subgroup = Exposure_v2,
#                  data = d,
#                  sm="OR", )
# 
# m$N <- d$N
# m$Cases <- d$Cases
# m$Controls <- d$Controls
# ## 
# forest_header <- paste0("Effect of Exposures\non ", var,"\n")
# pdf_string <- paste0("Figs/Age_Sex_SMK.",var,"_Risk.forest_plot.pdf") 
# 
# pdf(pdf_string, width = 12, height = 12)
# forest(x = m, 
#        common=F, 
#        random=F, 
#        hetstat=F, 
#        subgroup=k.w>=1, 
#        weight.study="same",  
#        level=0.95, 
#        xlim=c(0.5, 3), 
#        smlab=forest_header, 
#        smlab.pos=0, 
#        colgap=unit(7, "mm"),
#        xlab="Odds Ratio (95% CI)", 
#        squaresize=0.6, 
#        col.subgroup="black", 
#        colgap.left=unit(0.1,"cm"),
#        colgap.forest.left="3mm", 
#        colgap.forest.right="2mm", 
#        leftcols=c("studlab"), 
#        leftlabs = c("                     "),
#        rightcols=c("OR","CI95","Pval", "N"),
#        rightlabs=c("OR","95% CI", "P", "N"),
#        col.inside="black", 
#        plotwidth=unit(6.5, "cm"), 
#        print.subgroup.name=F)
# 
# dev.off()
# }

### By Exposures Age|Sex| Smk ####
exposures <- c("Age", "Male Sex", "Ever Smoker")
for(exposure in exposures){
#for(var in flag_vars){
  ## 
  cat(exposure,"\n")
  d <- d_all %>% 
    filter(abs(SE)<=10 &
            # Outcome==var &
             Exposure_v2==exposure 
           )
  m <- metagen(TE = Beta,
               lower = lSE,
               upper = uSE,
               #title = "Age",
               studlab = Cohort,
               subgroup = Outcome,
               data = d,
               sm="OR")
  
  m$N <- d$N
  m$Cases <- d$Cases
  m$Controls <- d$Controls
  ## 
  forest_header <- paste0("Effect of ",exposure,"\non CHIP categories\n")
  pdf_string <- paste0("Figs/effect_of_",exposure,".CHIPtraits_Risk.forest_plot.pdf") 
  
  pdf(pdf_string, width = 10, height = 15)
  forest(x = m, 
         common=TRUE, 
         random=FALSE, 
         hetstat=TRUE, 
         subgroup=TRUE, 
         weight.study="same",  
         level=0.95, 
         xlim=c(0.5, 3), 
         smlab=forest_header, 
         smlab.pos=0, 
         colgap=unit(7, "mm"),
         xlab="Odds Ratio (95% CI)", 
         #squaresize=0.6, 
         col.subgroup="black", 
         colgap.left=unit(0.1,"cm"),
         colgap.forest.left="3mm", 
         colgap.forest.right="2mm", 
         leftcols=c("studlab"), 
         leftlabs = c("                     "),
         rightcols=c("OR","CI95","Pval", "N"),
         rightlabs=c("OR","95% CI", "P", "N"),
         col.inside="black",
         col.diamond = "blue", 
         col.diamond.lines = "blue", 
         col.square = "darkgray", 
         col.study = "black",
         col.label.right = "black",
         col.label.left = "black",
         squaresize = 0.9,
         plotwidth=unit(8, "cm"), 
         print.subgroup.name=F, 
         fontsize = 12, 
         spacing = 1,
         addrow.overall=TRUE,
         addrow.subgroups=TRUE,
         addrows.below.overall=10,
         subgroup.hetstat = TRUE,
         overall = FALSE, 
         overall.hetstat = FALSE, 
         test.subgroup.common = FALSE,
         text.common.w = "Overall",
         digits.pval.Q = 1,
         digits.tau2 = 1,
         digits.I2 = 0)
  
  dev.off()
#}
}

## Ancestry ####
flag_vars <- unique(d_all$Outcome)

for(var in flag_vars){
  ## 
  cat(var,"\n")
  d <- d_all %>% 
    filter(abs(SE)<=10 &
             Outcome==var &
             Exposure_v2%in%c("AFR", 
                              "AMR", 
                              "EAS", 
                              "SAS") )

  m <- metagen(TE = Beta,
               lower = lSE,
               upper = uSE,
               studlab = Cohort,
               subgroup = Exposure_v2,
               data = d,
               sm="OR")
  
  m$N <- d$N
  m$Cases <- d$Cases
  m$Controls <- d$Controls
  ## 
  forest_header <- paste0("Effect of Genetic Ancestry (Ref: EUR)\non ",var,"\n")
  pdf_string <- paste0("Figs/noN.effect_of_genAnc.",var,".Risk.forest_plot.pdf") 
  
  pdf(pdf_string, width = 10, height = 10)
p4 <-  forest(x = m, 
         common=TRUE, 
         random=FALSE, 
         hetstat=TRUE, 
         subgroup=TRUE, 
         weight.study="same",  
         level=0.95, 
         xlim=c(0.5, 3), 
         smlab=forest_header, 
         smlab.pos=0, 
         colgap=unit(7, "mm"),
         xlab="Odds Ratio (95% CI)", 
         #squaresize=0.6, 
         col.subgroup="black", 
         colgap.left=unit(0.1,"cm"),
         colgap.forest.left="3mm", 
         colgap.forest.right="2mm", 
         leftcols=c("studlab"), 
         leftlabs = c("                     "),
         rightcols=c("OR","CI95","Pval"),
         rightlabs=c("OR","95% CI", "P"),
         col.inside="black",
         col.diamond = "blue", 
         col.diamond.lines = "blue", 
         col.square = "darkgray", 
         col.study = "black",
         col.label.right = "black",
         col.label.left = "black",
         squaresize = 0.9,
         plotwidth=unit(8, "cm"), 
         print.subgroup.name=F, 
         fontsize = 12, 
         spacing = 1,
         addrow.overall=TRUE,
         addrow.subgroups=TRUE,
         addrows.below.overall=3,
         subgroup.hetstat = TRUE,
         overall = TRUE, 
         overall.hetstat = TRUE, 
         test.subgroup.common = TRUE,
         text.common.w = "Overall", 
         # details=TRUE, # A logical specifying whether details on statistical methods should be printed.
         digits.pval.Q = 1,
         digits.tau2 = 1,
         digits.I2 = 0)
    dev.off()
  #}
}

library(patchwork) 
library(gridExtra)
library(grid)
# Assuming p1, p2, p3, and p4 are your forest plots generated with the forest() function
exposures <- c("Age", "Male Sex", "Ever Smoker")
  d <- d_all %>% 
    filter(abs(SE)<=10 &
             # Outcome==var &
             Exposure_v2=="Age" 
    )
  m1 <- metagen(TE = Beta,
               lower = lSE,
               upper = uSE,
               studlab = Cohort,
               subgroup = Outcome,
               data = d,
               sm="OR")
  m1$N <- d$N
  
# Save plots as objects
p1 <- grid.grabExpr(forest(x = m1, 
                        common=TRUE, 
                        random=FALSE, 
                        hetstat=TRUE, 
                        subgroup=TRUE, 
                        weight.study="same",  
                        level=0.95, 
                        xlim=c(0.5, 3), 
                        smlab=forest_header, 
                        smlab.pos=0, 
                        colgap=unit(7, "mm"),
                        xlab="Odds Ratio (95% CI)", 
                        #squaresize=0.6, 
                        col.subgroup="black", 
                        colgap.left=unit(0.1,"cm"),
                        colgap.forest.left="3mm", 
                        colgap.forest.right="2mm", 
                        leftcols=c("studlab"), 
                        leftlabs = c("                     "),
                        rightcols=c("OR","CI95","Pval", "N"),
                        rightlabs=c("OR","95% CI", "P", "N"),
                        col.inside="black",
                        col.diamond = "blue", 
                        col.diamond.lines = "blue", 
                        col.square = "darkgray", 
                        col.study = "black",
                        col.label.right = "black",
                        col.label.left = "black",
                        squaresize = 0.9,
                        plotwidth=unit(8, "cm"), 
                        print.subgroup.name=F, 
                        fontsize = 12, 
                        spacing = 1,
                        addrow.overall=TRUE,
                        addrow.subgroups=TRUE,
                        addrows.below.overall=3,
                        subgroup.hetstat = TRUE,
                        overall = TRUE, 
                        overall.hetstat = TRUE, 
                        test.subgroup.common = TRUE,
                        text.common.w = "Overall", 
                        # details=TRUE, # A logical specifying whether details on statistical methods should be printed.
                        digits.pval.Q = 1,
                        digits.tau2 = 1,
                        digits.I2 = 0))
# 
d <- d_all %>% 
  filter(abs(SE)<=10 &
           Exposure_v2=="Male Sex")
m2 <- metagen(TE = Beta,
              lower = lSE,
              upper = uSE,
              studlab = Cohort,
              subgroup = Outcome,
              data = d,
              sm="OR")
m2$N <- d$N

p2 <- grid.grabExpr(forest(x = m2, 
                        common=TRUE, 
                        random=FALSE, 
                        hetstat=TRUE, 
                        subgroup=TRUE, 
                        weight.study="same",  
                        level=0.95, 
                        xlim=c(0.5, 3), 
                        smlab=forest_header, 
                        smlab.pos=0, 
                        colgap=unit(7, "mm"),
                        xlab="Odds Ratio (95% CI)", 
                        #squaresize=0.6, 
                        col.subgroup="black", 
                        colgap.left=unit(0.1,"cm"),
                        colgap.forest.left="3mm", 
                        colgap.forest.right="2mm", 
                        leftcols=c("studlab"), 
                        leftlabs = c("                     "),
                        rightcols=c("OR","CI95","Pval", "N"),
                        rightlabs=c("OR","95% CI", "P", "N"),
                        col.inside="black",
                        col.diamond = "blue", 
                        col.diamond.lines = "blue", 
                        col.square = "darkgray", 
                        col.study = "black",
                        col.label.right = "black",
                        col.label.left = "black",
                        squaresize = 0.9,
                        plotwidth=unit(8, "cm"), 
                        print.subgroup.name=F, 
                        fontsize = 12, 
                        spacing = 1,
                        addrow.overall=TRUE,
                        addrow.subgroups=TRUE,
                        addrows.below.overall=3,
                        subgroup.hetstat = TRUE,
                        overall = TRUE, 
                        overall.hetstat = TRUE, 
                        test.subgroup.common = TRUE,
                        text.common.w = "Overall", 
                        # details=TRUE, # A logical specifying whether details on statistical methods should be printed.
                        digits.pval.Q = 1,
                        digits.tau2 = 1,
                        digits.I2 = 0))
#
d <- d_all %>% 
  filter(abs(SE)<=10 &
           Exposure_v2=="Ever Smoker")
m3 <- metagen(TE = Beta,
              lower = lSE,
              upper = uSE,
              studlab = Cohort,
              subgroup = Outcome,
              data = d,
              sm="OR")
m3$N <- d$N

p3 <- grid.grabExpr(forest(x = m3, 
                        common=TRUE, 
                        random=FALSE, 
                        hetstat=TRUE, 
                        subgroup=TRUE, 
                        weight.study="same",  
                        level=0.95, 
                        xlim=c(0.5, 3), 
                        smlab=forest_header, 
                        smlab.pos=0, 
                        colgap=unit(7, "mm"),
                        xlab="Odds Ratio (95% CI)", 
                        #squaresize=0.6, 
                        col.subgroup="black", 
                        colgap.left=unit(0.1,"cm"),
                        colgap.forest.left="3mm", 
                        colgap.forest.right="2mm", 
                        leftcols=c("studlab"), 
                        leftlabs = c("                     "),
                        rightcols=c("OR","CI95","Pval", "N"),
                        rightlabs=c("OR","95% CI", "P", "N"),
                        col.inside="black",
                        col.diamond = "blue", 
                        col.diamond.lines = "blue", 
                        col.square = "darkgray", 
                        col.study = "black",
                        col.label.right = "black",
                        col.label.left = "black",
                        squaresize = 0.9,
                        plotwidth=unit(8, "cm"), 
                        print.subgroup.name=F, 
                        fontsize = 12, 
                        spacing = 1,
                        addrow.overall=TRUE,
                        addrow.subgroups=TRUE,
                        addrows.below.overall=3,
                        subgroup.hetstat = TRUE,
                        overall = TRUE, 
                        overall.hetstat = TRUE, 
                        test.subgroup.common = TRUE,
                        text.common.w = "Overall", 
                        # details=TRUE, # A logical specifying whether details on statistical methods should be printed.
                        digits.pval.Q = 1,
                        digits.tau2 = 1,
                        digits.I2 = 0))
# m4
d <- d_all %>% 
  filter(abs(SE)<=10 &
           Outcome=="Overall CHIP" &
           Exposure_v2%in%c("AFR", 
                            "AMR", 
                            "EAS", 
                            "SAS") )

m4 <- metagen(TE = Beta,
              lower = lSE,
              upper = uSE,
              studlab = Cohort,
              subgroup = Outcome,
              data = d,
              sm="OR")
m4$N <- d$N

p4 <- grid.grabExpr(forest(x = m4, 
                        common=TRUE, 
                        random=FALSE, 
                        hetstat=TRUE, 
                        subgroup=TRUE, 
                        weight.study="same",  
                        level=0.95, 
                        xlim=c(0.5, 3), 
                        smlab=forest_header, 
                        smlab.pos=0, 
                        colgap=unit(7, "mm"),
                        xlab="Odds Ratio (95% CI)", 
                        #squaresize=0.6, 
                        col.subgroup="black", 
                        colgap.left=unit(0.1,"cm"),
                        colgap.forest.left="3mm", 
                        colgap.forest.right="2mm", 
                        leftcols=c("studlab"), 
                        leftlabs = c("                     "),
                        rightcols=c("OR","CI95","Pval"),
                        rightlabs=c("OR","95% CI", "P"),
                        col.inside="black",
                        col.diamond = "blue", 
                        col.diamond.lines = "blue", 
                        col.square = "darkgray", 
                        col.study = "black",
                        col.label.right = "black",
                        col.label.left = "black",
                        squaresize = 0.9,
                        plotwidth=unit(8, "cm"), 
                        print.subgroup.name=F, 
                        fontsize = 12, 
                        spacing = 1,
                        addrow.overall=TRUE,
                        addrow.subgroups=TRUE,
                        addrows.below.overall=3,
                        subgroup.hetstat = TRUE,
                        overall = TRUE, 
                        overall.hetstat = TRUE, 
                        test.subgroup.common = TRUE,
                        text.common.w = "Overall", 
                        # details=TRUE, # A logical specifying whether details on statistical methods should be printed.
                        digits.pval.Q = 1,
                        digits.tau2 = 1,
                        digits.I2 = 0))

# Combine the plots in a 2x2 grid and save as PDF
# pdf("multiple_forest_plots.pdf", 
#     width = 18, height = 20)
# grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
# dev.off()

library(ggpubr)
# Use ggpubr::as_ggplot to convert grobs to ggplot objects
pdf("multiple_forest_plots.cow.pdf", width = 18, height = 30)
plot_grid(as_ggplot(p1), 
          as_ggplot(p2), 
          as_ggplot(p3), 
          as_ggplot(p4), 
          labels = c("a)", 
                     "b)", 
                     "c)", 
                     "d)"), 
          label_size = 14, 
          ncol = 2, nrow = 2)
dev.off()


### Save the Combined Plot as PNG and PDF ###

# Save as PNG
# ggsave("~/FigS1.AoU.multipanel_plot.png", multi_plot, width = 16, height = 10, dpi = 300)
