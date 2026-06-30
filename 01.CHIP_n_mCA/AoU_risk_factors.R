##############################
#### Forest plot
##############################
library(data.table) # version 1.14.6
library(meta) # version 6.2-1
library(grid) # version 4.2.2
library(scales) # version 1.2.1
setwd("~/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/NG_v2/")

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
