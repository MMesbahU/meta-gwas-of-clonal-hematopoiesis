
######### plot h2 Enrichment
library(data.table)
# http://dougspeed.com/bldldak/
## added: echo -e "65 LDAK_Weightings\n66 Base_Category" >> bldnames 
bldnames <- fread("/Volumes/mesbah/gwas/sumher/ref/bldnames")
names(bldnames) <- c("Sr", "Category")

# chip_enr <- fread("/Volumes/mesbah/gwas/sumher/sumher.pop05_eur05.ldak.lifted_hg37.eur_metaGWAS.CHIP.GWAMA.hg37_dbSNP.enrich")
chip_enr <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/SumHer/EUR_sumher/sumher.pop0.0631_sample0.0639.eur_metaGWAS.CHIP.enrich")
chip_enr$Category <- bldnames$Category
chip_enr$Category <- factor(chip_enr$Category, levels = c(chip_enr$Category))
barplot(chip_enr$Enrichment[1:64], names.arg =  chip_enr$Category[1:64], las=2)
chip_enr$CHIP_Category <- "Overall CHIP"

# tet2_enr <- fread("/Volumes/mesbah/gwas/sumher/sumher.pop05_eur05.ldak.lifted_hg37.eur_metaGWAS.TET2.GWAMA.hg37_dbSNP.enrich")
tet2_enr <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/SumHer/EUR_sumher/sumher.pop0.0131_sample0.0137.eur_metaGWAS.TET2.enrich")
tet2_enr$Category <- bldnames$Category
tet2_enr$Category <- factor(tet2_enr$Category, levels = c(tet2_enr$Category))
tet2_enr$CHIP_Category <- "TET2"
barplot(tet2_enr$Enrichment[1:64], 
        names.arg = tet2_enr$Category[1:64], 
        las=2)

# dnmt3a_enr <- fread("/Volumes/mesbah/gwas/sumher/sumher.pop05_eur05.ldak.lifted_hg37.eur_metaGWAS.DNMT3A.GWAMA.hg37_dbSNP.enrich")
dnmt3a_enr <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/SumHer/EUR_sumher/sumher.pop0.0359_sample0.0349.eur_metaGWAS.DNMT3A.enrich")
dnmt3a_enr$Category <- bldnames$Category
dnmt3a_enr$Category <- factor(dnmt3a_enr$Category, levels = c(dnmt3a_enr$Category))
barplot(dnmt3a_enr$Enrichment[1:64], names.arg =  dnmt3a_enr$Category[1:64], las=2)
dnmt3a_enr$CHIP_Category <- "DNMT3A"

## combine all 
her_enrichment <- as.data.frame(rbind(chip_enr[1:64],
                                      tet2_enr[1:64], 
                                      dnmt3a_enr[1:64]))
names(her_enrichment)[3] <- "Share_SD" 
names(her_enrichment)[6] <- "Enrichment_SD" 
## plot
library(ggplot2)
library(ggpubr)
library(cowplot)
theme_set(theme_cowplot())
library(dplyr)
# Standard deviation of the mean as error bar
# Standard deviation of the mean as error bar
png("~/Documents/Project/CHIP_GWAS/rerun/Figures/fig2d.eur_enrichment.Sept12.png",
    width = 14, height = 8, units = 'in', res = 300)

her_enrichment %>% filter(Enrichment_SD>=0 & 
                            ((Enrichment_SD+Enrichment)>=0 & 
                               (Enrichment-Enrichment_SD)>=0 ) ) %>% 
  ggplot(., aes(x=Category, y=Enrichment, fill=CHIP_Category)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Enrichment - Enrichment_SD, 
                    ymax= Enrichment + Enrichment_SD), 
                width=0.2,
                position=position_dodge(0.9)) + 
  geom_hline(yintercept=1, lty=2, linetype="dashed") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1,  hjust=1),
        legend.title = element_blank(), legend.position = "right") +
  scale_y_continuous(breaks = c(0, 1,seq(from = 2, to = 20, by = 2))) +
  xlab("") + ylab(expression("Estimated Enrichment"))
dev.off()

# her_enrichment %>% filter(Enrichment_SD<=20 & Enrichment_SD>=-20 & 
#                             ((Enrichment_SD+Enrichment)>=1 & (Enrichment-Enrichment_SD)>=1 ) ) %>% 
#   ggplot(., aes(x=Category, y=Enrichment, fill=CHIP_Category)) + 
#   geom_bar(stat="identity", position=position_dodge()) +
#   geom_errorbar(aes(ymin=Enrichment - Enrichment_SD, 
#                     ymax= Enrichment + Enrichment_SD), 
#                 width=0.2,
#                 position=position_dodge(0.9)) + 
#   geom_hline(yintercept=1, lty=2, linetype="dashed") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 1,  hjust=1),
#         legend.title = element_blank(), legend.position = "right") +
#   scale_y_continuous(breaks = c(seq(from = 0, to = 15, by = 1) )) +
#   xlab("") + ylab("Estimated Enrichment")
# 

# her_enrichment %>% filter(Enrichment_SD<=5 & Enrichment_SD>=-5) %>% 
#   ggplot(., aes(x=Category, y=Enrichment, 
#                 ymin=Enrichment - Enrichment_SD, 
#                 ymax= Enrichment + Enrichment_SD, 
#                                        color=CHIP_Category)) +
#   geom_point(show.legend=TRUE, 
#                   position=position_jitter(width=0.5)) + 
#   geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
#   coord_flip() +  # flip coordinates (puts labels on y axis)
#   xlab("Annotation") + ylab("Enrichment") +
#   theme(legend.position="top") +
#   facet_wrap(~Category, scales = "free")

###########