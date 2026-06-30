## old file in : /medpop/esp/mesbah/GWAS_CHIP/Latest_Meta_Jun2021/meta_gwas/ukb200k/Meta_CHIP.Jun27_2021/jun27_2021.work_log.sh

	# Random Metal
git clone https://github.com/explodecomputer/random-metal.git
	# Get TopMEd gwas
gsutil -m cp gs://ukbb_v2/projects/muddin/META/inputTOPMed/*_2020_saige_results.tsv.gz topmed/
	# Get UKB200k 
gsutil -m cp gs://ukbb_v2/projects/muddin/chip_gwas/UKB200k_23Jun2021/additive/Chr*.regenie.gz ukb200k/
	# get MGB GSA data
gsutil -m cp gs://ukbb_v2/projects/muddin/chip_gwas/mgbb/GSA/CHIP_Feb28_2021/additive/Chr*_additive_has{CHIP,VAF10,DNMT3A,TET2,ASXL1}.regenie.gz mgb/

# Header
echo -e 'CHIP\nExtendedCHIP\nDNMT3A\nTET2\nASXL1' > ../Pheno_list.list

while read phenos; do echo -e "Name\tChr\tPos\tRef\tAlt\tTrait\tCohort\tModel\tEffect\tLCI_Effect\tUCI_Effect\tPval\tAAF\tNum_Cases\tCases_Ref\tCases_Het\tCases_Alt\tNum_Controls\tControls_Ref\tControls_Het\tControls_Alt\tREGENIE_BETA\tREGENIE_SE\tINFO\tMAC" > Chr1_22.has${phenos}.ukb200k.tsv; done <../Pheno_list.list


while read phenos; do for chr in {1..22}; do zcat Chr${chr}_additive_has${phenos}.regenie.gz | sed '1d' | tr ';' '\t' | sed -e 's:REGENIE_BETA=::g' -e 's:REGENIE_SE=::g' -e 's:INFO=::g' -e 's:MAC=::g' >> Chr1_22.has${phenos}.ukb200k.tsv; done; done <../Pheno_list.list


### August 15, 2021
# gsutil ls gs://ukbb_v2/projects/muddin/chip_gwas/UKB200k_11Aug2021/additive/Chr1_additive_has*.regenie.gz
gsutil -m cp gs://ukbb_v2/projects/muddin/chip_gwas/UKB200k_11Aug2021/additive/Chr{1..22}_additive_has*.regenie.gz .

echo -e 'CHIP\nExpandedCHIP\nDNMT3A\nTET2\nASXL1\nJAK2\nexpASXL1\nexpDNMT3A\nexpTET2' > Aug11.Pheno_list.list

while read phenos; do echo -e "Name\tChr\tPos_hg37\tRef\tAlt\tTrait\tCohort\tModel\tEffect\tLCI_Effect\tUCI_Effect\tPval\tAAF\tNum_Cases\tCases_Ref\tCases_Het\tCases_Alt\tNum_Controls\tControls_Ref\tControls_Het\tControls_Alt\tREGENIE_BETA\tREGENIE_SE\tINFO\tMAC" | gzip -c > chr1_22.has${phenos}.11Aug2021_ukb200k.tsv.gz; done <Aug11.Pheno_list.list


while read phenos; do for chr in {1..22}; do zcat Chr${chr}_additive_has${phenos}.regenie.gz | sed '1d' | tr ';' '\t' | sed -e 's:REGENIE_BETA=::g' -e 's:REGENIE_SE=::g' -e 's:INFO=::g' -e 's:MAC=::g' | gzip -c >> chr1_22.has${phenos}.11Aug2021_ukb200k.tsv.gz; done; done <Aug11.Pheno_list.list

## Convert hg19 to hg38
while read phenotypes; do qsub -R y -pe smp 2 -binding linear:2 -l h_rt=08:00:00 -l h_vmem=30G -wd /broad/hptmp/mesbah/ukb_chip/meta_gwas/ukb200k/aug11_2021_hg38/tmpdir -N prep.ukb200k.${phenotypes} /broad/hptmp/mesbah/ukb_chip/meta_gwas/ukb200k/run_ukb_lifted_2_GRCh38.11Aug2021.sh /broad/hptmp/mesbah/ukb_chip/meta_gwas/ukb200k/ukb_lifted_2_GRCh38.11Aug2021.R /broad/hptmp/mesbah/ukb_chip/meta_gwas/ukb200k/aug11_2021/chr1_22.has${phenotypes}.11Aug2021_ukb200k.tsv.gz /medpop/esp/mesbah/GWAS_CHIP/Latest_Meta_Jun2021/meta_gwas/lift_over/c2.varIDs_38.varIDs_37.tsv.gz 2 /broad/hptmp/mesbah/ukb_chip/meta_gwas/ukb200k/aug11_2021_hg38/hg38_chr1_22.has${phenotypes}.11Aug2021_ukb200k.tsv.gz; done </broad/hptmp/mesbah/ukb_chip/meta_gwas/ukb200k/aug11_2021/Aug11.Pheno_list.list

########################### 28 Oct, 2021
########################### 
######### ukb150k vs ukb50k GWAS
use Google-Cloud-SDK
use Tabix
gsutil -m cp gs://ukbb_v2/projects/muddin/chip_gwas/ukb50k_vs_150k/UKB50k_Oct2021/additive/Chr*.regenie.gz ukb50k/

gsutil -m cp gs://ukbb_v2/projects/muddin/chip_gwas/ukb50k_vs_150k/UKB150k_Oct2021/additive/Chr*.regenie.gz ukb150k/

# 
echo -e 'CHIP\nExpandedCHIP\nDNMT3A\nTET2' > Oct28.Pheno_list.list


while read phenos; do echo -e "Name\tChr\tPos_hg37\tRef\tAlt\tTrait\tCohort\tModel\tEffect\tLCI_Effect\tUCI_Effect\tPval\tAAF\tNum_Cases\tCases_Ref\tCases_Het\tCases_Alt\tNum_Controls\tControls_Ref\tControls_Het\tControls_Alt\tREGENIE_BETA\tREGENIE_SE\tINFO\tMAC" > ukb50k/chr1_22.has${phenos}.28Oct2021_ukb50k.tsv; done <Oct28.Pheno_list.list

while read phenos; do for files in $(ls -lhv /broad/hptmp/mesbah/ukb_chip/ukb150_vs_ukb50k/ukb50k/Chr*_has${phenos}.regenie.gz | awk '{print $NF}'); do zcat ${files} | sed '1d' | tr ';' '\t' | sed -e 's:REGENIE_BETA=::g' -e 's:REGENIE_SE=::g' -e 's:INFO=::g' -e 's:MAC=::g' >> /broad/hptmp/mesbah/ukb_chip/ukb150_vs_ukb50k/ukb50k/chr1_22.has${phenos}.28Oct2021_ukb50k.tsv; done; done <Oct28.Pheno_list.list &

while read phenos; do echo -e "Name\tChr\tPos_hg37\tRef\tAlt\tTrait\tCohort\tModel\tEffect\tLCI_Effect\tUCI_Effect\tPval\tAAF\tNum_Cases\tCases_Ref\tCases_Het\tCases_Alt\tNum_Controls\tControls_Ref\tControls_Het\tControls_Alt\tREGENIE_BETA\tREGENIE_SE\tINFO\tMAC" > ukb150k/chr1_22.has${phenos}.28Oct2021_ukb150k.tsv; done <Oct28.Pheno_list.list

while read phenos; do for my_files in $(ls -lhv /broad/hptmp/mesbah/ukb_chip/ukb150_vs_ukb50k/ukb150k/Chr*_has${phenos}.regenie.gz | awk '{print $NF}'); do zcat ${my_files} | sed '1d' | tr ';' '\t' | sed -e 's:REGENIE_BETA=::g' -e 's:REGENIE_SE=::g' -e 's:INFO=::g' -e 's:MAC=::g' >> /broad/hptmp/mesbah/ukb_chip/ukb150_vs_ukb50k/ukb150k/chr1_22.has${phenos}.28Oct2021_ukb150k.tsv; done; done <Oct28.Pheno_list.list &

bgzip /broad/hptmp/mesbah/ukb_chip/ukb150_vs_ukb50k/ukb50k/chr1_22.has${phenos}.28Oct2021_ukb50k.tsv
tabix -s 2 -b 3 -e 3 /broad/hptmp/mesbah/ukb_chip/ukb150_vs_ukb50k/ukb50k/chr1_22.has${phenos}.28Oct2021_ukb50k.tsv.gz
 
bgzip /broad/hptmp/mesbah/ukb_chip/ukb150_vs_ukb50k/ukb150k/chr1_22.has${phenos}.28Oct2021_ukb150k.tsv
tabix -s 2 -b 3 -e 3 /broad/hptmp/mesbah/ukb_chip/ukb150_vs_ukb50k/ukb150k/chr1_22.has${phenos}.28Oct2021_ukb150k.tsv.gz


## Topmed 2019 P<= 5e-5
zcat /medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/topmed2019.all_samples.hasCHIP.tsv.gz | head -1 > topmed.chip.p5e5.tsv
zcat /medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/topmed2019.all_samples.hasCHIP.tsv.gz | awk '$12<=5e-5' >> topmed.chip.p5e5.tsv &

zcat /medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/topmed2019.all_samples.hasDNMT3A.tsv.gz | head -1 > topmed.dnmt3a.p5e5.tsv
zcat /medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/topmed2019.all_samples.hasDNMT3A.tsv.gz | awk '$12<=5e-5' >> topmed.dnmt3a.p5e5.tsv

zcat /medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/topmed2019.all_samples.hasTET2.tsv.gz | head -1 > topmed.tet2.p5e5.tsv
zcat /medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/topmed2019.all_samples.hasTET2.tsv.gz | awk '$12<=5e-5' >> topmed.tet2.p5e5.tsv &

## 
use R-4.0
R
data.table::setDTthreads(2)
library(data.table)
# Lifover <- fread("/medpop/esp/mesbah/GWAS_CHIP/inputUKB/UKB_liftedGRCh38.txt.gz", header=T, fill=TRUE) 
Lifover <- fread("/medpop/esp2/mesbah/variantIDs/hg38_rsid_hg19_liftedukb.txt.gz", header=T, fill=TRUE) 
names(Lifover) <- c("SNPid_38", "RSID", "ID37")
Lifover.nodup <- subset(Lifover, !duplicated(Lifover$RSID))
rm(Lifover)

topmed_chip <- fread("/broad/hptmp/mesbah/ukb_chip/ukb150_vs_ukb50k/topmed.chip.p5e5.tsv", header=T)
topmed_chip$SNPid_38 <- paste(topmed_chip$CHR, topmed_chip$POS, topmed_chip$REF,topmed_chip$ALT, sep=":")
topmed_chip <- merge(topmed_chip, Lifover.nodup, by="SNPid_38")

topmed_dnmt3a <- fread("/broad/hptmp/mesbah/ukb_chip/ukb150_vs_ukb50k/topmed.dnmt3a.p5e5.tsv", header=T)
topmed_dnmt3a$SNPid_38 <- paste(topmed_dnmt3a$CHR, topmed_dnmt3a$POS, topmed_dnmt3a$REF, topmed_dnmt3a$ALT, sep=":")
topmed_dnmt3a <- merge(topmed_dnmt3a, Lifover.nodup, by="SNPid_38")

topmed_tet2 <- fread("/broad/hptmp/mesbah/ukb_chip/ukb150_vs_ukb50k/topmed.tet2.p5e5.tsv", header=T)
topmed_tet2$SNPid_38 <- paste(topmed_tet2$CHR, topmed_tet2$POS, topmed_tet2$REF, topmed_tet2$ALT, sep=":")
topmed_tet2 <- merge(topmed_tet2, Lifover.nodup, by="SNPid_38")

### load UKB
Ukb50k.chip <- fread("/broad/hptmp/mesbah/ukb_chip/ukb150_vs_ukb50k/ukb50k/chr1_22.hasCHIP.28Oct2021_ukb50k.tsv", header=T)
Ukb50k.chip <- subset(Ukb50k.chip, !duplicated(Ukb50k.chip$Name) )
topmed.ukb50k_chip <- merge(Ukb50k.chip, topmed_chip, by.x="Name", by.y="RSID")
r_chip50k <- cor(topmed.ukb50k_chip$BETA/topmed.ukb50k_chip$SE, topmed.ukb50k_chip$REGENIE_BETA/topmed.ukb50k_chip$REGENIE_SE)

Ukb150k.chip <- fread("/broad/hptmp/mesbah/ukb_chip/ukb150_vs_ukb50k/ukb150k/chr1_22.hasCHIP.28Oct2021_ukb150k.tsv", header=T)
Ukb150k.chip <- subset(Ukb150k.chip, !duplicated(Ukb150k.chip$Name) )
topmed.ukb150k_chip <- merge(Ukb150k.chip, topmed_chip, by.x="Name", by.y="RSID")
r_chip150k <- cor(topmed.ukb150k_chip$BETA/topmed.ukb150k_chip$SE, topmed.ukb150k_chip$REGENIE_BETA/topmed.ukb150k_chip$REGENIE_SE)

	# plot:
png("/broad/hptmp/mesbah/ukb_chip/ukb150_vs_ukb50k/chip.ukb50k_vs_150k.png", width=8, height=5, units= "in", res=300, pointsize = 5)
par(mfrow=c(1,2), mar= c(5, 5.5, 4, 5.5))

plot(y=topmed.ukb50k_chip$BETA/topmed.ukb50k_chip$SE, x=topmed.ukb50k_chip$REGENIE_BETA/topmed.ukb50k_chip$REGENIE_SE, ylab="TopMed Z-stat", xlab="UKB50k Z-stat", main=paste0("CHIP: UKB50k, r= ", round(r_chip50k,3) ))

plot(y=topmed.ukb150k_chip$BETA/topmed.ukb150k_chip$SE, x=topmed.ukb150k_chip$REGENIE_BETA/topmed.ukb150k_chip$REGENIE_SE, ylab="TopMed Z-stat", xlab="UKB150k Z-stat", main=paste0("UKB150k, r= ", round(r_chip150k,3) ))

dev.off()


	# DNMT3A
Ukb50k.dnmt3a <- fread("/broad/hptmp/mesbah/ukb_chip/ukb150_vs_ukb50k/ukb50k/chr1_22.hasDNMT3A.28Oct2021_ukb50k.tsv", header=T)
Ukb50k.dnmt3a <- subset(Ukb50k.dnmt3a, !duplicated(Ukb50k.dnmt3a $Name) )
topmed.ukb50k_dnmt3a <- merge(Ukb50k.dnmt3a, topmed_dnmt3a, by.x="Name", by.y="RSID")
r_dnmt3a50k <- cor(topmed.ukb50k_dnmt3a$BETA/topmed.ukb50k_dnmt3a$SE, topmed.ukb50k_dnmt3a$REGENIE_BETA/topmed.ukb50k_dnmt3a$REGENIE_SE)


Ukb150k.dnmt3a <- fread("/broad/hptmp/mesbah/ukb_chip/ukb150_vs_ukb50k/ukb150k/chr1_22.hasDNMT3A.28Oct2021_ukb150k.tsv", header=T)
Ukb150k.dnmt3a <- subset(Ukb150k.dnmt3a, !duplicated(Ukb150k.dnmt3a $Name) )
topmed.ukb150k_dnmt3a <- merge(Ukb150k.dnmt3a, topmed_dnmt3a, by.x="Name", by.y="RSID")
r_dnmt3a150k <- cor(topmed.ukb150k_dnmt3a$BETA/topmed.ukb150k_dnmt3a$SE, topmed.ukb150k_dnmt3a$REGENIE_BETA/topmed.ukb150k_dnmt3a$REGENIE_SE)

	# plot
png("/broad/hptmp/mesbah/ukb_chip/ukb150_vs_ukb50k/dnmt3a.ukb50k_vs_150k.png", width=8, height=5, units= "in", res=300, pointsize = 5)
par(mfrow=c(1,2), mar= c(5, 5.5, 4, 5.5))

plot(x=topmed.ukb50k_dnmt3a$REGENIE_BETA/topmed.ukb50k_dnmt3a$REGENIE_SE, y=topmed.ukb50k_dnmt3a$BETA/topmed.ukb50k_dnmt3a$SE,  ylab="TopMed Z-stat", xlab="UKB50k Z-stat", main=paste0("DNMT3A: UKB50k, r= ", round(r_dnmt3a50k,3) ))

plot(x=topmed.ukb150k_dnmt3a$REGENIE_BETA/topmed.ukb150k_dnmt3a$REGENIE_SE, y= topmed.ukb150k_dnmt3a$BETA/topmed.ukb150k_dnmt3a$SE,  ylab="TopMed Z-stat", xlab="UKB150k Z-stat", main=paste0("UKB150k, r= ", round(r_dnmt3a150k,3) ))
dev.off()

rb_dnmt3a150k <- cor(topmed.ukb150k_dnmt3a$BETA, topmed.ukb150k_dnmt3a$REGENIE_BETA)
rb_dnmt3a50k <- cor(topmed.ukb50k_dnmt3a$BETA, topmed.ukb50k_dnmt3a$REGENIE_BETA)

png("/broad/hptmp/mesbah/ukb_chip/ukb150_vs_ukb50k/dnmt3a_beta.ukb50k_vs_150k.png", width=8, height=5, units= "in", res=300, pointsize = 5)
par(mfrow=c(1,2), mar= c(5, 5.5, 4, 5.5))

plot(x=topmed.ukb50k_dnmt3a$REGENIE_BETA, y=topmed.ukb50k_dnmt3a$BETA,  ylab="TopMed Beta", xlab="UKB50k Beta", main=paste0("DNMT3A: UKB50k, r= ", round(rb_dnmt3a50k,3) ))

plot(x=topmed.ukb150k_dnmt3a$REGENIE_BETA, y= topmed.ukb150k_dnmt3a$BETA,  ylab="TopMed Beta", xlab="UKB150k Beta", main=paste0("UKB150k, r= ", round(rb_dnmt3a150k,3) ))
dev.off()

	# TET2
Ukb50k.tet2 <- fread("/broad/hptmp/mesbah/ukb_chip/ukb150_vs_ukb50k/ukb50k/chr1_22.hasTET2.28Oct2021_ukb50k.tsv", header=T)
Ukb50k.tet2 <- subset(Ukb50k.tet2, !duplicated(Ukb50k.tet2 $Name) )
topmed.ukb50k_tet2 <- merge(Ukb50k.tet2, topmed_tet2, by.x="Name", by.y="RSID")
r_tet250k <- cor(topmed.ukb50k_dnmt3a$BETA/topmed.ukb50k_dnmt3a$SE, topmed.ukb50k_dnmt3a$REGENIE_BETA/topmed.ukb50k_dnmt3a$REGENIE_SE)


Ukb150k.tet2 <- fread("/broad/hptmp/mesbah/ukb_chip/ukb150_vs_ukb50k/ukb150k/chr1_22.hasTET2.28Oct2021_ukb150k.tsv", header=T)
Ukb150k.tet2 <- subset(Ukb150k.tet2, !duplicated(Ukb150k.tet2 $Name) )
topmed.ukb150k_tet2 <- merge(Ukb150k.tet2, topmed_tet2, by.x="Name", by.y="RSID")
r_tet2150k <- cor(topmed.ukb50k_dnmt3a$BETA/topmed.ukb50k_dnmt3a$SE, topmed.ukb50k_dnmt3a$REGENIE_BETA/topmed.ukb50k_dnmt3a$REGENIE_SE)


	# plot
png("/broad/hptmp/mesbah/ukb_chip/ukb150_vs_ukb50k/tet2.ukb50k_vs_150k.png", width=6.5, height=5.5, units= "in", res=300, pointsize = 5)
par(mfrow=c(1,2), mar= c(5, 5.5, 4, 5.5))

plot(x=topmed.ukb50k_dnmt3a$BETA/topmed.ukb50k_dnmt3a$SE, y=topmed.ukb50k_dnmt3a$REGENIE_BETA/topmed.ukb50k_dnmt3a$REGENIE_SE, xlab="TopMed Z-stat", xlab="UKB50k Z-stat", main=paste0("DNMT3A: UKB50k, r= ", r_dnmt3a50k) )

plot(x= topmed.ukb150k_dnmt3a$BETA/topmed.ukb150k_dnmt3a$SE, y=topmed.ukb150k_dnmt3a$REGENIE_BETA/topmed.ukb150k_dnmt3a$REGENIE_SE, xlab="TopMed Z-stat", xlab="UKB150k Z-stat", main=paste0("UKB150k, r= ", r_dnmt3a150k) )
dev.off()

##### 12 Jan 2022
## prep for plot
while read namess; do zcat chr1_22.annot.GWAMA.meta_${namess}.allsamples_topmed_ukbb_mgbb_BioVU.hg38.eaf001_min2Studies.allPval.tsv.gz | awk '{print $6"\t"$1"\t"$2"\t"$12}' > for_plot.chr1_22.annot.GWAMA.meta_${namess}.allsamples_topmed_ukbb_mgbb_BioVU.hg38.eaf001_min2Studies.allPval.tsv; done < <(echo -e "CHIP\nDNMT3A\nTET2") &

## jobs
while read mytraits; do qsub -R y -wd /broad/hptmp/mesbah/ukb_chip.v2/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=20G -l h_rt=10:00:00 -N plot_all_meta4.${mytraits}.12Jan22 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/02.Metal_MetaAnalysis/02.plotGWAS/01.run.plot.metaGWAS.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/02.Metal_MetaAnalysis/02.plotGWAS/01.plot.metaGWAS.R /broad/hptmp/mesbah/ukb_chip.v2/gwas_summary/for_plot.chr1_22.annot.GWAMA.meta_${mytraits}.allsamples_topmed_ukbb_mgbb_BioVU.hg38.eaf001_min2Studies.allPval.tsv /broad/hptmp/mesbah/ukb_chip.v2/gwas_summary ${mytraits}; done < <(echo -e "CHIP\nDNMT3A\nTET2")

