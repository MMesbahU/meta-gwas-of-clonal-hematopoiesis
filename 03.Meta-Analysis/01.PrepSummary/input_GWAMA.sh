## MGBB53k 
while read traits; do echo -e "/medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/topmed2019.all_samples.has${traits}.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/inputUKB/allUKB/lifted_hg38.chr1_22.has${traits}.11Aug2021_ukb200k_allsamples.INFO.tsv.gz\n/broad/hptmp/mesbah/gwas/mgbb53k/mgbb53k_gwama.chr1_22.has${traits}.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/BioVU/BioVU.saige_has${traits}_results_merged_subset.tsv.gz" > /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/01.PrepSummary/gwama_input.has${traits}.TOPMed2019UKB200kMGBB53kBioVU54k.txt; done < <(echo -e "CHIP\nDNMT3A\nTET2\nASXL1")

##hg37:  UKB 250 vs 200k 
while read traits; do echo -e "/broad/hptmp/mesbah/gwas/ukb450k/ukb250vs200k/chr1_22.ukb200k.has${traits}.tsv.gz\n/broad/hptmp/mesbah/gwas/ukb450k/250k.step2/chr1_22.ukb250k.has${traits}.tsv.gz" > /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/01.PrepSummary/gwama_input.has${traits}.UKB200k_vs_UKB250k.hg37.txt; done < <(echo -e "CHIP\nDNMT3A\nTET2\nASXL1")

## hg38: UKB450k, TOPMED, MGBB53k, BioVU
# /broad/hptmp/mesbah/gwas/ukb450k/450k.step2/lifted_hg38.chr1_22.hasASXL1.ukb450k.INFO.tsv.gz
while read traits; do echo -e "/medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/topmed2019.all_samples.has${traits}.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/inputUKB/allUKB/lifted_hg38.chr1_22.has${traits}.11Aug2021_ukb200k_allsamples.INFO.tsv.gz\n/broad/hptmp/mesbah/gwas/mgbb53k/mgbb53k_gwama.chr1_22.has${traits}.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/BioVU/BioVU.saige_has${traits}_results_merged_subset.tsv.gz" > /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/01.PrepSummary/gwama_input.has${traits}.TOPMed2019UKB450kMGBB53kBioVU54k.txt; done < <(echo -e "CHIP\nDNMT3A\nTET2\nASXL1")

## meta UKB450k, TOPMED, MGBB53k, BioVU
while read traits; do echo -e "/medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/topmed2019.all_samples.has${traits}.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/inputUKB/allUKB/lifted_hg38.chr1_22.has${traits}.11Aug2021_ukb200k_allsamples.INFO.tsv.gz\n/broad/hptmp/mesbah/gwas/mgbb53k/mgbb53k_gwama.chr1_22.has${traits}.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/BioVU/BioVU.saige_has${traits}_results_merged_subset.tsv.gz" > /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/01.PrepSummary/gwama_input.has${traits}.TOPMed2019metaUKB450kMGBB53kBioVU54k.txt; done < <(echo -e "CHIP\nDNMT3A\nTET2\nASXL1")


########################### Jan 31, 2023
	
	### Prep AoU Summary
while read traits; do echo -e "SNPID\tREF\tALT\tAAF\tBETA\tSE\tN\tP" > /medpop/esp2/mesbah/datasets/allofus/gwas_ch/${traits}_ch.chr1_22.tsv; for chr in {1..22}; do zcat /medpop/esp2/mesbah/datasets/allofus/gwas_ch/${traits}*step2_chr${chr}_*.regenie.gz | awk 'NR>1{print $3"\t"$4"\t"$5"\t"$6"\t"$10"\t"$11"\t"$8"\t"10^-$13}' | sed 's:\_:\::g' >> /medpop/esp2/mesbah/datasets/allofus/gwas_ch/${traits}_ch.chr1_22.tsv; done; done < <(echo -e "CHIP\nDNMT3A\nTET2") &

gzip /medpop/esp2/mesbah/datasets/allofus/gwas_ch/{CHIP,TET2,DNMT3A}_ch.chr1_22.tsv &

	### TOPMED 2021
while read traits; do zcat /medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2021/topmed_has${traits}*.tsv.gz | awk '(NR==1){print "SNPID\tREF\tALT\tAAF\tBETA\tSE\tN\tP"}(NR>1){print $1"\t"$4"\t"$5"\t"$7"\t"$10"\t"$11"\t"$9"\t"$13}' | gzip -c > /medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2021/topmed2021_has${traits}.tsv.gz; done < <(echo -e "CHIP\nDNMT3A\nTET2") &

#############

