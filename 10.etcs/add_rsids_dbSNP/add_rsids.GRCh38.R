## Aug 15, 2021: 
## while read phenotypes; do qsub -R y -pe smp 2 -binding linear:2 -l h_rt=08:00:00 -l h_vmem=30G -wd /broad/hptmp/mesbah/ukb_chip/meta_gwas/ukb200k/aug11_2021_hg38/tmpdir -N prep.ukb200k.${phenotypes} /broad/hptmp/mesbah/ukb_chip/meta_gwas/ukb200k/run_ukb_lifted_2_GRCh38.11Aug2021.sh /broad/hptmp/mesbah/ukb_chip/meta_gwas/ukb200k/ukb_lifted_2_GRCh38.11Aug2021.R /broad/hptmp/mesbah/ukb_chip/meta_gwas/ukb200k/aug11_2021/chr1_22.has${phenotypes}.11Aug2021_ukb200k.tsv.gz /medpop/esp/mesbah/GWAS_CHIP/Latest_Meta_Jun2021/meta_gwas/lift_over/c2.varIDs_38.varIDs_37.tsv.gz 2 /broad/hptmp/mesbah/ukb_chip/meta_gwas/ukb200k/aug11_2021_hg38/hg38_chr1_22.has${phenotypes}.11Aug2021_ukb200k.tsv.gz; done </broad/hptmp/mesbah/ukb_chip/meta_gwas/ukb200k/aug11_2021/Aug11.Pheno_list.list
##
###########################
############################


ARGs <- commandArgs(TRUE)

gwas_data <- ARGs[1]

liftOverfile <- ARGs[2]

numThreads <- as.numeric(ARGs[3])

output_file <- ARGs[4]

######
data.table::setDTthreads(numThreads)
library(data.table)
########
# zcat /medpop/esp2/mesbah/Meta_GWAS/Aug2021/eur_chip_maf1.29kTopmed_167kUKB_10kMGB_54kBioVU.tsv.gz | head -2
# SNP	ALLELE1	ALLELE2	FreqAllele1	Effect	StdErr	Pvalue	N
# chr1:10177:A:AC	A	AC	0.601	0.0099338	0.0223773	0.6571	167713

GWAS <- fread(gwas_data, header=T, stringsAsFactors=F, nThread=numThreads, fill=TRUE)
#
gc()
str(GWAS)

## Load rsids
# chr1:POS:REF:ALT     chr1    POS     ID      REF     ALT
liftedData <- fread(liftOverfile, header=T, stringsAsFactors=F, nThread=numThreads, fill=TRUE)
names(liftedData) <- c("SNP", "CHR", "POS", "rsid","REF", "ALT")

gc()
str(liftedData)
##
dat <- merge(liftedData, GWAS, by = "SNP")
# Convert LogP to P-val
# dat$P <- 10^(-dat$LOG10P)
str(dat)
# Write
fwrite(dat, output_file, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE, compress = "gzip", nThread=numThreads)


