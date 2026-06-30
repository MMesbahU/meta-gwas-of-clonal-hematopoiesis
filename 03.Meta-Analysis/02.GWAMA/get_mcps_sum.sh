# summary
echo -e "CHIP\thttps://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90435001-GCST90436000/GCST90435341/GCST90435341.tsv\nDNMT3A\thttps://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90435001-GCST90436000/GCST90435342/GCST90435342.tsv\nTET2\thttps://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90435001-GCST90436000/GCST90435343/GCST90435343.tsv\nASXL1\thttps://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90435001-GCST90436000/GCST90435344/GCST90435344.tsv\nPPM1D\thttps://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90435001-GCST90436000/GCST90435345/GCST90435345.tsv\nTP53\thttps://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90435001-GCST90436000/GCST90435346/GCST90435346.tsv\nSF3B1\thttps://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90435001-GCST90436000/GCST90435347/GCST90435347.tsv\nSRSF2\thttps://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90435001-GCST90436000/GCST90435348/GCST90435348.tsv\nSF3B1_SRSF2\thttps://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90435001-GCST90436000/GCST90435352/GCST90435352.tsv" > mcps.file_list.txt

#
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90435001-GCST90436000/GCST90435341/GCST90435341.tsv-meta.yaml


## download files
while read lines; do wget -O /broad/hptmp/mesbah/dataset/mcps/$(echo $lines | awk '{print $1}').$(basename $(echo $lines | awk '{print $2}') ) $(echo $lines | awk '{print $2}') && gzip /broad/hptmp/mesbah/dataset/mcps/$(echo $lines | awk '{print $1}').$(basename $(echo $lines | awk '{print $2}') ); done <mcps.file_list.txt 1>>my_log.txt 2>>my_err.txt &

## prepare GWAMA input
while read files; do zcat CHIP.GCST90435341.tsv.gz | awk 'NR==1{print "SNPID\tREF\tALT\tAAF\tBETA\tSE\tP\tN"}(NR>1){print "chr"$1":"$2":"$4":"$3"\t"$4"\t"$3"\t"$7"\t"$5"\t"$6"\t"10^-$8"\t136401"}' > /broad/hptmp/mesbah/dataset/mcps/gwama.$(basename ${files} ".gz"); done < <(ls -lhrt /broad/hptmp/mesbah/dataset/mcps/{CHIP,DNMT3A,TET2,ASXL1,TP53,PPM1D,SF3B1_SRSF2}.*.tsv.gz | awk '{print $NF}') &

