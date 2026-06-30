echo -e "variant_id\tCHROM\tPOS\tREF\tALT" > sort.overall_chip_variants.tsv; awk 'NR>1'  overall_chip_variants.tsv | sort -V | uniq | sed 's:\::\t:g'| awk '{print $1"_"$2"_"$3"_"$4"\t"$1"\t"$2"\t"$3"\t"$4}' >> sort.overall_chip_variants.tsv

