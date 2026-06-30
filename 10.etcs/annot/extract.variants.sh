use Bcftools

while read trait; do echo -e "CHR\tPOS\tMarkerID\tREF\tATL\tEffectAllele\tOtherAllele\tEAF\tBETA\tSE\tP\tN\tDirection\tHet_P\tFunc.refGene\tGene.refGene\tAAChange.refGene\tExonicFunc.refGene" > tsv.p5e6eaf001nstd2.metagwas.has${trait}.AoU250k_UKBtopImp200k250k_topmed74k_mgbb53k_biovu54k.Oct2023.tsv; bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/EffectAllele\t%INFO/OtherAllele\t%INFO/EAF\t%INFO/BETA\t%INFO/SE\t%INFO/Pvalue\t%INFO/N_METAL\t%INFO/Effect_Direction\t%INFO/Q_Pvalue\t%INFO/Func.refGene\t%INFO/Gene.refGene\t%INFO/AAChange.refGene\t%INFO/ExonicFunc.refGene\n' annovar.p5e6eaf001nstd2.metagwas.has${trait}.AoU250k_UKBtopImp200k250k_topmed74k_mgbb53k_biovu54k.Oct2023.hg38_multianno.vcf.gz >> tsv.p5e6eaf001nstd2.metagwas.has${trait}.AoU250k_UKBtopImp200k250k_topmed74k_mgbb53k_biovu54k.Oct2023.tsv; done < <(echo -e "CH\nDNMT3A\nTET2")


