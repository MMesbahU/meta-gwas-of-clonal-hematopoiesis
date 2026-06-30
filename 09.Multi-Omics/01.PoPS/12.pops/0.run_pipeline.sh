
cd /broad/hptmp/btruong/pops/chip_mesbah/data

trait=DNMT3A
trait=CHIP
trait=TET2


# $(NF-1)>0.01 && $(NF-1)<0.99 && 

zcat /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/multi_ancestry_summary/lifted_hg37.metaGWAS.${trait}.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.hg37_dbSNP.eaf001_min2Studies.tsv.gz | awk -v FS='\t' '{
  if (FNR==1) {print "SNP\tCHR\tPOS\tP\tN"; next}
  print $3,$1,$2,$11,$12
  
  }' OFS="\t" > ${trait}.snploc


############################ RUN MAGMA

trait=DNMT3A
trait=CHIP
trait=TET2
traitlist="DNMT3A CHIP TET2"

geneloc=/broad/hptmp/btruong/HDP/pops/NCBI37.3.gene.loc
geneloc=/medpop/esp2/projects/software/PoPS/data/gene_loc.txt


for trait in ${traitlist}; do
  qsub -v wdir=/broad/hptmp/btruong/pops/chip_mesbah/data -v trait=${trait} -v geneloc=${geneloc} $WDIR/scripts/pops/1.1.magma.sh
done

########################## run POPS

for trait in ${traitlist}; do
  qsub -v wdir=/broad/hptmp/btruong/pops/chip_mesbah/data -v trait=${trait} $WDIR/scripts/pops/1.2.pops.sh
done


