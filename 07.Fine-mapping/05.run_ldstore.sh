#!/bin/bash

# error: no match in rsid
# bgenix -g /broad/hptmp/mesbah/RefSeq/mgbb53k/region/chr22_27865160_29554205.region.bgen -list  > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/mgbb/chr22_27865160_29554205.region.bgen_ids.txt
# solution: have to match chromosome e.g. chr22 in bgen is only 22. 
# added: awk 'BEGIN{OFS=" "} NR==1{print; next} {sub(/^chr/,"", $2); print }'

# bash /medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/07.Fine-mapping/05.run_ldstore.sh "/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/mgbb/chr22_27865160_29554205_chr22:28707610:T:C.hasCH.chr1_22.MultiANC.z" /broad/hptmp/mesbah/RefSeq/mgbb53k/ldstore /broad/hptmp/mesbah/RefSeq/mgbb53k/region/chr22_27865160_29554205.region.bgen /broad/hptmp/mesbah/RefSeq/mgbb53k/region/chr22_27865160_29554205.region.bgen.bgi 50305

set -euo pipefail

    #############################################
    # STEP 2: LDSTORE2 (recommended)
    #############################################
# ldstore z file
# rsid chromosome position allele1 allele2
# mkdir -p /broad/hptmp/mesbah/RefSeq/mgbb53k/ldstore
gwas_Z=${1}
ldstore_out=${2}
ldstore_Z=${ldstore_out}/ldstoreZ.$(basename ${gwas_Z})
bgen=${3}
bgi=${4}
n_samples=${5}
ldstore_master=${ldstore_out}/ldstore.$(basename ${gwas_Z} ".z").master
ldstore_bdose=${ldstore_out}/ldstore.$(basename ${gwas_Z} ".z").bdose
ldstore_bcor=${ldstore_out}/ldstore.$(basename ${gwas_Z} ".z").bcor
LDstore2="/medpop/esp2/mesbah/tools/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64" 
##
cut -d ' ' -f1-5 ${gwas_Z} | awk 'BEGIN{OFS=" "} NR==1{print; next} {sub(/^chr/,"", $2); print }' > ${ldstore_Z}
## Master file LDstore
cat > ${ldstore_master} << EOF
z;bgen;bgi;bcor;bdose;n_samples
${ldstore_Z};${bgen};${bgi};${ldstore_bcor};${ldstore_bdose};${n_samples}
EOF
## --n-threads ${thr} --memory ${mem}
${LDstore2} \
	--in-files ${ldstore_master} \
        --write-bdose \
        --write-bcor
##
rm ${ldstore_bdose}

#############################################
