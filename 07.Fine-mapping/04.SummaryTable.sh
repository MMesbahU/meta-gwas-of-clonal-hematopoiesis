#!/bin/bash

find finemap_output -name "*.snp" | while read f; do
    locus=$(basename $(dirname $f))

    awk -v l=$locus '
    NR>1{
        print l,$1,$2
    }' $f

done > finemap_pip_summary.txt