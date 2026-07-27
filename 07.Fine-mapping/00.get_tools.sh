#!/bin/bash

# FINEMAP
# wget http://www.christianbenner.com/finemap_v1.4.2_x86_64.tgz
# tar -xzf finemap_v1.4_x86_64.tgz
FINEMAP="/medpop/esp2/mesbah/tools/finemap_v1.4.2_x86_64/finemap_v1.4.2_x86_64"
${FINEMAP} --help

# bgenix
use .bgen-1.1.4

bgenix -help

# qctool2
QCTOOL2="/medpop/esp2/mesbah/tools/qcTool/bin/qctool_v2.0.7"

# LDstore2
## wget https://bitbucket.org/gavinband/ldstore/downloads/ldstore_v2.0_x86_64.tgz
# wget http://www.christianbenner.com/ldstore_v2.0_x86_64.tgz
# tar -xzf ldstore_v2.0_x86_64.tgz
LDstore2="/medpop/esp2/mesbah/tools/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64" # --help

# plink2
PLINK2="/medpop/esp2/mesbah/tools/plink2_linux_x86_64_20200831"



