#! /opt/hpc/bin/bash

brainno=$1

source /sonas-hs/it/hpc/home/easybuild/lmod-setup.sh

module load foss/2016a
module load IntelPython/2.7.12

python brain_region.py $brainno
