#!/bin/bash

# $1 path with file parameters
# $2 path source package

# source LSST
source /opt/lsst/software/stack/loadLSST.bash
setup lsst_distrib

# start main script loop on raw image lsst with some number visit
export OPENBLAS_NUM_THREADS=1
export PATH_SCRIPT=$2
cd $PATH_SCRIPT
ls -l 
source init_stage1.bash
#echo $PATH
#echo $PYTHONPATH
echo "SGE_TASK_ID:"
echo $SGE_TASK_ID

# lancement du job lsst
echo "path file parameters:"
echo $1
script_lsst_step.py $1

