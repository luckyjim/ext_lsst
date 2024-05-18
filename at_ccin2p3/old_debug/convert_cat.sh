#!/bin/sh
. /cvmfs/euclid.in2p3.fr/CentOS7/EDEN-2.1/etc/profile.d/euclid.sh
 
export PYTHONPATH=$PYTHONPATH:/sps/euclid/Users/colley/wd_euc/SIM_EXT/SIM_EXT/python

convert_cat.py $1 $2
