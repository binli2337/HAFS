#!/bin/sh
# prepare GFS forcing files for MOM6
#
# in the run directory, use this script to generate GFS forcing files starting from $CDATE 
# For example, if run_dir is HAFS_rt_regional_static_C192s1n4_atm_ocn/2020082512/13L/ocn_prep
# cd $run_dir
# ./gfs_global_forcing.sh $CDATE
# The generated forcing file can be used to make 126-hr forecast
# 
cp /scratch1/NCEPDEV/stmp2/Bin.Li/hafs_mom6/test_2020082512/DATM_INPUT/gfs_global_2020082512.nc .
exit
