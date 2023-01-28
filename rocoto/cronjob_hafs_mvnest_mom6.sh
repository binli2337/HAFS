#!/bin/sh
set -x
date

# NOAA RDHPCS Hera
HOMEhafs=/scratch1/NCEPDEV/hwrf/save/${USER}/HAFS
dev="-s sites/hera.ent -f"
PYTHON3=/apps/intel/intelpython3/bin/python3
HOMEhafs=${HOMEhafs:-/lfs/h2/emc/hur/noscrub/${USER}/save/HAFS}
source ${HOMEhafs}/ush/hafs_pre_job.sh.inc

cd ${HOMEhafs}/rocoto
EXPT=$(basename ${HOMEhafs})
opts="-t -f"
scrubopt="config.scrub_work=no config.scrub_com=no"

#===============================================================================
# Example hafs moving nest experiments

# ./run_hafs.py ${opts} 2020082512 13L HISTORY \
#     config.EXPT=${EXPT} config.SUBEXPT=${EXPT}_C96_regional_1mvnest_storm \
#     config.NHRS=12 ${scrubopt} \
#     ../parm/hafs_C96_regional_1mvnest_storm.conf

 ./run_hafs.py ${opts} 2020082512 13L HISTORY \
     config.EXPT=${EXPT} config.SUBEXPT=${EXPT}_C96_regional_1mvnest_storm_mom6 \
     config.NHRS=12 ${scrubopt} \
     forecast.ccpp_suite_regional=FV3_HAFS_v1_gfdlmp_tedmf_nonsst \
     ../parm/hafs_C96_regional_1mvnest_storm_mom6.conf \
     ../parm/hafs_mom6.conf

# ./run_hafs.py ${opts} 2020082512 13L HISTORY \
#     config.EXPT=${EXPT} config.SUBEXPT=${EXPT}_C512_regional_1mvnest_storm \
#     config.NHRS=12 ${scrubopt} \
#     ../parm/hafs_C512_regional_1mvnest_storm.conf

# ./run_hafs.py ${opts} 2020082512 13L HISTORY \
#     config.EXPT=${EXPT} config.SUBEXPT=${EXPT}_C192_regional_1mvnest_atm_ocn \
#     config.NHRS=12 ${scrubopt} \
#     ../parm/hafs_C192_regional_1mvnest_storm.conf \

# ./run_hafs.py ${opts} 2020082512 13L HISTORY \
#     config.EXPT=${EXPT} config.SUBEXPT=${EXPT}_C512_regional_1mvnest_atm_ocn \
#     config.NHRS=12 ${scrubopt} \
#     ../parm/hafs_C512_regional_1mvnest_storm.conf \
#     ../parm/hafsv0p3_hycom.conf

# ./run_hafs.py ${opts} 2020082512 13L HISTORY \
#     config.EXPT=${EXPT} config.SUBEXPT=${EXPT}_C192_global_1mvnest_storm \
#     config.NHRS=12 ${scrubopt} \
#     ../parm/hafs_C192_global_1mvnest_storm.conf

# ./run_hafs.py ${opts} 2020082512 13L HISTORY \
#     config.EXPT=${EXPT} config.SUBEXPT=${EXPT}_C768_global_1mvnest_storm \ 
#     config.NHRS=12 ${scrubopt} \
#     ../parm/hafs_C768_global_1mvnest_storm.conf

#===============================================================================

date

echo 'cronjob done'
