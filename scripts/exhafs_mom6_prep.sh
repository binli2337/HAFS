#!/bin/sh

set -xe

NCP=${NCP:-'/bin/cp'}
NLN=${NLN:-'/bin/ln -sf'}
NDATE=${NDATE:-ndate}

TOTAL_TASKS=${TOTAL_TASKS:-2016}
NCTSK=${NCTSK:-12}
NCNODE=${NCNODE:-24}
OMP_NUM_THREADS=${OMP_NUM_THREADS:-2}
APRUNC=${APRUNC:-"aprun -b -j1 -n${TOTAL_TASKS} -N${NCTSK} -d${OMP_NUM_THREADS} -cc depth"}

#nest_grids=${nest_grids:-1}

CDATE=${CDATE:-${YMDH}}
cyc=${cyc:-00}
STORM=${STORM:-FAKE}
STORMID=${STORMID:-00L}

ymd=`echo $CDATE | cut -c 1-8`
month=`echo $CDATE | cut -c 5-6`
day=`echo $CDATE | cut -c 7-8`
hour=`echo $CDATE | cut -c 9-10`
CDATEprior=`${NDATE} -6 $CDATE`
PDY_prior=`echo ${CDATEprior} | cut -c1-8`
cyc_prior=`echo ${CDATEprior} | cut -c9-10`

WGRIB2=${WGRIB2:-wgrib2}
#CHGRESCUBEEXEC=${CHGRESCUBEEXEC:-${EXEChafs}/hafs_chgres_cube.x}

#ENSDA=${ENSDA:-NO}
#FGAT_MODEL=${FGAT_MODEL:-gfs}
#FGAT_HR=${FGAT_HR:-00}
mom6_intercom=${WORKhafs}/intercom/mom6

# Set options specific to the deterministic/ensemble forecast

# Generate the ICs and BC hour 0
if [ $gtype = regional ] ; then
  export REGIONAL=1
  export HALO=4
else
# for gtype = uniform, stretch, or nest
  export REGIONAL=0
fi
#  export INIDIR=${COMgfs}/gfs.$PDY/$cyc
export INIDIR=${COMgfs}/mom6.$PDY/$cyc/mom6

OUTDIR=${OUTDIR:-${WORKhafs}/intercom/mom6}
DATA=${DATA:-${WORKhafs}/mom6}
mkdir -p ${OUTDIR} ${DATA}
cd $DATA
#
#mom6 ic
#${USHhafs}/mom6_ic.sh $CDATE
#./mom6_ic.sh $CDATE
ln -sf /scratch1/NCEPDEV/hwrf/noscrub/Bin.Li/hafs-ic-bc/rtofs_ssh_ic*nc ${OUTDIR}
ln -sf /scratch1/NCEPDEV/hwrf/noscrub/Bin.Li/hafs-ic-bc/rtofs_ts_ic*nc ${OUTDIR}
#
#mom6 lbc
#${USHhafs}/mom6_lbc.sh $CDATE
#./mom6_lbc.sh $CDATE
ln -sf /scratch1/NCEPDEV/hwrf/noscrub/Bin.Li/hafs-ic-bc/obc*nc ${OUTDIR}
#mv obc_ts_*.nc ${OUTDIR}/obc_ts_*nc.nc
#mv obc_uv_*.nc ${OUTDIR}/obc_uv_*nc.nc
#mv obc_ssh_*.nc ${OUTDIR}/obc_ssh_*nc.nc
#
# gfs forcing for mom6
#${USHhafs}/gfs_global_forcing.sh $CDATE
#./gfs_global_forcing.sh $CDATE
# final output forcing files are in OUTDIR
cp /scratch1/NCEPDEV/stmp2/Bin.Li/hafs_mom6_tmp/test_2020082512_126h/DATM_INPUT/gfs_global_2020082512.nc ${OUTDIR}
#mv gfs_global_forcing$CDATE.nc ${OUTDIR}/gfs_global_forcing$CDATE.nc
exit
