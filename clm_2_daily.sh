#!/bin/sh

#SBATCH -p debug
#SBATCH -N 1
#SBATCH -t 00:03:00
#SBATCH -L cscratch1

export VARNAME='FPSN'
export RUNNAME='redpcp'
export INDIR="/global/cscratch1/sd/twhilton/daily_CLM_output/input"
# export INDIR="/global/cscratch1/sd/twhilton/archive/IDE_$RUNNAME/lnd/hist/"
export OUTDIR="/global/cscratch1/sd/twhilton/daily_CLM_output/tmp"
export FNAME_TMP="$OUTDIR/IDE_${RUNNAME}_6hrly_${VARNAME}.nc"
export FNAME_FINAL="$OUTDIR/IDE_${RUNNAME}_daily_${VARNAME}.nc"

module load nco
echo 'starting' `date`
echo "$FNAME_TMP $FNAME_FINAL"
rm -fv $FNAME_TMP $FNAME_FINAL
ncrcat -v $VARNAME $INDIR/IDE_${RUNNAME}.clm2.h1.0001-01-0* $FNAME_TMP
ncra -O --mro -d time,1,,4,4 $FNAME_TMP $FNAME_FINAL
rm -fv $FNAME_TMP
echo 'done' `date`
