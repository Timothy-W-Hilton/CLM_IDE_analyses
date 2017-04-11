#!/bin/sh

#SBATCH -p regular
#SBATCH -N 1
#SBATCH -t 03:00:00
#SBATCH -L cscratch1

export VARNAME="$1"
export RUNNAME="$2"
export HISTSUFFIX="$3"
# export INDIR="/global/cscratch1/sd/twhilton/daily_CLM_output/input"
export INDIR="/global/cscratch1/sd/twhilton/archive/IDE_$RUNNAME/lnd/hist/"
export OUTDIR="/global/cscratch1/sd/twhilton/daily_CLM_output/output"
export FNAME_TMP="$OUTDIR/IDE_${RUNNAME}_6hrly_${VARNAME}.nc"
export FNAME_FINAL="$OUTDIR/IDE_${RUNNAME}_daily_${VARNAME}.nc"

module load nco
echo 'starting' `date`
echo "$FNAME_TMP $FNAME_FINAL"
rm -fv $FNAME_TMP $FNAME_FINAL
# production: next line processes all data
ncrcat -v $VARNAME $INDIR/IDE_${RUNNAME}.clm2.${HISTSUFFIX}.* $FNAME_TMP
# debug: next line only uses first ten days of data
# ncrcat -v $VARNAME $INDIR/IDE_${RUNNAME}.clm2.${HISTSUFFIX}.0001-01-0* $FNAME_TMP
ncra -O --mro -d time,1,,4,4 $FNAME_TMP $FNAME_FINAL
rm -fv $FNAME_TMP
echo 'done' `date`
