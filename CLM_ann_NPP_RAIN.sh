#!/bin/sh

# calculate annual rain and NPP from CLM output

export WORKDIR=$CSCRATCH/ann_RAIN_NPP
mkdir -pv $WORKDIR
cd $WORKDIR

export CTL="$CSCRATCH/archive/IDE_ctl/lnd/hist"
export IDE="$CSCRATCH/archive/IDE_redpcp/lnd/hist"
export OUTFILE=$WORKDIR/NPP_RAIN.nc

echo "concatenating CTL NPP, RAIN"
rm -fv $OUTFILE
ncrcat -h -v RAIN,FPSN $CTL/*h0*.nc $OUTFILE
ncrename -h -v RAIN,RAINctl $OUTFILE
ncrename -h -v FPSN,NPPctl $OUTFILE
echo "concatenating IDE NPP, RAIN"
ncrcat -A -h -v RAIN,FPSN $IDE/*h0*.nc $OUTFILE
ncrename -h -v RAIN,RAINide $OUTFILE
ncrename -h -v FPSN,NPPide $OUTFILE
echo "done concatenating"

echo "calculating nsecs"
ncap2 -A -s "nsecs=(time_bounds(:,1)-time_bounds(:,0))*60*60*24" -v $OUTFILE $OUTFILE
echo "done calculating nsecs"

echo "integrating NPP, RAIN"
ncap2 -A -s "umol2g=12.001*1e-6;NPPctl=NPPctl*nsecs*umol2g;NPPide=NPPide*nsecs*umol2g;RAINctl=RAINctl*nsecs;RAINide=RAINide*nsecs" $OUTFILE $OUTFILE
echo "done integrating"

echo "calculating annual sums "
ncra -O -v NPPide,NPPctl,RAINide,RAINctl --op_typ=ttl --mro -d time,,,12,12 $OUTFILE $OUTFILE
echo "done"
ncap2 -A -s "time={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}" -v $OUTFILE $OUTFILE
ncatted -h -a units,NPPide,o,c,'gC m-2 yr-1' $OUTFILE
ncatted -h -a units,NPPctl,o,c,'gC m-2 yr-1' $OUTFILE
ncatted -h -a units,RAINide,o,c,'mm yr-1' $OUTFILE
ncatted -h -a units,RAINctl,o,c,'mm yr-1' $OUTFILE
ncatted -h -a units,time,o,c,'year of simulation' $OUTFILE
