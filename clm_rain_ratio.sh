#!/bin/sh

# look at year 0001 February CLM RAIN output to test whether the
# IDE/CTL ratio is constant in time.  It is for the PCP data that
# drive CLM.  It appears to vary some in time after it goes through
# CLM.  This seems curious.

export WORKDIR=$CSCRATCH/ratio_test
cd $WORKDIR

export CTL="$CSCRATCH/archive/IDE_ctl/lnd/hist"
export IDE="$CSCRATCH/archive/IDE_redpcp/lnd/hist"

echo "concatenating CTL PCP"
ncrcat -v RAIN $CTL/*h2.0001-02*.nc ctl.nc
echo "done concatenating"

echo "concatenating IDE PCP"
ncrcat -v RAIN $IDE/*h2.0001-02*.nc ide.nc
echo "done concatenating"

ncrename -v RAIN,RAINctl ctl.nc
ncrename -v RAIN,RAINide ide.nc
mv -v ctl.nc clm_rain.nc
ncks -A -v RAINide ide.nc clm_rain.nc
rm -fv ide.nc
echo "calculating ratio"
ncap2 -A -s "ratio=RAINide/RAINctl" clm_rain.nc

ncra -y min -v ratio clm_rain.nc clm_minmax.nc
ncrename -v ratio,ratio_min clm_minmax.nc
ncra -y max -v ratio clm_rain.nc max.nc
ncrename -v ratio,ratio_max max.nc
ncks -A -v ratio_max max.nc clm_minmax.nc
rm -fv max.nc
ncap2 -A -s "diff=ratio_max-ratio_min" clm_minmax.nc
