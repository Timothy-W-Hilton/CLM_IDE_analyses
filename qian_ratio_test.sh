#!/bin/sh


# look at year 1972 Qian et al (2006) preciptation output to test
# whether the IDE/CTL ratio is constant in time.  It appears to be
# constant - the difference between the minimum and maximum ratio
# values in each grid cell (variable diff in the output file
# qian_minmax.nc) maxes out at 1e-6, seven orders of magnitudes less
# than the ratio values.

export WORKDIR=$CSCRATCH/ratio_test
cd $WORKDIR

export DQIANIDE="$CSCRATCH/Qian_pcp_reduced"
export DQIAN=/project/projectdirs/ccsm1/inputdata/atm/datm7/atm_forcing.datm7.Qian.T62.c080727/Precip6Hrly

for f in $DQIAN/*1972*.nc
do
    echo "$(basename $f)"
    ncks -O --mk_rec_dmn time $f $(basename $f)
done
echo "concatenating Qian PCP"
ncrcat *1972*.nc qian72.nc
rm -fv *1972*.nc
echo "done concatenating"

for f in $DQIANIDE/*1972*.nc
do
    echo "$(basename $f)"
    ncks -O --mk_rec_dmn time $f $(basename $f)
done
echo "concatenating Qian IDE PCP"
ncrcat ./*1972*.nc qianide72.nc
rm -f *1972*.nc


ncrename -v PRECTmms,PRECTmmsIDE qianide72.nc
ncks -A -v PRECTmmsIDE qianide72.nc qian72.nc
rm -fv qianide72.nc
echo "calculating ratio"
ncap2 -s "ratio=PRECTmmsIDE/PRECTmms" qian72.nc ratio.nc
rm -fv qian72.nc

# ncra -y rms -v prs_sfc in.nc foo.nc

ncra -y min -v ratio ratio.nc qian_minmax.nc
ncrename -v ratio,ratio_min qian_minmax.nc
ncra -y max -v ratio ratio.nc max.nc
ncrename -v ratio,ratio_max max.nc
ncks -A -v ratio_max max.nc qian_minmax.nc
rm -fv max.nc
ncap2 -A -s "diff=ratio_max-ratio_min" qian_minmax.nc
# ncap2 -s 'min=min(ratio)' ratio.nc qian_minmax.nc
# ncap2 -A -s 'max=max(ratio)' ratio.nc qian_minmax.nc
