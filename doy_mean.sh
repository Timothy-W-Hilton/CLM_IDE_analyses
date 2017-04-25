#!/bin/sh

VAR=$1

for run in "redpcp" "ctl"; do
    echo "working on $run $var"
    for nday in {1..365}; do
	doy=`printf "%03d" $nday`
	echo "doy: ${doy}"
	ncra -h -O -F -d time,${doy},,365 IDE_${run}_daily_${VAR}.nc IDE_${run}_daily_${VAR}_doy${doy}.nc
    done
    echo "concatenating"
    ncrcat -O -h IDE_${run}_daily_${VAR}_doy*.nc IDE_${run}_daily_${VAR}_doymean.nc
    echo "removing temporary files"
    rm -fv IDE_${run}_daily_${VAR}_doy[0-9]*.nc
done
