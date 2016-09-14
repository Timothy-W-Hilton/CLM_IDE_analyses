import os
import spinup_diagnostics as sd
from nco import Nco
from datetime import datetime

if __name__ == "__main__":
    data_dir = os.path.join(os.getenv('CSCRATCH'), 'archive',
                            'CLM_f05_g16_control', 'lnd', 'hist')
    ctl_clm_run = sd.CLM_spinup_analyzer(data_dir=data_dir,
                                         CASE='CLM_f05_g16_control')
    ctl_clm_run.gather_filenames(glob_pat='*h1*.nc')

    nco = Nco()
    # t0 = datetime.now()
    # print "concatenating netcdf files"
    # nco.ncrcat(input=ctl_clm_run.all_files,
    #            options='-v FPSN',
    #            output="out.nc")
    # print "done concatenating (t = {})".format(datetime.now - t0)

    # ncra --mro -O -F -d time,6,,12,3 -n 150,4,1 1850.nc 1850_2009_JJA.nc
    t0 = datetime.now()
    print "concatenating netcdf files (start: {})".format(
        t0.strftime("%H:%M:%S"))
    varname = "H2OSOI"
    if len(ctl_clm_run.all_files) is not 0:
        nco.ncra(input=ctl_clm_run.all_files,
                 output="out_avg_{}.nc".format(varname),
                 options="--mro -O -F -d time,1,,24,24 -v {}".format(varname))
    print "done concatenating (t = {})".format(datetime.now() - t0)
