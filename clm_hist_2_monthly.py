import os
from nco import Nco
from datetime import datetime
import tempfile

import spinup_diagnostics as sd

if __name__ == "__main__":
    cases = ['IDE_ctl', 'IDE_redpcp']
    data_dirs = [os.path.join(os.getenv('CSCRATCH'), 'archive',
                              this_case, 'lnd', 'hist')
                 for this_case in cases]
    clm_runs = [sd.CLM_spinup_analyzer(data_dir=this_data,
                                       CASE=this_case)
                for this_data, this_case in zip(data_dirs, cases)]

    nco = Nco()

    tmpdir = tempfile.mkdtemp(dir=os.path.join(os.getenv('CSCRATCH')))
    print "placing monthly averages in {}".format(tmpdir)
    for this_year in range(1, 5):
        for this_month in range(1, 13):
            for this_run in clm_runs:
                pat = "*h2*{:04d}-{:02d}*.nc".format(this_year, this_month)
                this_run.gather_filenames(glob_pat=pat)
                if len(this_run.all_files) != 0:
                    out_fname = os.path.join(
                        tmpdir,
                        '{}_{:04d}-{:02d}_h2avg.nc'.format(this_run.CASE,
                                                           this_year,
                                                           this_month))
                    print 'calculating {}'.format(out_fname)
                    nco.ncra(input=this_run.all_files, output=out_fname)
