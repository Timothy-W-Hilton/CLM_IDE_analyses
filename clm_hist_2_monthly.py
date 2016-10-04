import os
from nco import Nco
from datetime import datetime
import tempfile

import spinup_diagnostics as sd


def month_year_iter(start_month, start_year, end_month, end_year):
    """posted to:
    http://stackoverflow.com/questions/5734438/how-to-create-a-month-iterator
    """
    ym_start = 12 * start_year + start_month - 1
    ym_end = 12 * end_year + end_month - 1
    for ym in range(ym_start, ym_end):
        y, m = divmod(ym, 12)
        yield y, m+1

if __name__ == "__main__":
    which_hist = 'h1'
    cases = ['IDE_ctl', 'IDE_redpcp']
    cases = ['IDE_redpcp']
    data_dirs = [os.path.join(os.getenv('CSCRATCH'), 'archive',
                              this_case, 'lnd', 'hist')
                 for this_case in cases]
    clm_runs = [sd.CLM_spinup_analyzer(data_dir=this_data,
                                       CASE=this_case)
                for this_data, this_case in zip(data_dirs, cases)]

    nco = Nco()

    tmpdir = tempfile.mkdtemp(prefix='monthly_mean_tmp',
                              dir=os.path.join(os.getenv('CSCRATCH')))
    print "placing monthly averages in {}".format(tmpdir)
    m0, y0 = (12, 1)  # start month, start year
    m1, y1 = (6, 4)  # end month, end year

    for this_run in clm_runs:
        out_files = []
        for this_year, this_month in month_year_iter(m0, y0, m1, y1):
            pat = "*{}*{:04d}-{:02d}*.nc".format(which_hist,
                                                 this_year,
                                                 this_month)
            this_run.gather_filenames(glob_pat=pat)
            if len(this_run.all_files) != 0:
                out_fname = os.path.join(
                    tmpdir,
                    '{}_{:04d}-{:02d}_{}avg.nc'.format(this_run.CASE,
                                                       this_year,
                                                       this_month,
                                                       which_hist))
                print 'calculating {}'.format(out_fname)
                nco.ncra(input=this_run.all_files, output=out_fname)
                out_files.append(out_fname)
        print "concatenating {} monthly averages".format(this_run.CASE)
        nco.ncrcat(input=out_files,
                   output=os.path.join(
                       tmpdir,
                       '{}_{}avg.nc'.format(this_run.CASE, which_hist)))
