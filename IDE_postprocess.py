import os
import netCDF4
import numpy as np
import itertools
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from spinup_diagnostics import CLM_spinup_analyzer as CSA
import IDE_locations
from datetime import datetime, timedelta


def month_from_t0_offset(ndays, t0):
    return (timedelta(days=float(ndays)) + t0)


class MonthlyParser(object):
    """parses monthly mean CLM data
    """
    def __init__(self, casename, datadir, fname, varname):
        self.datadir = datadir
        self.casename = casename
        self.fname = fname
        self.data_file = os.path.join(self.datadir, self.fname)
        self.varname = varname

    def parse(self, t0=datetime(2001, 1, 1, 0, 0, 0), loc=None):
        """parse a variable and time stamp from netCDF
        """
        nc = netCDF4.Dataset(os.path.join(self.datadir, self.fname))
        self.time = nc.variables['time'][:]
        self.tunits = nc.variables['time'].units
        self.t0 = t0
        self.calc_moy()
        if loc is None:
            self.data = nc.variables[self.varname][...]
        else:
            data = nc.variables[self.varname][:, loc.clm_y, loc.clm_x]
            self.data = pd.DataFrame(index=self.dt,
                                     data={self.varname: data,
                                           'time': self.time,
                                           'case': self.casename})
        self.vunits = nc.variables[self.varname].units
        nc.close()

    def calc_moy(self):
        """calculate month of year from timestamp
        """
        # doing date conversions this way because numpy.timedelta64
        # cannot handle fractional days (!!!), and with pd.TimeStamp:
        # "Since pandas represents timestamps in nanosecond
        # resolution, the timespan that can be represented using a
        # 64-bit integer is limited to approximately 584 years:"
        # http://pandas-docs.github.io/pandas-docs-travis/timeseries.html#timestamp-limitations
        time_list = list(self.time)
        self.dt = [month_from_t0_offset(t, self.t0) for t in time_list]
        self.moy = np.array([this.month for this in self.dt])

    def get_data(self):
        return self.data

    def get_time(self):
        return self.time

    def get_moy(self):
        return self.moy


# if __name__ == "__main__":
#     d_archive = os.path.join(os.getenv('CSCRATCH'), 'monthly_means')
#     locs = IDE_locations.CLMf05g16_get_spatial_info()
#     d = locs[0]
#     locs = locs[1:]
#     runs = {'ctl': CSA(os.path.join(d_archive, '),
#                        'ide_ctl'),
#             'ide': CSA(os.path.join(d_archive, 'IDE_redpcp',
#                                     'lnd', 'hist'),
#                        'ide_redpcp')}
#     mavg = []

#     s_per_tstep = 60 * 60 * 6  # tsteps are 6 hours

#     data = {}
#     for this_loc in locs:
#         vars = {'RAIN': {}, 'QVEGE': {}, 'QVEGT': {}}
#         for run_name, run in runs.iteritems():
#             for this_var in vars.keys():
#                 run.gather_filenames(glob_pat="*h2*.nc")
#                 vars[this_var][run_name] = run.parse_var(
#                     this_var,
#                     lon_idx=this_loc.clm_x,
#                     lat_idx=this_loc.clm_y)
#                 vars[this_var][run_name].data *= s_per_tstep
#                 vars[this_var][run_name].units = (
#                     vars[this_var][run_name].units.replace('/s', ''))
#         data[this_loc.name] = vars

#     # for this_site, this_data in data.items():
#     #     fig, ax = plt.subplots(ncols=1, nrows=3, figsize=(7, 10))
#     #     for this_var, this_ax in zip(this_data.values(), ax):
#     #         this_ax.plot(this_var['ctl'].data.squeeze(), 'b-', label='CTL')
#     #         this_ax.plot(this_var['ide'].data.squeeze(), 'k--', label='IDE')
#     #         this_ax.set_ylabel('{} ({})'.format(this_var['ctl'].varname,
#     #                                             this_var['ctl'].units))
#     #     ax[-1].set_xlabel('day of simulation')
#     #     ax[0].set_title('{}'.format(this_site))
#     #     ax[0].legend(loc='best')
#     #     fig.tight_layout()
#     #     fig.savefig(os.path.join(os.getenv('CSCRATCH'),
#     #                              'plots_4year',
#     #                              '{}.pdf'.format(this_site.replace(' ', ''))))
#     #     # plt.show()


if __name__ == "__main__":
    data_dir = os.path.join(os.getenv('CSCRATCH'), 'monthly_means')
    runs = ['ctl', 'redpcp']
    vars = ['FPSN', 'WT']
    locs = IDE_locations.CLMf05g16_get_spatial_info()
    data = dict.fromkeys(runs, None)
    for this_run in data.keys():
        data[this_run] = {}
        for this_var in vars:
            mp = MonthlyParser(this_run, data_dir,
                               'IDE_{}_h1avg.nc'.format(this_run), this_var)
            mp.parse(loc=locs[1])
            data[this_run][this_var] = mp


    # x, y = (locs[1].clm_x, locs[1].clm_y)
    # plt.figure()
    # plt.plot(data['ctl']['FPSN'].time, data['ctl']['FPSN'].data[:, y, x],
    #          data['redpcp']['FPSN'].time, data['redpcp']['FPSN'].data[:, y, x])
    # plt.figure()
    # plt.plot(data['ctl']['WT'].time, data['ctl']['WT'].data[:, y, x],
    #          data['redpcp']['WT'].time, data['redpcp']['WT'].data[:, y, x])

    for v in vars:
        df = pd.concat([data[r][v].data for r in runs])
        with sns.axes_style("white"):
            g = sns.FacetGrid(df, hue='case', palette="Set1",
                              size=5, hue_kws={"marker": ["^", "v"],
                                               "linestyle": ['-', '-']})
        g.map(plt.plot, "time", v)
        g.add_legend()
