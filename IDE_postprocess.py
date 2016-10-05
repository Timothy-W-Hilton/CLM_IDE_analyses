import os
import sys
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


class Vividict(dict):
    """
    http://stackoverflow.com/questions/635483/what-is-the-best-way-to-implement-nested-dictionaries-in-python
    """
    def __missing__(self, key):
        value = self[key] = type(self)()
        return value


class CLMSiteContainer(object):
    """container for site-specific CLM data
    """
    def __init__(self, loc):
        """itialize run, var, loc

        """
        self.runs = []
        self.vars = []
        self.loc = loc
        self.data = Vividict()

    def add_var(self, run, var, df):
        self.data[run][var] = df

    def get_var(self, run, var, df):
        return data[run][var]


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
                                     data={'value': data,
                                           'date': self.dt,
                                           'time': self.time,
                                           'case': self.casename,
                                           'var': self.varname,
                                           'loc': loc.name})
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
    runs = ['CTL', 'IDE']
    vars = ['FPSN', 'WT', 'EFLX_LH_TOT_R']
    sp_info = IDE_locations.CLMf05g16_get_spatial_info()
    domain = sp_info[0]
    locs = sp_info[1:]
    data = Vividict()
    for this_loc in locs:
        sys.stdout.write('reading {}\n'.format(this_loc.name))
        for this_run in runs:
            sys.stdout.write('    {}: '.format(this_run))
            for this_var in vars:
                sys.stdout.write(' {} '.format(this_var))
                sys.stdout.flush()
                mp = MonthlyParser(this_run,
                                   data_dir,
                                   'IDE_{}_h1avg.nc'.format(this_run),
                                   this_var)
                mp.parse(loc=this_loc)
                data[this_run][this_loc.name][this_var] = mp
            sys.stdout.write('\n')
            sys.stdout.flush()

    for v in vars[0:1]:
        df = pd.concat([data[r][loc.name][v].data
                        for r in runs for loc in locs])
        with sns.axes_style("white"):
            g = sns.FacetGrid(df, hue='case', palette="Set1",
                              col='loc', col_wrap=3,
                              size=5, aspect=2,
                              hue_kws={"marker": ["^", "v"],
                                       "linestyle": ['-', '-']})
        g.map(plt.plot, 'date', 'value')
        g.fig.get_axes()[0].legend(loc='best')

    # alldata = pd.concat([data[r][v].data for r in runs for v in vars])

    # with sns.axes_style("white"):
    #     g = sns.FacetGrid(alldata, col='loc',
    #                       hue='case', palette="Set1",
    #                       size=2, aspect=1)
    # g.map(plt.plot, 'date', 'value')
