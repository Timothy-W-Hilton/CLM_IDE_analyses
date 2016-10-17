import os
import sys
import calendar
import netCDF4
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
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
        self.vunits = ''

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
        self.lname = nc.variables[self.varname].long_name
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


def format_factorgrid(g, var_sname, var_lname, var_units, x_var_name):
    g.set_axis_labels(x_var='year of simulation',
                      y_var='{var} ({units})'.format(
                          var=var_sname,
                          units=var_units))
    g.set_titles(template='{col_name}')
    g.fig.get_axes()[0].legend(loc='best')
    plt.figtext(0.5, 0.99,
                '{shortname}: {longname}'.format(
                    shortname=var_sname,
                    longname=var_lname),
                horizontalalignment='center')
    return g


def carbon_umol_m2_s_2_g_m2_yr(Cin, ndays):
    """convert umol m-2 s-1 to g m-2 yr-1"""
    C_mol_wt = 12.001
    mol_per_umol = 1e-6
    secs_per_day = 60 * 60 * 24
    nsecs = ndays * secs_per_day
    return (Cin * C_mol_wt * mol_per_umol * nsecs)


def calc_dvar(df, varname):
    """calculate mean monthly IDE-CTL difference in a variable
    """
    mean_var = df[df['var'] == varname][
        ['value', 'case', 'loc', 'month']].groupby(
            ['case', 'loc', 'month']).mean()
    dvar = mean_var.ix['CTL'] - mean_var.ix['IDE']
    dvar['ndays'] = map(lambda x: calendar.monthrange(2001, x)[1],
                        dvar.index.get_level_values('month'))
    dvar['gC'] = carbon_umol_m2_s_2_g_m2_yr(dvar['value'], dvar['ndays'])
    annual_total = dvar.groupby(level='loc').sum()
    return dvar, annual_total


if __name__ == "__main__":
    calculate_data = False
    data_dir = os.path.join(os.getenv('CSCRATCH'), 'monthly_means')
    runs = ['CTL', 'IDE']
    # TODO: H2OSOI & other variables with soil depth dimension
    h1vars = ['FPSN', 'WT', 'EFLX_LH_TOT_R', 'FCTR', 'FGEV', 'FIRA', 'FSH',
            'FSH_V', 'H2OSNO', 'QBOT', 'QCHARGE', 'QDRAI',
            'QINFL', 'QVEGT', 'TBOT', 'ZWT']
    h2vars = ['BTRAN', 'FSDS', 'QOVER', 'QRUNOFF', 'QVEGE', 'RAIN']
    sp_info = IDE_locations.CLMf05g16_get_spatial_info()
    domain = sp_info[0]
    locs = sp_info[1:]

    cal_wet_to_dry = ['McLaughlin NRS',
                      'Sierra Foothill Research Extension Center',
                      'Younger Lagoon',
                      'Sedgewick NRS',
                      'Loma Ridge Global Change Experiment',
                      'Box Springs']
    cal_locs = cal_wet_to_dry
    for i, s in enumerate(cal_wet_to_dry):
        cal_locs[i] = next((this for this in locs if this.name == s), None)

    if calculate_data:
        data = Vividict()
        units = Vividict()
        df_list = Vividict()
        for which_hist in [1, 2]:
            vars = locals()['h{}vars'.format(which_hist)]
            for this_loc in cal_locs:
                sys.stdout.write('reading {}\n'.format(this_loc.name))
                for this_run in runs:
                    sys.stdout.write('    {}: '.format(this_run))
                    for this_var in vars:
                        sys.stdout.write(' {} '.format(this_var))
                        sys.stdout.flush()
                        mp = MonthlyParser(this_run,
                                           data_dir,
                                           'IDE_{}_h{}avg.nc'.format(
                                               this_run, which_hist),
                                           this_var)
                        mp.parse(loc=this_loc)
                        data[this_run][this_loc.name][this_var] = mp
                        units[this_run][this_loc.name][this_var] = mp.vunits
                    sys.stdout.write('\n')
                    sys.stdout.flush()
        # plt.rcParams['figure.figsize']=(10,10)
        for v in (h1vars + h2vars):
            sys.stdout.write('{} '.format(v))
            sys.stdout.flush()
            df = pd.concat([data[r][loc.name][v].data
                            for r in runs for loc in cal_locs])
            df['fyear'] = df['time'] / 365.0
            df['month'] = df.index.month
            df['units'] = data[r][loc.name][v].vunits
            df['var_lname'] = data[r][loc.name][v].lname
            df_list[v] = df
        sys.stdout.write('\n')
        sys.stdout.flush()
        all_vars = pd.concat(df_list).reset_index(drop=True)
        all_vars.to_csv('monthly_vals.txt')
    else:
        all_vars = pd.read_csv('./monthly_vals.txt')

    sys.stdout.write('plotting ')
    for v in (h1vars + h2vars):
        sys.stdout.write('{} '.format(v))
        sys.stdout.flush()
        df = all_vars[all_vars['var'] == v]
        with sns.axes_style("white"):
            g = sns.factorplot(data=df, x='month', y='value',
                               col='loc', hue='case',
                               kind='point', col_wrap=3, legend=False,
                               palette="Dark2")
            format_factorgrid(g,
                              df['var'].iloc[0],
                              df['var_lname'].iloc[0],
                              df['units'].iloc[0],
                              'month')
            g.savefig(os.path.join(os.getenv('CSCRATCH'), 'plots',
                                   '{var}_bymonth.pdf'.format(var=v)))
            plt.close(g.fig)
        with sns.axes_style("white"):
            g = sns.FacetGrid(df, hue='case', palette="Dark2",
                              col='loc', col_wrap=3,
                              size=3, aspect=1.5,
                              hue_kws={"marker": ["^", "v"],
                                       "linestyle": ['-', '-']})
            g.map(plt.plot, 'fyear', 'value')
            format_factorgrid(g,
                              df['var'].iloc[0],
                              df['var_lname'].iloc[0],
                              df['units'].iloc[0],
                              'year of simulation')
            g.savefig(os.path.join(os.getenv('CSCRATCH'), 'plots',
                                   '{var}_timeseries.pdf'.format(var=v)))
            plt.close(g.fig)
    sys.stdout.write('\n')
    sys.stdout.flush()
