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

import LomaRidgeTools


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

    def parse(self, t0=datetime(2001, 1, 1, 0, 0, 0), loc=None, lev=None):
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
            if lev is None:
                data = nc.variables[self.varname][:, loc.clm_y, loc.clm_x]
            else:
                data = nc.variables[self.varname][:, lev,
                                                  loc.clm_y, loc.clm_x]
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


class DailyVarCsvMaker(object):
    """makes a daily-resolution CSV file for a CLM variable
    """
    def __init__(self, varname, runname):
        """class constructor

        ARGS:
        varname (str): CLM variable name
        runname (str): CLM run name
        """
        self.varname = varname
        self.runname = runname

    def make_csv(self, runs=['IDE_ctl', 'IDE_redpcp'], lev=None):
        """produce a CSV file for specified variable

        ARGS:
        runs (iterable): list or tuple of strings containing CLM run
            names to process
        lev (int): optional vertical level to extract to CSV.  For
            variables without vertical levels use None. Default is
            None.
        """
        sp_info = IDE_locations.CLMf05g16_get_spatial_info()
        locs = sp_info[1:]

        df_all = None
        for this_site in locs:
            for this_run in runs:
                print "parsing {} data".format(this_site.name)
                mp = MonthlyParser(this_run,
                                   os.path.join('/', 'global',
                                                'cscratch1',
                                                'sd',
                                                'twhilton/',
                                                'daily_CLM_output',
                                                'output'),
                                   '{runname}_daily_{varname}.nc'.format(
                                       runname=this_run,
                                       varname=self.varname),
                                   self.varname)
                mp.parse(loc=this_site, lev=lev)
                mp.data['doy'] = mp.data.index.dayofyear
                if df_all is None:
                    df_all = mp.data
                else:
                    df_all = pd.concat((df_all, mp.data))
        if lev is not None:
            suffix = "_{:02d}".format(lev)
        else:
            suffix = ""
        df_all.to_csv('./{varname}_daily_all.csv.gz'.format(
            varname=self.varname + suffix),
                      compression='gzip')


def plot_site_annual_rain_gpp(all_vars, locs):
    """plot annual rain vs annual GPP for a site
    """
    SECS_PER_DAY = 24 * 60 * 60
    site_data = all_vars[all_vars['var'].isin(['RAIN', 'FPSN'])].copy()
    site_data['year'] = np.floor(site_data['fyear']).astype('int')
    # site_data = site_data.ix[site_data['year'] < 7]
    site_data.loc[:, 'ndays'] = map(lambda x: calendar.monthrange(2001, x)[1],
                                    site_data['month'])
    rain_idx = site_data['var'] == 'RAIN'
    site_data.loc[rain_idx, 'annual_'] = (site_data.loc[rain_idx, 'value'] *
                                          site_data.loc[rain_idx, 'ndays'] *
                                          SECS_PER_DAY)
    fpsn_idx = site_data['var'] == 'FPSN'
    site_data.loc[fpsn_idx, 'annual_'] = carbon_umol_m2_s_2_g_m2_yr(
        site_data.loc[fpsn_idx, 'value'],
        site_data.loc[fpsn_idx, 'ndays'])
    grp = site_data.groupby(['year', 'loc', 'var', 'case']).sum()
    totals = grp['annual_'].reset_index()

    g = dict(list(totals.groupby('var')))
    anntot = pd.merge(g['FPSN'], g['RAIN'], on=['year', 'case', 'loc'],
                      suffixes=g.keys())
    anntot.loc[anntot['case'] == 'IDE', 'case'] = "drought"
    anntot.loc[anntot['case'] == 'CTL', 'case'] = "control"
    anntot.loc[anntot['loc'] == ('Loma Ridge '
                                 'Global Change'
                                 ' Experiment'), 'loc'] = "Loma Ridge"
    anntot.loc[anntot['loc'] == 'McLaughlin NRS', 'loc'] = "McLaughlin"
    anntot.loc[anntot['loc'] == 'Sedgewick NRS', 'loc'] = "Sedgwick"
    anntot.loc[anntot['loc'] == ('Sierra Foothill'
                                 ' Research Extension'
                                 ' Center'), 'loc'] = "Sierra Foothill"
    anntot.drop(['varFPSN', 'varRAIN'], axis=1, inplace=True)
    anntot = anntot.rename(columns={'locRAIN': 'loc'})
    # reorder from wet to dry
    lomaridge_gpp = LomaRidgeTools.annual_total_main(
        os.path.join('/', 'project', 'projectdirs', 'm2319', 'Data',
                     'LomaRidgeGlobalChangeExperiment', 'Grass_v3_4.mat'))

    sns.set_context("talk")
    with sns.axes_style("white"):
        g = sns.FacetGrid(anntot, col='loc', hue='case',
                          col_wrap=3,
                          col_order=locs,
                          margin_titles=True, size=6,
                          hue_kws={"marker": ["^", "v"]})
    g.map(plt.scatter, 'annual_RAIN', 'annual_FPSN', s=100)
    g.set_titles(template='{col_name}')
    g.set_axis_labels(y_var=r'annual FPSN (g C m$^{{-2}}$ yr$^{{-1}}$)',
                      x_var='annual rain (mm)')
    g.add_legend()
    loma_idx = ["Loma" in x for x in locs].index(True)
    g.axes[loma_idx].scatter(lomaridge_gpp['RAIN'],
                             lomaridge_gpp['GPP_gC_GRASS'],
                             marker='x', s=80,
                             label='VPRM RE - obs NEE')
    g.axes[loma_idx].legend()
    plt.savefig(os.path.join(os.getenv('CSCRATCH'), 'plots',
                             './rain_vs_pcp.pdf'))
    plt.close()
    return site_data, totals


def plot_CLM_variable(df, varname, ann_diff=None, locs=None):
    """plot monthly time series, monthly point/whisker plots

    ARGS:
    df (pandas.DataFrame): containing monthly variable values, units, variable
       long name
    varname (string): CLM variable name (i.e. the short, all-caps name)
    """
    sns.set_context("talk")
    with sns.axes_style("white"):
        g = sns.factorplot(data=df, x='month', y='value',
                           col='loc', hue='case',
                           kind='point', col_wrap=3,
                           legend=False, legend_out=False,
                           col_order=locs,
                           palette="Dark2",
                           markers=["o", "x"])
        format_factorgrid(g,
                          df['var'].iloc[0],
                          df['var_lname'].iloc[0],
                          df['units'].iloc[0],
                          'month',
                          ann_diff)
        g.savefig(os.path.join(os.getenv('CSCRATCH'), 'plots',
                               '{var}_bymonth.pdf'.format(var=varname)))
        plt.close(g.fig)
    with sns.axes_style("white"):
        g = sns.FacetGrid(df, hue='case', palette="Dark2",
                          col='loc', col_wrap=3,
                          size=3, aspect=1.5,
                          hue_kws={"marker": ["o", "x"],
                                   "linestyle": ['-', '-']},
                          col_order=locs,
                          legend_out=False)
        g.map(plt.plot, 'fyear', 'value')
        g.set_titles(template='{col_name}')
        g = format_factorgrid(g,
                              df['var'].iloc[0],
                              df['var_lname'].iloc[0],
                              df['units'].iloc[0],
                              'year of simulation')
        g.savefig(os.path.join(os.getenv('CSCRATCH'), 'plots',
                               '{var}_timeseries.pdf'.format(var=varname)))
        plt.close(g.fig)


def format_factorgrid(g, var_sname, var_lname, var_units, x_var_name,
                      ann_diff=None):
    fpsn_units_raw = "umol/m2s"
    fpsn_units = "$\mu$mol m$^{{-2}}$ s$^{{-1}}$"
    g.set_axis_labels(x_var=x_var_name,
                      y_var='{var} ({units})'.format(
                          var=var_sname,
                          units=var_units.replace(fpsn_units_raw, fpsn_units)))
    g.set_titles(template='{col_name}')
    new_labels = {'IDE': 'drought', 'CTL': 'control'}
    for k in g._legend_data.keys():
        g._legend_data[new_labels[k]] = g._legend_data.pop(k)
    g.add_legend()
    plt.figtext(0.5, 0.99,
                '{shortname}: {longname}'.format(
                    shortname=var_sname,
                    longname=var_lname),
                horizontalalignment='center')
    if ann_diff is not None:
        for this_site, this_ax in zip(ann_diff.index, g.axes):
            this_ax.annotate(s=(r'$\Delta${v}: {gC:0.0f} '
                                'gC m$^{{-2}}$ yr$^{{-1}}$'
                                ' ({pct:0.0f}\%)').format(
                                    v=var_sname,
                                    gC=ann_diff.ix[this_site].d_gC,
                                    pct=ann_diff.ix[this_site].pct),
                             xy=(0.5, 0.01),
                             xycoords="axes fraction",
                             ha='center')
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
    ide = mean_var.ix['IDE']
    ctl = mean_var.ix['CTL']
    dvar = ide.join(ctl, lsuffix='IDE', rsuffix='CTL')
    dvar['d'] = dvar['valueCTL'] - dvar['valueIDE']
    dvar['ndays'] = map(lambda x: calendar.monthrange(2001, x)[1],
                        dvar.index.get_level_values('month'))
    dvar['gC_IDE'] = carbon_umol_m2_s_2_g_m2_yr(dvar['valueIDE'],
                                                dvar['ndays'])
    dvar['gC_CTL'] = carbon_umol_m2_s_2_g_m2_yr(dvar['valueCTL'],
                                                dvar['ndays'])
    dvar['d_gC'] = dvar['gC_CTL'] - dvar['gC_IDE']
    annual_total = dvar.groupby(level='loc').sum()
    annual_total.eval('pct = (-100 * (gC_CTL - gC_IDE) / gC_CTL)',
                      inplace=True)
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
                      'Sierra Foothill R.E.C.',
                      'Younger Lagoon',
                      'Sedgewick NRS',
                      'Loma Ridge G.C.E.',
                      'Box Springs',
                      'Mammoth Lakes']
    # cal_wet_to_dry = ['McLaughlin NRS',
    #                   'Sierra Foothill Research Extension Center',
    #                   'Younger Lagoon',
    #                   'Sedgewick NRS',
    #                   'Loma Ridge Global Change Experiment',
    #                   'Box Springs']
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

    if False:
        # all_vars = all_vars.loc[all_vars['loc'] != "Mammoth Lakes", :].copy()
        sys.stdout.write('plotting ')
        for v in (h1vars + h2vars):  # ('FPSN', ):
            sys.stdout.write('{} '.format(v))
            sys.stdout.flush()
            df = all_vars[all_vars['var'] == v]
            mon_diff, ann_diff = calc_dvar(all_vars, v)
            ann_diff = ann_diff.reindex([x.name for x in cal_locs])
            if v != "FPSN":
                ann_diff = None
            plot_CLM_variable(df, v, ann_diff,
                              [s.name for s in cal_wet_to_dry])
        sys.stdout.write('\n')
        sys.stdout.flush()

    order = ['McLaughlin', 'Sierra Foothill', 'Younger Lagoon',
             'Sedgwick', 'Loma Ridge', 'Box Springs', 'Mammoth Lakes']
    site_data, anntot_long = plot_site_annual_rain_gpp(
        all_vars, order)

    # fpsn = all_vars.ix[all_vars['var'] == 'FPSN']
    # ndays = map(lambda x: calendar.monthrange(2001, x)[1],
    #             fpsn['month'].values)
    # fpsn.loc[:, 'gC'] = np.vectorize(carbon_umol_m2_s_2_g_m2_yr)(fpsn['value'],
    #                                                              ndays)
    # bar = fpsn.join(fpsn.loc[:, ('gC', 'fyear')].cumsum(), rsuffix='sum')
    # g = sns.factorplot(data=bar, x='fyear', y='gCsum',
    #                    col='loc', hue='case',
    #                    kind='point', col_wrap=3,
    #                    legend=False, legend_out=False,
    #                    palette="Dark2",
    #                    markers=["o", "x"])
