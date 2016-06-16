"""tools to check up on CLM spinups.  Currently it is mostly geared
toward parsing time series of a few key variables and making sure the
spinup achieved steady state """

import matplotlib
matplotlib.use('AGG')

import netCDF4
import numpy as np
import os
import glob
import calendar
import matplotlib.pyplot as plt
import pandas as pd

from clm_domain import CLM_Domain, Location


class ARM_data(object):
    """container class for ARM Southern Great Plains Ameriflux data
    """
    def __init__(self, fname):
        self.fname = fname
        self.data = None

    def parse_data(self):
        self.data = pd.read_csv(self.fname, header=2, na_values=('-9999', ))
        self.data = self.data[['TIMESTAMP_START', 'H', 'LE']]
        #extract year from date in format YYYYMMDDHHMM
        self.data['YEAR'] = np.floor(self.data['TIMESTAMP_START'] / 1e8)

    def annual_mean(self):
        self.am = self.data.groupby('YEAR').aggregate(np.mean)


class CLM_var(object):
    """container class for a CLM output variable data and metadata
    """

    def __init__(self,
                 varname="",
                 varname_long=None,
                 data=None,
                 time=None,
                 units="",
                 depth=None,
                 missing_value=1.0e36):
        """
        ARGS:
        varname (string): the name of the variable in the netCDF output
        varname_long (string): name of variable suitable for plotting
        data (array-like): the data values
        time (array-like): variable time stamps
        units (string): the units for the variables
        depth (float): soil depth in meters (if applicable)
        """
        self.varname = varname
        self.varname_long = varname_long
        self.data = data
        self.time = time
        self.units = units
        self.depth = depth
        self.missing_value = missing_value

    def get_long_name(self):
        """returns long name if provided; short name otherwise
        """
        if self.varname_long is not None:
            return self.varname_long
        else:
            return self.varname

    def set_data(self, data):
        self.data = data

    def set_tstamp(self, tstamp):
        self.time = tstamp

    def get_plot_filename(self, location):
        """generate a filename for a diagnostic plot

        ARGS:
        location: a Location instance or a string

        constructs a file name containing the variable short name,
        depth (if applicable), and location
        """
        var_str = self.varname
        if self.depth is not None:
            var_str = var_str + '_{:1.2f}m'.format(self.depth)
        fname = 'spinup_{}_{}.pdf'.format(var_str,
                                          location.name.replace(" ", ""))
        return fname

    def annual_mean(self):
        """return pandas dataframe containing the annual mean
        """
        if self.data.squeeze().ndim != 1:
            raise ValueError(('currently can only calculate annual mean '
                              'for a scalar time series.'))
        df = pd.DataFrame({self.varname: self.data.squeeze()})
        # I'm not using np.nanmean because I don't think CLM should
        # produce any Nans, so I want it to throw an error if it
        # encounters one
        #
        # extract a year, month, day from a numpy datetime64 by
        #  (dt64.astype(object).year)
        yr = [t.astype(object).year for t in self.time]
        am = df.groupby(yr).aggregate(np.mean)
        return am

    def monthly_mean(self):
        """return pandas dataframe containing the annual mean
        """
        if self.data.squeeze().ndim != 1:
            raise ValueError(('currently can only calculate annual mean '
                              'for a scalar time series.'))
        df = pd.DataFrame({'tstamp': self.time,
                           self.varname: self.data.squeeze()})
        # I'm not using np.nanmean because I don't think CLM should
        # produce any Nans, so I want it to throw an error if it
        # encounters one
        month = [t.month for t in df.tstamp]
        mm = df.groupby(month).aggregate(np.mean)
        return mm


class CLM_spinup_analyzer(object):
    """class to read CLM spinup output and plot key variables
    """

    def __init__(self,
                 data_dir,
                 CASE):
        """
        Args:
        data_dir: full path to the directory containing the spinup run data
        CASE (string): the "case" name of the CESM run
        """
        self.data_dir = data_dir
        self.CASE = CASE
        self.gather_filenames()

    def gather_filenames(self):
        """gather all netCDF (*.nc) filenames from self.data_dir and place
        them in self.data_filenames
        """
        self.all_files = sorted(glob.glob(os.path.join(self.data_dir, "*.nc")))

    def build_fname(self, year, month):
        """assemble the filename for the
        """
        return "{}.clm2.h0.{:04d}-{:02d}.nc".format(self.CASE, year, month)

    def parse_var(self,
                  varname,
                  lon_idx=None,
                  lat_idx=None,
                  soil_lev=None,
                  mask_missing=True,
                  verbose=False):
        """parse a specified variable from all spinup files.

        Optionally, restrict parsing to a specified grid cell.

        Args:
        varname (string): name of the netCDF variable to parse.  Must
        be present in all files in self.data_filenames
        lon_idx (integer): optional; if specified only a single grid
        cell will be parsed.  Ignored if lat_idx is not also specified.
        lat_idx (integer): optional; if specified only a single grid
        cell will be parsed
        soil_lev (integer): optional, if specified, indicates that the
        variable is a soil variable and a soil level is needed
        mask_missing ({True}|False): if True, check the requested
        netCDF variable for an attribute named "missing_value".
        If present, parse it and return a numpy masked array with
        values equal (within floating point tolerance) to
        missing_value masked.
        verbose (True|{False}): if True, display a status message as
        each file is read

        RETURNS:
        a CLM_var instance containing the parsed data

        It is assumed that CLM output files contain one and only one
        timestep.
        """

        nt = len(self.all_files)
        nlon = 1
        nlat = 1
        nc = netCDF4.Dataset(self.all_files[0])
        soil_depth = None
        if soil_lev is not None:
            soil_depth = nc.variables['levgrnd'][soil_lev]

        var = CLM_var(varname,
                      varname_long=nc.variables[varname].long_name,
                      units=nc.variables[varname].units,
                      time=None,  # ignore for now
                      depth=soil_depth,
                      missing_value=nc.variables[varname].missing_value)

        if (lon_idx is None) or (lat_idx is None):
            nlon = len(nc.dimensions['lon'])
            nlat = len(nc.dimensions['lat'])
            lon_idx = np.arange(nlon)
            lat_idx = np.arange(nlat)
            nc.close()
        data = np.zeros([nt, nlat, nlon])
        tstamp = np.empty(nt, dtype='datetime64[s]')

        missing_value = None
        for i in range(nt):
            if verbose:
                print 'parsing {} from {}'.format(
                    varname,
                    os.path.basename(self.all_files[i]))
            nc = netCDF4.Dataset(self.all_files[i])
            # parse data
            if soil_lev is None:
                # each file contains one timestep so the time index is
                # always 0
                data[i, ...] = nc.variables[varname][0, lat_idx, lon_idx]
            else:
                data[i, ...] = nc.variables[varname][0, soil_lev,
                                                     lat_idx, lon_idx]
            # parse time stamp
            t_YYYYMMDD_str = str(nc.variables['mcdate'][0])
            t_year = np.int(t_YYYYMMDD_str[0:-4])
            t_month = np.int(t_YYYYMMDD_str[-4:-2])
            t_day = np.int(t_YYYYMMDD_str[-2:])
            t_secs = np.int(nc.variables['mcsec'][0])
            # TODO: update for numpy 1.11 "timezone-naive" datetime64
            # http://docs.scipy.org/doc/numpy/reference/arrays.datetime.html#changes-with-numpy-1-11
            this_tstamp_date = np.datetime64(
                "{:04d}-{:02d}-{:02d}T00:00".format(t_year,
                                                    t_month,
                                                    t_day))
            tstamp[i] = this_tstamp_date + np.timedelta64(t_secs, 's')
            # parse missing value
            try:
                missing_value = nc.variables[varname].missing_value
            except AttributeError:
                if verbose:
                    print "no missing_value for {} in {}".format(
                        os.path.basename(self.all_files[i]), varname)
            nc.close()
        if missing_value is not None:
            data = np.ma.masked_values(data, missing_value)
        var.set_data(data)
        var.set_tstamp(tstamp)
        return var


class AnnualCyclePlotter(object):
    """produce 6-panel plot of CLM variable annual cycle for spinup diagnosis

    can plot either monthly mean or annual mean time series

    The variables plotted are sensible heat (FSH), latent heat (QSOIL
    + QVEGE + QVEGT), liquid runoff (QRUNOFF), total water storage
    (WT), total leaf area index (TLAI), and GPP (FPSN).  CLM does not
    report latent heat flux as such so I use [soil evaporation (QSOIL)
    + canopy evaporation (QVEGE) + canopy transpiration (QVEGT)].
    """

    def __init__(self, spinup_container, location):
        """populate fields spinup, location, var_list
        """
        self.location = location
        self.spinup = spinup_container
        self.var_list_parse = ['FSH', 'QSOIL', 'QVEGE', 'QVEGT',
                               'QRUNOFF', 'WT', 'TLAI', 'FPSN']
        self.var_list_plot = ['FSH', 'LE', 'QRUNOFF',
                              'WT', 'TLAI', 'FPSN']


    def get_ARM_data(self):
        """calculate Ameriflux observed annual mean from US-ARM
        """
        usarm = ARM_data('./AMF_US-ARM_BASE_HH_6-1.csv')
        usarm.parse_data()
        usarm.annual_mean()
        nyears = usarm.am.shape[0]
        usarm.am.index = np.arange(51 - nyears, 51)
        return usarm.am


    def get_data(self):
        """read FSH, QSOIL, QVEGE, QVEGT, QRUNOFF, WT, TLAI, and FPSN
        variables from the spinup run and calculate annual means.
        Populates class fields t with timestamps and fields FSH,
        QSOIL, QVEGE, QVEGT, QRUNOFF, WT, TLAI, and FPSN with annual
        mean values for those variables.
        """

        for this_var in self.var_list_parse:
            soil_lev = None  # none of these are soil variables
            var_obj = self.spinup.parse_var(this_var,
                                            lon_idx=self.location.clm_x,
                                            lat_idx=self.location.clm_y,
                                            soil_lev=soil_lev)
            setattr(self, this_var, var_obj)

        C_g_per_mol = 12.01
        mol_per_umol = 1e-6
        s_per_year = 60 * 60 * 24 * 365
        self.FPSN.data = (self.FPSN.data *
                          C_g_per_mol * mol_per_umol * s_per_year)
        self.FPSN.units = 'g C m-2 yr-1'

        # populate latent heat
        self.LE = CLM_var(varname='LE',
                          varname_long='latent heat flux',
                          data=mm_s_2_w_m2_s(self.QSOIL.data +
                                             self.QVEGE.data +
                                             self.QVEGT.data),
                          time=self.QSOIL.time,
                          units='W/m^2',
                          missing_value=self.QSOIL.missing_value)
        s_per_day = 24.0 * 60.0 * 60.0
        self.QRUNOFF.data = self.QRUNOFF.data * s_per_day
        self.QRUNOFF.units = 'mm/day'

    def _setup_plot(self):
        self.fig, self.ax = plt.subplots(nrows=3, ncols=2, figsize=(8.5, 11))
        self.fig.text(0.5, 0.95, '{}'.format(self.location.name),
                      verticalalignment='top',
                      horizontalalignment='center',
                      size='large')

    def plot(self, cycle='annual'):
        """draw the six panel plot

        ARGS:
        cycle ({"annual"} | "monthly"): if annual (default), plots a
        time series of annual means.  If "monthly", plots the monthly
        means calulated across the entire time series.
        """
        lw = 2.0  # line width
        self._setup_plot()
        for this_var, this_ax in zip(self.var_list_plot, self.ax.flatten()):
            if cycle is "monthly":
                mean_vals = getattr(self, this_var).monthly_mean()
                label_str = 'monthly mean'
            elif cycle is "annual":
                # the last year of the spinup only went through
                # October so year 51 is not a full year annual average
                mean_vals = getattr(self, this_var).annual_mean()[:50]
                label_str = 'annual mean '
            else:
                raise ValueError('cycle must be "annual" or "monthly"')
            this_ax.plot(mean_vals.index, mean_vals.values, label=label_str,
                         linewidth=lw)
            if cycle is "annual":
                this_ax.plot(mean_vals.index,
                             pd.rolling_mean(mean_vals.values, window=10),
                             '--', label='10-year running mean',
                             linewidth=lw)
                if self.location.name == "ARM Southern Great Plains":
                    armobs = self.get_ARM_data()
                    self.ax[0, 0].plot(armobs.index, armobs.H, 'k:',
                                       label='observed',
                                       linewidth=lw)
                    self.ax[0, 1].plot(armobs.index, armobs.LE, 'k:',
                                       linewidth=lw)
                    # put a dummy line in the axes from which the
                    # legend is drawn
                    dummy = self.ax[-1, 0].plot([], [], 'k:', label='observed',
                                                linewidth=lw)
            # don't place parenthetical part of long name in plot title
            t_str = getattr(self, this_var).varname_long.split("(")[0]
            this_ax.set_title(t_str)
            this_ax.set_ylabel(getattr(self, this_var).units)
        self.ax[-1, 0].legend(loc='upper left',
                              bbox_to_anchor=(0.2, -0.15),
                              ncol=2)
        for this_col in range(2):
            self.ax[-1, this_col].set_xlabel('month')


def mm_s_2_w_m2_s(mm_s):
    """convert latent heat flux (or evapotranspiration) from units of
    mm s-1 to W m-2 s-1
    """
    # 1 mm m-2 s-1 equals 1 kg m-2 s-1
    # (http://www.colorado.edu/geography/class_homepages/geog_3511_s11/notes/Notes_9.pdf)
    #
    # LE = E * lambda_v
    #
    # lambda_v is latent heat of evaporation of water = 2.501x1e6 J kg-1
    #
    # therefore LE = Q (mm m-2 s-1) * lambda_v
    lambda_v = 2.501*1e6
    return mm_s * lambda_v


def parse_CLM_f05_g15(spinup_container, varname, soil_lev=None, location=None):
    """wrapper function to read data from CLM_f05_g16 spinup run
    """
    if location is None:
        clm_x = None
        clm_y = None
    else:
        clm_x = location.clm_x
        clm_y = location.clm_y
    data = spinup_container.parse_var(varname,
                                      soil_lev=soil_lev,
                                      lon_idx=clm_x,
                                      lat_idx=clm_y,
                                      verbose=True)
    return data


def plot_laugh_test_map(fpsn, santacruz):
    """plot one month of spinup run GPP and Santa Cruz, CA

    This is a quick 'laugh test' just to make sure the GPP spatial
    placement actually looks like the continents and that the grid
    indices of Santa Cruz are correct
    """
    fig, ax = plt.subplots()
    ax.pcolormesh(fpsn[100, ...], cmap=plt.get_cmap('Greens'))
    ax.scatter(santacruz.clm_x, santacruz.clm_y, c='black', marker='x', s=40)
    return fig, ax


def plot_monthly_timeseries(spinup_run, var, location):
    """12-panel plot: timeseries of each month from all spinup years

    ARGS:
    spinup_run (CLM_spinup_analyzer): object of class
    CLM_spinup_analyzer describing the spinup run
    var (CLM_var): object of class CLM_var containing the data and
    metadata to be plotted
    location (Location): object of class Location describing the
    location of the data to be plotted
    """
    months = range(1, 13)
    nrows = np.int(np.ceil(len(months) / 2))
    fig, ax = plt.subplots(ncols=2, nrows=nrows, figsize=(8.5, 11))

    if var.depth is not None:
        depth_str = "{:0.2f} m".format(var.depth)
    else:
        depth_str = ""
    fig.text(0.5, 0.95, '{} {} {} ({})'.format(location.name,
                                               var.varname_long,
                                               depth_str,
                                               var.units),
             verticalalignment='top',
             horizontalalignment='center',
             size='large')
    for i, this_month in enumerate(months):
        ax_idx = np.unravel_index(i, ax.shape)
        idx = np.array([this_file.find('{:02d}.nc'.format(this_month))
                        for this_file in spinup_run.all_files]) > 0
        ax[ax_idx].plot(var.data[idx, ...].squeeze())
        ax[ax_idx].set_ylabel(var.varname)
        ax[ax_idx].xaxis.set_major_formatter(plt.NullFormatter())
        ax[ax_idx].annotate(calendar.month_abbr[this_month],
                            xycoords='axes fraction',
                            xy=(0, 0),
                            bbox=dict(boxstyle="round", fc="0.8"))
    for this_col in range(2):
        ax[nrows-1, this_col].xaxis.set_major_formatter(plt.ScalarFormatter())
        ax[nrows-1, this_col].set_xlabel('year of spinup')
    return(fig, ax)


def CLMf05g16_get_spatial_info():
    CLM_f05_g16 = CLM_spinup_analyzer(os.path.join('/', 'global',
                                                   'cscratch1', 'sd',
                                                   'twhilton',
                                                   'archive',
                                                   'CLM_f05_g16',
                                                   'lnd',
                                                   'hist'),
                                      'CLM_f05_g16')
    domain_f05_g16 = CLM_Domain(CLM_f05_g16.all_files[0])
    santacruz = Location((-122.03089741, ), (36.9741, ), 'Santa Cruz')
    santacruz.get_clm_xy(domain_f05_g16)
    mclaughlin = Location((-122.431667, ), (38.873889, ), 'McLaughlin NRS')
    mclaughlin.get_clm_xy(domain_f05_g16)
    sierra_foothills = Location((-121.311623, ), (39.249407, ),
                                'Sierra Foothill Research Extension Center')
    sierra_foothills.get_clm_xy(domain_f05_g16)
    loma_ridge = Location((-117.7, ), (33.74, ),
                          'Loma Ridge Global Change Experiment')
    loma_ridge.get_clm_xy(domain_f05_g16)
    ARM_SGP = Location((-97.4888, ), (36.6058, ),
                       'ARM Southern Great Plains')
    ARM_SGP.get_clm_xy(domain_f05_g16)
    return (CLM_f05_g16, santacruz, mclaughlin,
            sierra_foothills, loma_ridge, ARM_SGP)


def plot_CLMf05g16_monthly_timeseries_main():
    """plot FPSN', 'H2OSOI', 'H2OSOI', 'WT', 'TLAI' at Santa Cruz, McLaughlin
    """
    (CLM_f05_g16, santacruz, mclaughlin,
     sierra_foothills, loma_ridge, ARM_SGP) = CLMf05g16_get_spatial_info()

    for this_loc in (santacruz, mclaughlin):
        varname = ['FPSN', 'H2OSOI', 'H2OSOI', 'WT', 'TLAI']
        soil_lev = [None, 7, 2, None, None]
        for this_varname, this_lev in zip(varname, soil_lev):
            var = parse_CLM_f05_g15(CLM_f05_g16,
                                    this_varname,
                                    this_lev,
                                    location=this_loc)
            fig, ax = plot_monthly_timeseries(CLM_f05_g16, var, this_loc)
            fig.savefig(
                os.path.join(os.getenv('HOME'),
                             'plots',
                             'CLM_f05_g16',
                             var.get_plot_filename(this_loc)))
            plt.close(fig)


def test_tstamp_parse():
    (CLM_f05_g16, santacruz, mclaughlin,
     sierra_foothills, loma_ridge, ARM_SGP) = CLMf05g16_get_spatial_info()
    var = parse_CLM_f05_g15(CLM_f05_g16,
                            'FPSN',
                            None,
                            location=santacruz)
    return var


def annual_cycle_plots_main():
    (CLM_f05_g16, santacruz, mclaughlin,
     sierra_foothills, loma_ridge, ARM_SGP) = CLMf05g16_get_spatial_info()
    for this_site in (santacruz, mclaughlin,
                      sierra_foothills, loma_ridge, ARM_SGP):
        # for this_site in (ARM_SGP, ):
        print 'plotting summary for {}'.format(this_site.name)
        acp = AnnualCyclePlotter(CLM_f05_g16, this_site)
        acp.get_data()
        for this_cycle in ['annual', 'monthly']:
            acp.plot(this_cycle)
            acp.fig.savefig(os.path.join(os.getenv('HOME'),
                                         'plots',
                                         'CLM_f05_g16',
                                         '{}_spinup_{}_means.pdf'.format(
                        this_site.name.replace(' ', ''),
                        this_cycle)))

if __name__ == "__main__":
    # plot_CLMf05g16_monthly_timeseries_main()
    foo = annual_cycle_plots_main()
