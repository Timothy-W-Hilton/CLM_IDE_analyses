"""tools to check up on CLM spinups.  Currently it is mostly geared
toward parsing time series of a few key variables and making sure the
spinup achieved steady state """

import netCDF4
import numpy as np
import os
import glob
import calendar
import matplotlib.pyplot as plt
import pandas as pd

from clm_domain import CLM_Domain, Location


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
        self.tstamp = tstamp

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
        df = pd.DataFrame({'tstamp': self.tstamp,
                           self.varname: self.data.squeeze()})
        # I'm not using np.nanmean because I don't think CLM should
        # produce any Nans, so I want it to throw an error if it
        # encounters one
        am = df.groupby([t.year for t in df.tstamp]).aggregate(np.mean)
        return am


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
        tstamp = pd.Series(np.zeros([nt]))

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
            this_tstamp_date = pd.to_datetime(
                np.datetime64("{:04d}-{:02d}-{:02d}".format(t_year,
                                                            t_month,
                                                            t_day)))
            tstamp[i] = this_tstamp_date + pd.Timedelta(t_secs, 's')
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


class AnnualMeanPlotter(object):
    """produce 6-panel plot of CLM variable annual means for spinup diagnosis

    The variables plotted are sensible heat (FSH), latent heat (QSOIL
    + QVEGE + QVEGT), liquid runoff (QRUNOFF), water table depth
    (ZWT), total leaf area index (TLAI), and GPP (FPSN).  CLM does not
    report latent heat flux as such so I use [soil evaporation (QSOIL)
    + canopy evaporation (QVEGE) + canopy transpiration (QVEGT)].
    """

    def __init__(self, spinup_container, location):
        """populate fields spinup, location, var_list
        """
        self.location = location
        self.spinup = spinup_container
        self.var_list = ['FSH', 'QSOIL', 'QVEGE', 'QVEGT',
                         'QRUNOFF', 'ZWT', 'TLAI', 'FPSN']
    def get_data(self):
        """read FSH, QSOIL, QVEGE, QVEGT, QRUNOFF, ZWT, TLAI, and FPSN
        variables from the spinup run and calculate annual means.
        Populates class fields t with timestamps and fields FSH,
        QSOIL, QVEGE, QVEGT, QRUNOFF, ZWT, TLAI, and FPSN with annual
        mean values for those variables.
        """

        for this_var in self.var_list:
            soil_lev = None  # none of these are soil variables
            var_obj = self.spinup.parse_var(this_var,
                                            lon_idx=self.location.clm_x,
                                            lat_idx=self.location.clm_y,
                                            soil_lev=soil_lev)
            setattr(self, this_var, var_obj)
        # populate latent heat
        self.Q = CLM_var(varname='Q',
                         varname_long='latent heat flux',
                         data=(self.QSOIL.data +
                               self.QVEGE.data +
                               self.QVEGT.data),
                         units=self.QSOIL.units,
                         missing_value=self.QSOIL.missing_value)

    def _setup_plot(self):
        self.fig, self.ax = plt.subplots(nrows=3, ncols=2, figsize=(8.5, 11))
        self.fig.text(0.5, 0.95, '{}'.format(self.location.name),
                      verticalalignment='top',
                      horizontalalignment='center',
                      size='large')

    def plot(self):
        """draw the six panel plot
        """
        self._setup_plot()
        am = self.FSH.annual_mean()
        self.ax[0, 0].plot(am.index, am.values)
        self.ax[0, 0].set_title(self.FSH.varname_long)
        self.ax[0, 0].set_ylabel(self.FSH.units)
        # self.ax[0, 1].plot(self.t, self.QSOIL + self.QVEGE + self.QVEGT)
        # self.ax[1, 0].plot(self.t, self.QRUNOFF)
        # self.ax[1, 1].plot(self.t, self.ZWT)
        # self.ax[2, 0].plot(self.t, self.TLAI)
        # self.ax[2, 1].plot(self.t, self.FPSN)



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
    return (CLM_f05_g16, santacruz, mclaughlin)


def plot_CLMf05g16_monthly_timeseries_main():
    """plot FPSN', 'H2OSOI', 'H2OSOI', 'ZWT', 'TLAI' at Santa Cruz, McLaughlin
    """
    (CLM_f05_g16, santacruz, mclaughlin) = CLMf05g16_get_spatial_info()

    for this_loc in (santacruz, mclaughlin):
        varname = ['FPSN', 'H2OSOI', 'H2OSOI', 'ZWT', 'TLAI']
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
    (CLM_f05_g16, santacruz, mclaughlin) = CLMf05g16_get_spatial_info()
    var = parse_CLM_f05_g15(CLM_f05_g16,
                            'FPSN',
                            None,
                            location=santacruz)
    return var
if __name__ == "__main__":
    # plot_CLMf05g16_monthly_timeseries_main()
    CLM_f05_g16, santacruz, mclaughlin = CLMf05g16_get_spatial_info()
    amp = AnnualMeanPlotter(CLM_f05_g16, santacruz)
    amp.get_data()
    amp.plot()
