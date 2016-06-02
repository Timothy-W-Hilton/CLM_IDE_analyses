"""tools to check up on CLM spinups.  Currently it is mostly geared
toward parsing time series of a few key variables and making sure the
spinup achieved steady state """

import netCDF4
import numpy as np
import os
import glob
import matplotlib.pyplot as plt
import calendar

from clm_domain import CLM_Domain, Location


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
        mask_missing ({True}|False): if True, check the requested
            netCDF variable for an attribute named "missing_value".
            If present, parse it and return a numpy masked array with
            values equal (within floating point tolerance) to
            missing_value masked.
        verbose (True|{False}): if True, display a status message as
            each file is read

        It is assumed that CLM output files contain one and only one
        timestep.
        """

        nt = len(self.all_files)
        nlon = 1
        nlat = 1
        if (lon_idx is None) or (lat_idx is None):
            nc = netCDF4.Dataset(self.all_files[0])
            nlon = len(nc.dimensions['lon'])
            nlat = len(nc.dimensions['lat'])
            lon_idx = np.arange(nlon)
            lat_idx = np.arange(nlat)
            nc.close()
        data = np.zeros([nt, nlat, nlon])

        # TODO: intialize array for time here

        missing_value = None
        for i in range(nt):
            if verbose:
                print 'parsing {}'.format(os.path.basename(self.all_files[i]))
            nc = netCDF4.Dataset(self.all_files[i])
            data[i, ...] = nc.variables[varname][0, lat_idx, lon_idx]
            # TODO: parse time here
            try:
                missing_value = nc.variables[varname].missing_value
            except AttributeError:
                if verbose:
                    print "no missing_value for {} in {}".format(
                        os.path.basename(self.all_files[i]), varname)
            nc.close()
        if missing_value is not None:
            data = np.ma.masked_values(data, missing_value)
        return data


def parse_CLM_f05_g15():
    """wrapper function to read data from CLM_f05_g16 spinup run
    """
    CLM_f05_g16 = CLM_spinup_analyzer(os.path.join('/', 'global',
                                                   'cscratch1', 'sd',
                                                   'twhilton',
                                                   'archive',
                                                   'CLM_f05_g16',
                                                   'lnd',
                                                   'hist'),
                                      'CLM_f05_g16')
    fpsn = CLM_f05_g16.parse_var('FPSN', verbose=True)
    return (CLM_f05_g16, fpsn)


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


def plot_monthly_timeseries(spinup_run, data, santacruz):
    months = range(1, 13)
    fig, ax = plt.subplots(ncols=2, nrows=np.int(np.ceil(len(months) / 2)))
    ax[0, 0].set_title('Santa Cruz GPP (umol m-2 s-1)')
    for i, this_month in enumerate(months):
        ax_idx = np.unravel_index(i, ax.shape)
        idx = np.array([this_file.find('{:02d}.nc'.format(this_month))
                        for this_file in spinup_run.all_files]) > 0
        ax[ax_idx].plot(fpsn[idx, santacruz.clm_y, santacruz.clm_x])
        ax[ax_idx].set_ylabel('GPP')
        ax[ax_idx].xaxis.set_major_formatter(plt.NullFormatter())
        ax[ax_idx].annotate(calendar.month_abbr[this_month],
                            xycoords='axes fraction',
                            xy=(0, 0),
                            bbox=dict(boxstyle="round", fc="0.8"))
    ax[ax_idx].xaxis.set_major_formatter(plt.ScalarFormatter())
    ax[ax_idx].set_xlabel('year of spinup')
    return(fig, ax)

if __name__ == "__main__":

    CLM_f05_g16, fpsn = parse_CLM_f05_g15()
    domain_f05_g16 = CLM_Domain(CLM_f05_g16.all_files[0])
    santacruz = Location((-122.03089741, ), (36.9741, ))
    santacruz.clm_y, santacruz.clm_x = domain_f05_g16.find_nearest_xy(
        santacruz.lon,
        santacruz.lat)
    plot_monthly_timeseries(CLM_f05_g16, fpsn, santacruz)
