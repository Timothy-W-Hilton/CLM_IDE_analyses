"""plot CLM LE, ZWT, TWS against PCP

LE: latent heat flux
ZWT: total water table depth
TWS: total water storage
PCP: precipitation
"""
import matplotlib
matplotlib.use('AGG')

import os
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from datetime import datetime
from RegionTester.region_tester import InUSState
from clm_domain import CLM_Domain

class CalMask(object):
    """mask to determine if a lat/lon points are within California
    """
    def __init__(self, lon, lat):
        """populate self.lon, self.lat

        ARGS:
            lon (array-like): 2D array with shape (nlon, nlat) of
               longitude values in the range (-180, 180)
            lat (array-like): 2D array with shape (nlon, nlat) of
               latitude values in the range (-90, 90)
           """
        self.lon = lon
        self.lat = lat

    def mask(self):
        California = InUSState('')
        California.get_state_shape('California')
        iscal = np.empty(self.lon.shape, dtype=bool)
        iscal[...] = False
        for (x, y), this_lon in np.ndenumerate(lon):
            iscal[x, y] = California.point_inside(lon[x, y], lat[x, y])
        return iscal


class USAMask(object):
    """mask to determine if a lat/lon points are within the USA

    currently unimplemented
    """
    def __init__(self):
        warnings.warn('not yet implemented')


class QianMonthlyPCPData(object):
    """
    read total monthly PCP for Qian et al (2006) atmosphere

    The monthly total values are calculated by QianTotaller
    """

    def __init__(self, fname):
        """populate field monthly_pcp

        ARGS:
            fname (string): full path to file containing monthly PCP
                from QianTotaller """
        self.fname = fname

    def read_nc(self):
        """read pcp data from netcdf file
        """
        nc = netCDF4.Dataset(self.fname, 'r')
        self.pcp_all = nc.variables['PRECTmms'][:]
        self.lon = nc.variables['LONGXY'][...]
        self.lat = nc.variables['LATIXY'][...]
        nc.close()

    def extract(self, mask):
        """extract pixels of interest from Qian data

        ARGS:
            mask (arraylike): an optional mask.  If supplied, it must
                have shape (nlon, lat) matching the size of the pcp
                grid.  Must containin boolean values: where True the
                PCP data are to be read, where False the PCP data are
                ignored. """
        if mask is None:
            x = slice(None)
            y = slice(None)
        else:
            x, y = np.where(mask)

    def interpolate(self, dlon, dlat):
        # verify that lon does not change with latitude, and lat does
        # not change with longitue
        lon_not_regular = any([np.diff(self.lon[:, i]).any()
                               for i in range(self.lon.shape[1])])
        lat_not_regular = any([np.diff(self.lat[i, :]).any()
                               for i in range(self.lat.shape[0])])
        if lat_not_regular or lon_not_regular:
            raise ValueError('latitude and longitude must be a regular grid')

        ntime = self.pcp_all.shape[0]
        nlon = dlon.shape[0]
        nlat = dlat.shape[1]
        self.pcp = np.zeros((ntime, nlon, nlat))

        for t in np.arange(ntime)[:5]:
            if np.mod(t, 100) == 0:
                print t, datetime.now()
            finterp = interpolate.RectBivariateSpline(self.lat[:, 0],
                                                      self.lon[0, :],
                                                      self.pcp_all[t, ...])
            self.pcp[t, ...] = finterp.ev(dlon, dlat)


def check_results(qmd, dlon, dlat):
    fig, ax = plt.subplots(2, 1, figsize=(8.5, 11))
    cmap=plt.get_cmap('Blues')
    cm = ax[0].pcolormesh(qmd.lon, qmd.lat, qmd.pcp_all[0, ...],
                          cmap=cmap)
    for this_ax in ax:
        this_ax.set_xlabel('lon E')
        this_ax.set_ylabel('lat N')
    cb = plt.colorbar(cm, ax=ax[0])
    cm = ax[1].pcolormesh(dlon, dlat, qmd.pcp[0, ...],
                          cmap=cmap)
    cb = plt.colorbar(cm, ax=ax[1])
    fig.savefig(os.path.join(os.getenv('HOME'), 'plots',
                             'qian_pcpinterp_test.png'))
    plt.close(fig)


if __name__ == "__main__":
    pcp_ncfile = os.path.join(os.getenv('SCRATCH'),
                              'qian_pcp_monthly_totals.nc')
    qmd = QianMonthlyPCPData(pcp_ncfile)
    qmd.read_nc()
    d = CLM_Domain(fname=os.path.join('/', 'global',
                                      'cscratch1', 'sd',
                                      'twhilton',
                                      'archive',
                                      'CLM_f05_g16',
                                      'lnd',
                                      'hist',
                                      'CLM_f05_g16.clm2.h0.0050-01.nc'))
    qmd.interpolate(d.get_lon(), d.get_lat())
    check_results(qmd, d.get_lon(), d.get_lat())
