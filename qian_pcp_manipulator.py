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
from numpy import ma
import matplotlib.pyplot as plt
from scipy import interpolate
from datetime import datetime
from mpl_toolkits.basemap import Basemap
import warnings
from timutils import colormap_nlevs
from timutils.io import delete_if_exists
from RegionTester.region_tester import InUSState
from clm_domain import CLM_Domain
from IDE_locations import CLMf05g16_get_spatial_info


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
        for (x, y), this_lon in np.ndenumerate(self.lon):
            iscal[x, y] = California.point_inside(self.lon[x, y],
                                                  self.lat[x, y])
        return np.logical_not(iscal)


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
        self.lon = None
        self.lat = None
        self.pcp_all = None
        self.pcp = None

    def read_nc(self):
        """read pcp data from netcdf file
        """
        nc = netCDF4.Dataset(self.fname, 'r')
        self.pcp_all = nc.variables['PRECTmms'][:] * 21600
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
        """regrid Qian data

        The `Qian et al. (2006) <http://dx.doi.org/10.1175/JHM540.1>`_
        atmospheric data are on a 1 x 2 degree grid.  Use bivariate
        spline approximation over a rectangular mesh to regrid the
        data to a specified lon, lat grid.

        The regridded values are placed in self.pcp.

        PARAMETERS:
        lon, lat: ndarray; of arbitrary longitudes and latitudes.  Must
           contain the same number of elements.

        """
        # verify that lon does not change with latitude, and lat does
        # not change with longitue
        lon_not_regular = any([np.diff(self.lon[:, i]).any()
                               for i in range(self.lon.shape[1])])
        lat_not_regular = any([np.diff(self.lat[i, :]).any()
                               for i in range(self.lat.shape[0])])
        if lat_not_regular or lon_not_regular:
            raise ValueError('latitude and longitude must be a regular grid')

        self.dlon = dlon
        self.dlat = dlat
        ntime = self.pcp_all.shape[0]
        nlon = dlon.shape[0]
        nlat = dlat.shape[1]
        self.pcp = np.zeros((ntime, nlon, nlat))

        for t in np.arange(ntime):
            if np.mod(t, 100) == 0:
                print t, datetime.now()
            finterp = interpolate.RectBivariateSpline(
                self.lon[0, :],
                self.lat[:, 0],
                np.transpose(self.pcp_all[t, ...]))
            self.pcp[t, ...] = finterp.ev(dlon, dlat)
            self.pcp[self.pcp < 0] = 0.0

    def get_IDE_reduction(self):
        """calculate IDE precipitation reduction

        `International Drought Experiment
        (IDE) <https://wp.natsci.colostate.edu/droughtnet/activities/international-drought-experiment/>`_
        protocol is to reduce the precipitation at each site to 1% of
        the long-term average annual precipitation.
        get_IDE_reduction_pct() calculates that redution for each Qian
        pcp data pixel as ((1st percentile) / (50th percentile)), with
        the percentiles calculated from the compete time series of
        annual total pcp for each grid cell.

        RETURNS:
        an NxM numpy array of values in the range [0.0, 1.0].  NxM are
           the number of rows and columns, respectively, in the pcp data.
        """
        pctl = ma.masked_invalid(np.percentile(a=self.pcp, q=(1, 50), axis=0))
        frac = pctl[0, ...] / pctl[1, ...]
        return frac

    def pcp_reduction_to_netcdf(self, outfile):
        """write pcp reduction fractions to netcdf file

        the netcdf file has three variables: lat, lon, frac

        PARAMETERS:
        outfile: string
           full path to the netcdf file to be written.  The file is
           overwritten if it already exists

        Returns
        bool
            True if successful, False otherwise
        """
        if self.pcp is None:
            raise ValueError('pcp not yet parsed')

        frac = self.get_IDE_reduction()

        delete_if_exists(outfile)
        nc = netCDF4.Dataset(outfile, mode='w', format="NETCDF3_CLASSIC")
        nc.createDimension("lat", self.dlat.shape[1])
        nc.createDimension("lon", self.dlon.shape[0])
        latvar = nc.createVariable("lat", "f8", ("lat", "lon"))
        lonvar = nc.createVariable("lon", "f8", ("lat", "lon"))
        fracvar = nc.createVariable("frac", "f8", ("lat", "lon"))
        latvar[:] = self.dlat[...]
        latvar.units = 'deg N'
        lonvar[:] = self.dlon[...]
        lonvar.units = 'deg E'
        fracvar[:] = frac[...]
        fracvar.units = 'fraction'
        fracvar.description = (
            'fractional precipitation'
            'reduction for the International Drought Experiment (IDE).'
            'Calculated as ((1st percentile / 50th percentile)) of Qian et'
            'al (2006) 1948-2004 precipitation')
        nc.close()

    def show_reduction_pct(self, d, locations=None):
        """d: CLM_domain object
        locations: list of Location objects
        """
        frac = self.get_IDE_reduction()
        fig = plt.figure(figsize=(12, 6))
        ax1 = plt.subplot2grid((60, 11), (0, 0), colspan=5, rowspan=50)
        ax2 = plt.subplot2grid((60, 11), (0, 6), colspan=5, rowspan=50)
        ax3 = plt.subplot2grid((60, 11), (52, 0), colspan=11, rowspan=8)

        cmap, norm = colormap_nlevs.setup_colormap(0.0, 1.0, nlevs=11,
                                                   cmap=plt.get_cmap('YlGnBu'),
                                                   extend='neither')
        mworld = setup_worldmap(ax1)
        mcal = setup_calmap(ax2)
        cm = mworld.pcolormesh(d.get_lon(),
                               d.get_lat(),
                               ma.masked_invalid(frac),
                               cmap=cmap,
                               norm=norm,
                               latlon=True)
        cm = mcal.pcolormesh(d.get_lon(), d.get_lat(), ma.masked_invalid(frac),
                             cmap=cmap,
                             norm=norm,
                             latlon=True)
        if locations is not None:
            for here in locations:
                pt = mcal.scatter(here.lon[0], here.lat[0], latlon=True,
                                  marker='*', s=100, c='r')
                ax2.annotate(s="{:0.2f}".format(frac[here.clm_y, here.clm_x]),
                             xy=mcal(here.lon[0], here.lat[0]))

        cb = plt.colorbar(cm, cax=ax3,
                          orientation='horizontal')
        cb.ax.set_xlabel('1948 - 2005 pcp: (1st percentile / 50th percentile)')
        fig.tight_layout()
        fig.savefig(os.path.join(os.getenv('HOME'), 'plots', 'maptest',
                                 'IDE_pct_map.png'))
        plt.close(fig)


def get_f05g16_pcp():
    pcp_ncfile = os.path.join(os.getenv('SCRATCH'),
                              'qian_pcp_annual_totals.nc')
    qd = QianMonthlyPCPData(pcp_ncfile)
    qd.read_nc()
    d = CLM_Domain(fname=os.path.join('/', 'global',
                                      'cscratch1', 'sd',
                                      'twhilton',
                                      'archive',
                                      'CLM_f05_g16',
                                      'lnd',
                                      'hist',
                                      'CLM_f05_g16.clm2.h0.0050-01.nc'))
    qd.interpolate(d.get_lon(), d.get_lat())
    return qd


def setup_calmap(ax):
    """basic map of California with parallels, meridians, coastlines
    """
    m = Basemap(llcrnrlon=-125, llcrnrlat=30,
                urcrnrlon=-112, urcrnrlat=43,
                projection='mill',
                resolution='h',
                ax=ax)
    m.drawcoastlines(linewidth=1.25)
    m.fillcontinents(color='0.8', zorder=0)
    m.drawparallels(np.arange(30, 44, 4), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-125, -110, 5), labels=[0, 0, 0, 1])
    m.drawstates()
    m.drawcountries()
    return m


def setup_worldmap(ax):
    """basic map of world with parallels, meridians, coastlines
    """
    m = Basemap(llcrnrlon=-180, llcrnrlat=-80,
                urcrnrlon=180, urcrnrlat=80,
                projection='mill',
                ax=ax)
    m.drawcoastlines(linewidth=1.25)
    m.fillcontinents(color='0.8', zorder=0)
    m.drawparallels(np.arange(-80, 81, 20), labels=[0, 1, 0, 0])
    m.drawmeridians(np.arange(0, 360, 60), labels=[0, 0, 0, 1])
    return m


def check_results(qd, dlon, dlat):
    fig, ax = plt.subplots(2, 1, figsize=(8.5, 11))
    cmap = plt.get_cmap('Blues')
    cm = ax[0].pcolormesh(qd.lon, qd.lat, qd.pcp_all[0, ...],
                          cmap=cmap)
    for this_ax in ax:
        this_ax.set_xlabel('lon E')
        this_ax.set_ylabel('lat N')
    cb = plt.colorbar(cm, ax=ax[0])
    cm = ax[1].pcolormesh(dlon, dlat, qd.pcp[0, ...],
                          cmap=cmap)
    cb = plt.colorbar(cm, ax=ax[1])
    fig.savefig(os.path.join(os.getenv('HOME'), 'plots',
                             'qian_pcpinterp_test.png'))
    plt.close(fig)

    print "T62 min, max: ", qd.pcp_all[0, ...].min(), qd.pcp_all[0, ...].max()
    print "0.5 deg min, max: ", qd.pcp[0, ...].min(), qd.pcp[0, ...].max()


def site_summary(qd, site):

    y, x = (site.clm_y, site.clm_x)
    pcp = qd.pcp[:, y, x]
    pctl01, pctl50 = np.percentile(a=qd.pcp[:, y, x], q=(1, 50))
    fig, ax = plt.subplots(1, 1, figsize=[8, 8])
    ax.scatter(np.arange(len(pcp)) + 1948, pcp)
    ax.set_xlabel('year')
    ax.set_ylabel('pcp (mm)')
    plt.axhline(pctl01, color="#1b9e77", linestyle='--', linewidth=2.0,
                label='1st pctl')
    plt.axhline(pctl50, color="#d95f02", linestyle='--', linewidth=2.0,
                label='50th pctl')
    ax.set_title(site.name)
    ax.legend(loc='best')
    fig.savefig(os.path.join(os.getenv('HOME'), 'plots', 'maptest',
                             "{}_pcp.pdf".format(site.name.replace(' ', ''))))
    plt.close(fig)


if __name__ == "__main__":

    qd = get_f05g16_pcp()
    qd.pcp_reduction_to_netcdf('./pcp_reduction_frac.nc')
    # (domain_f05_g16, santacruz, mclaughlin,
    #  sierra_foothills, loma_ridge, sedgewick,
    #  boxsprings, ARM_SGP, harvard, wlef) = CLMf05g16_get_spatial_info()
    # qd.show_reduction_pct(domain_f05_g16,
    #                       (santacruz, mclaughlin,
    #                        sierra_foothills, loma_ridge,
    #                        sedgewick, boxsprings))
    # for this_site in (santacruz, mclaughlin, sierra_foothills,
    #                   loma_ridge, sedgewick, boxsprings, ARM_SGP):
    #     print "plotting summary: {}".format(this_site.name)
    #     site_summary(qd, this_site)
    # how to index:
    # qd.pcp[:, santacruz.clm_y, santacruz.clm_x]
    # qd.pcp[:, ARM_SGP.clm_y, ARM_SGP.clm_x]
