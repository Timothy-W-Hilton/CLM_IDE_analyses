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
from timutils import colormap_nlevs, colorbar_from_cmap_norm
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

        for t in np.arange(ntime):
            if np.mod(t, 100) == 0:
                print t, datetime.now()
            finterp = interpolate.RectBivariateSpline(
                self.lon[0, :],
                self.lat[:, 0],
                np.transpose(self.pcp_all[t, ...]))
            self.pcp[t, ...] = finterp.ev(dlon, dlat)
            self.pcp[self.pcp < 0] = 0.0


    def show_reduction_pct(self, d):
        """d: CLM_domain object
        """
        pct = (np.percentile(a=self.pcp, q=1, axis=0) /
               np.percentile(a=self.pcp, q=50, axis=0))
        fig = plt.figure(figsize=(8, 8))
        ax1 = plt.subplot2grid((6, 11), (0, 0), colspan=5, rowspan=5)
        ax2 = plt.subplot2grid((6, 11), (0, 6), colspan=5, rowspan=5)
        ax3 = plt.subplot2grid((6, 11), (5, 0), colspan=11, rowspan=1)

        cmap, norm = colormap_nlevs.setup_colormap(0.0, 1.0, nlevs=11,
                                                   cmap=plt.get_cmap('Blues'),
                                                   extend='neither')
        mworld = setup_worldmap(ax1)
        mcal = setup_calmap(ax2)
        cm = mworld.pcolormesh(d.get_lon(), d.get_lat(), ma.masked_invalid(pct),
                               cmap=cmap,
                               norm=norm,
                               latlon=True)
        cm = mcal.pcolormesh(d.get_lon(), d.get_lat(), ma.masked_invalid(pct),
                             cmap=cmap,
                             norm=norm,
                             latlon=True)
        cb = plt.colorbar(cm, cax=ax3,
                          orientation='horizontal')
        fig.tight_layout()
        fig.savefig(os.path.join(os.getenv('HOME'), 'plots', 'maptest',
                                 'IDE_pct_map.png'))
        plt.close(fig)


def setup_calmap(ax):
    """basic map of world with parallels, meridians, coastlines
    """
    m = Basemap(llcrnrlon=-125, llcrnrlat=32,
                urcrnrlon=-114, urcrnrlat=43,
                projection='mill',
                resolution='h',
                ax=ax)
    m.drawcoastlines(linewidth=1.25)
    m.fillcontinents(color='0.8', zorder=0)
    m.drawparallels(np.arange(30, 44, 2), labels=[1,1,0,0])
    m.drawmeridians(np.arange(-114, -125, 2), labels=[0,0,0,1])
    m.drawstates()
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
    m.drawparallels(np.arange(-80,81,20),labels=[1,1,0,0])
    m.drawmeridians(np.arange(0,360,60),labels=[0,0,0,1])
    return m

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

    print "T62 min, max: ", qmd.pcp_all[0, ...].min(), qmd.pcp_all[0, ...].max()
    print "0.5 deg min, max: ", qmd.pcp[0, ...].min(), qmd.pcp[0, ...].max()




def site_summary(qmd, site):

    pcp = qmd.pcp[:, site.clm_y, site.clm_x]
    fig, ax = plt.subplots(1, 1, figsize=[8, 8])
    ax.scatter(np.arange(len(pcp)) + 1948, pcp)
    ax.set_xlabel('year')
    ax.set_ylabel('pcp (mm)')
    fig.savefig(os.path.join(os.getenv('HOME'), 'plots',
                             "{}_pcp.pdf".format(site.name.replace(' ', ''))))
    plt.close(fig)


if __name__ == "__main__":
    pcp_ncfile = os.path.join(os.getenv('SCRATCH'),
                              'qian_pcp_annual_totals.nc')
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

    (domain_f05_g16, santacruz, mclaughlin,
     sierra_foothills, loma_ridge, ARM_SGP) = CLMf05g16_get_spatial_info()
    qmd.show_reduction_pct(domain_f05_g16)

    #how to index:
    # qmd.pcp[:, santacruz.clm_y, santacruz.clm_x]
    # qmd.pcp[:, ARM_SGP.clm_y, ARM_SGP.clm_x]
