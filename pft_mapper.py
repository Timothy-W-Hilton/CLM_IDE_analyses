import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import netCDF4
from timutils import colormap_nlevs
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature


class IDEPaperMap(object):
    """setup a two-panel map with horizonal colorbar

    world map on left panel, California map on right panel
    """

    def __init__(self, vmin, vmax, ncolorlevs=10, colormap=None):
        self.vmin = vmin
        self.vmax = vmax
        if colormap is None:
            self.base_colormap = plt.get_cmap('Blues')
        else:
            self.base_colormap = colormap
        self.ncolorlevs = ncolorlevs
        self.mapworld = None
        self.mapcal = None
        self.axcbar = None
        self.fig = None
        self.projcal = ccrs.LambertConformal(central_longitude=262.5,
                                             central_latitude=38.5,
                                             standard_parallels=(38.5, 38.5))
        self.projglobe = ccrs.Miller()
        self.__setup_axes()
        self.__setup_maps()
        self.__get_cmap_norm()

    def __setup_axes(self):
        """create axes for two-panel plot with colorbar at bottom
        """
        self.fig = plt.figure(figsize=(12, 6))
        self.mapworld = plt.subplot2grid((60, 11), (0, 0),
                                         colspan=5, rowspan=50,
                                         projection=self.projglobe)
        self.mapcal = plt.subplot2grid((60, 11), (0, 6),
                                       colspan=5, rowspan=50,
                                       projection=self.projcal)
        self.axcbar = plt.subplot2grid((60, 11), (52, 0),
                                       colspan=11, rowspan=8)

    def __setup_maps(self):
        """basic maps with parallels, meridians, coastlines

        draw state boundaries in california panel
        """
        for this_ax in (self.mapcal, self.mapworld):
            this_ax.coastlines()
            this_ax.gridlines(draw_labels=False)
        self.mapcal.set_extent((-122, -114, 32, 43))
        states = NaturalEarthFeature(category='cultural',
                                     scale='50m',
                                     facecolor='none',
                                     name='admin_1_states_provinces_shp')
        self.mapcal.add_feature(states)

    def __get_cmap_norm(self):
        """get colormap, normalizer

        """
        self.cmap, self.norm = colormap_nlevs.setup_colormap(
            vmin=self.vmin,
            vmax=self.vmax,
            nlevs=self.ncolorlevs,
            cmap=self.base_colormap,
            extend='neither')

    def map_pft(self, pft_data):
        """plot grid cell percentage for one PFT on a map

        ARGS:
        pft_data: a PFTData object
        """
        for thismap in (self.mapworld, self.mapcal):
            cm = thismap.pcolormesh(pft_data.lon,
                                    pft_data.lat,
                                    pft_data.pft_data,
                                    transform=ccrs.PlateCarree(),
                                    cmap=self.cmap,
                                    norm=self.norm)
        plt.colorbar(cm, cax=self.axcbar, orientation='horizontal')
        self.axcbar.set_title('percentage')
        self.fig.suptitle('PFT: ' + pft_data.pft_name.replace('_', ' '))

    def draw_sites(self):
        """draw IDE site locations to map
        """
        sites = pd.read_csv(os.path.join('/Users', 'tim', 'work',
                                         'Code',
                                         'CLMHyperbolaFit',
                                         'IDE_sites.csv'))
        for s in sites.itertuples():
            self.mapcal.plot(s.lon, s.lat, marker='*',
                             markerfacecolor="None",
                             markersize=12,
                             transform=ccrs.Geodetic())


class PFTData(object):
    """reads CLM PFT percentages
    """

    def __init__(self, fname_pft_data, fname_pft_names, pft_idx):
        """PFTData constructor

        initializes fields fname_pft_data, fname_pft_names,
        pft_data, pft_names, pft_idx, lon, lat

        ARGS:
        fname_pft_data: full path to the netCDF file containing PFT
            percentages
        fname_pft_names: full path to the netCDF file containing PFT
            names
        pft_idx: index of the pft to parse (0 to 21, currently)
        """
        self.fname_pft_data = fname_pft_data
        self.fname_pft_names = fname_pft_names
        self.pft_idx = pft_idx

        self.__read_pft_names()
        self.__read()

    def __read_pft_names(self):
        """read pft names from netcdf file
        """
        nc = netCDF4.Dataset(self.fname_pft_names)
        pft_name = nc.variables['pftname'][self.pft_idx]
        nc.close()
        self.pft_name = ''.join(list(pft_name)).strip()

    def __read(self):
        """read a global map of PFT locations for one PFT
        """
        nc = netCDF4.Dataset(self.fname_pft_data)
        self.pft_data = nc.variables['PCT_PFT'][self.pft_idx, ...]
        self.lon = nc.variables['LONGXY'][...]
        self.lat = nc.variables['LATIXY'][...]
        nc.close()
        idx = self.lon > 180.0
        self.lon[idx] = -1.0 * ((-1.0 * self.lon[idx]) % 180.0)


if __name__ == "__main__":
    fname_pft_names = os.path.join('/Users', 'tim', 'work', 'Data',
                                   'CLM_Output',
                                   'pft-physiology.clm40.c130424.nc')
    fname_pft_data = '/Users/tim/work/Data/CLM_Output/CLM_PFTs.nc'

    map = IDEPaperMap(vmin=0.0, vmax=100.0, ncolorlevs=11)
    for this_pft_idx in np.arange(17):  # there are 17 PFTs
        print 'mapping PFT {:02d}'.format(this_pft_idx)
        try:
            this_pft = PFTData(fname_pft_data, fname_pft_names, this_pft_idx)
            map.map_pft(this_pft)
            map.draw_sites()
            map.fig.savefig('PFT{:02d}_pct_map.png'.format(this_pft_idx))
            plt.close(map.fig)
        except:
            print "error mapping PFT {:02d}".format(this_pft_idx)
            raise
