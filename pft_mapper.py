import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4
from timutils import midpt_norm, colorbar_from_cmap_norm
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature

class IDEPaperMap(object):
    """setup a two-panel map with horizonal colorbar

    world map on left panel, California map on right panel
    """

    def __init__(self):
        self.axworld = None
        self.axcal = None
        self.axcbar = None
        self.fig = None
        self.projcal = ccrs.LambertConformal(central_longitude=262.5,
                                             central_latitude=38.5,
                                             standard_parallels=(38.5, 38.5))
        self.projglobe = ccrs.Miller()

        self.__setup_axes()
        self.__setup_maps()

    def __setup_axes(self):
        """create axes for two-panel plot with colorbar at bottom
        """
        self.fig = plt.figure(figsize=(12, 6))
        self.axworld = plt.subplot2grid((60, 11), (0, 0),
                                        colspan=5, rowspan=50,
                                        projection=self.projglobe)
        self.axcal = plt.subplot2grid((60, 11), (0, 6),
                                      colspan=5, rowspan=50,
                                      projection=self.projcal)
        self.axcbar = plt.subplot2grid((60, 11), (52, 0),
                                       colspan=11, rowspan=8)

    def __setup_maps(self):
        """basic maps with parallels, meridians, coastlines

        draw state boundaries in california panel
        """
        for this_ax in (self.axcal, self.axworld):
            this_ax.coastlines()
            this_ax.gridlines(draw_labels=False)
        self.axcal.set_extent((-122, -114, 32, 43))
        states = NaturalEarthFeature(category='cultural',
                                     scale='50m',
                                     facecolor='none',
                                     name='admin_1_states_provinces_shp')
        self.axcal.add_feature(states)


class PFTReader(object):
    """reads CLM PFT percentages
    """

    def __init__(self, fname_pft_data, fname_pft_names):
        """PFTReader constructor

        initializes fields fname_pft_data, fname_pft_names, pft_data, pft_names

        ARGS:
        fname_pft_data: full path to the netCDF file containing PFT percentages
        """
        self.fname_pft_data = fname_pft_data
        self.fname_pft_names = fname_pft_names
        self.read_pft_names()
        self.pft_data = None

    def read_pft_names(self):
        """read pft names from netcdf file
        """
        nc = netCDF4.Dataset(self.fname_pft_names)
        pft_names = nc.variables['pftname'][:]
        nc.close()
        self.pft_names = [''.join(list(this_pft)).strip()
                          for this_pft in pft_names]

    def read(self, pft_idx):
        """read a global map of PFT locations for one PFT
        """
        nc = netCDF4.Dataset(self.fname_pft_data)
        self.pft_data = nc.variables['PCT_PFT'][pft_idx, ...]
        nc.close()


if __name__ == "__main__":
    fname_pft_names = os.path.join('/Users', 'tim', 'work', 'Data',
                                   'CLM_Output',
                                   'pft-physiology.clm40.c130424.nc')
    fname_pft_data = '/Users/tim/work/Data/CLM_Output/CLM_PFTs.nc'
    pfts = PFTReader(fname_pft_data, fname_pft_names)
