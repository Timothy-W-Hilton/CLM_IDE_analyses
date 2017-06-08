import sys
import os
import numpy as np
from numpy import ma
from PIL import Image  # to overlay site labels, crop
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from timutils import midpt_norm


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


class WorldCalMap(object):
    """plot a map of the world and a map of California with one colorbar
    """
    def __init__(self, figsize=(12, 6), interp=False):
        """class constructor
        """
        matplotlib.rcParams.update({'font.size': 14})
        self.figsize = figsize
        self.interp = interp

        self.fig = plt.figure(figsize=(12, 6))
        self.ax1 = plt.subplot2grid((60, 110), (0, 0), colspan=50, rowspan=50)
        self.ax2 = plt.subplot2grid((60, 110), (0, 53), colspan=50, rowspan=50)
        self.ax3 = plt.subplot2grid((60, 110), (52, 4),
                                    colspan=90, rowspan=10)

        self.mworld = setup_worldmap(self.ax1)
        self.mcal = setup_calmap(self.ax2)

        self.mworld.plot(x=(self.mcal.llcrnrlon, self.mcal.llcrnrlon,
                            self.mcal.urcrnrlon, self.mcal.urcrnrlon,
                            self.mcal.llcrnrlon),
                         y=(self.mcal.llcrnrlat, self.mcal.urcrnrlat,
                            self.mcal.urcrnrlat, self.mcal.llcrnrlat,
                            self.mcal.llcrnrlat),
                         latlon=True,
                         color="#d95f02",
                         linewidth=2.0)

    def plot(self, data, lon, lat,
             cmap_arg=plt.get_cmap('YlGnBu'),
             midpoint=0.0,
             bands_above=5,
             bands_below=5,
             vmin=0.0, vmax=1.0,
             locations=None,
             cbar_tstr=None,
             extend='both',
             site_labels=None):
        """ locations: list of Location objects
        """

        cmap, norm = midpt_norm.get_discrete_midpt_cmap_norm(
            vmin, vmax,
            midpoint=midpoint,
            bands_above_mdpt=bands_above,
            bands_below_mdpt=bands_below,
            this_cmap=cmap_arg,
            extend=extend)
        for this_map in (self.mworld, self.mcal):
            cm = this_map.pcolormesh(lon,
                                     lat,
                                     ma.masked_invalid(data),
                                     cmap=cmap,
                                     norm=norm,
                                     latlon=True)
        # draw california map boundary box on world map

        if locations is not None:
            try:
                for here in locations:
                    pt = self.mcal.scatter(here.lon[0], here.lat[0],
                                           latlon=True,
                                           marker='*', s=100, c='r')
                    if site_labels is "fraction":
                        self.ax2.annotate(s="{:0.2f}".format(data[here.clm_y,
                                                                  here.clm_x]),
                                          xy=self.mcal(here.lon[0],
                                                       here.lat[0]))
                    elif site_labels is "names":
                        self.ax2.annotate(s=here.name,
                                          xy=self.mcal(here.lon[0],
                                                       here.lat[0]))
            except IndexError:
                print('{sitename}: clm_x or clm_y'
                      ' exceeds domain bounds'.format(
                          sitename=here.name))
            except:
                print("Unexpected error:", sys.exc_info()[0])
                raise
        cb = plt.colorbar(cm, cax=self.ax3, orientation='horizontal')
        if cbar_tstr is not None:
            cb.ax.set_xlabel(cbar_tstr)
        self.fig.tight_layout()

        # self.fig.savefig(
        #     os.path.join(os.getenv('HOME'), 'plots', 'maptest',
        #                  'IDE_pct_map_interp{}.png'.format(
        #                      self.lat.size != self.dlat.size)))
        # plt.close(self.fig)

    def label_crop_save(self, fname_image, fname_labels=None, dpi=(300, 300)):
        """overlay site name labels, crop out whitespace

        fname_labels (str): full path to png image containing labels
            (default is ./site_labels.png'
        """

        if fname_labels is None:
            fname_labels = './site_labels.png'

        self.fig.savefig('tmp.png')

        map_image = Image.open("tmp.png")
        labels_image = Image.open(fname_labels)
        final_image = Image.new("RGBA", map_image.size)
        # 1080, 576, 11, 14
        # left upper right lower
        final_image = Image.alpha_composite(final_image, map_image)
        final_image = Image.alpha_composite(final_image, labels_image)
        final_image = final_image.crop((0, 13, 1075, 587))
        final_image.save(fname_image, dpi=dpi)
        map_image.close()
        labels_image.close()
        os.remove("tmp.png")
