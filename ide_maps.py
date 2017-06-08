"""demonstrate plotting gridded data on map of N America
"""
import numpy as np
import matplotlib.pyplot as plt
import netCDF4
from timutils import midpt_norm
from world_cal_maps import WorldCalMap


def parse_inflection_points(fname):
    """read CLM hyperbola fit inflection points from netCDF
    """
    nc = netCDF4.Dataset(fname)
    infl_rats = nc.variables['inflc_rat'][:]
    lon1d = nc.variables['lon'][:]
    lat1d = nc.variables['lat'][:]
    lat, lon = np.meshgrid(lat1d, lon1d)
    nc.close()
    return infl_rats, lon, lat


def define_colormap():
    my_cmap, my_norm = midpt_norm.get_discrete_midpt_cmap_norm(
        vmin=0.0,  # infl_rats.min(),
        vmax=2.0,  # infl_rats.max(),
        midpoint=1.0,
        bands_above_mdpt=6,
        bands_below_mdpt=6,
        this_cmap=plt.get_cmap(plt.get_cmap('BrBG')),  # cool
        extend='both')
    # get rid of the white band in the middle by replacing it and the
    # bands for all larger values with the colors returned from
    # running get_discrete_midpt_cmap_norm with bands_above_mdpt=7.
    # TODO: This is a hack and needs to be fixed.
    my_cmap.colors[5:] = [np.array([0.80868898,  0.92441369,
                                    0.90788159,  1.0]),
                          np.array([0.59477127,  0.84183007,
                                    0.80392158,  1.0]),
                          np.array([0.34625146,  0.69181086,
                                    0.65305654,  1.0]),
                          np.array([0.13587082,  0.52433681,
                                    0.49296426,  1.0]),
                          np.array([0.00322953,  0.37093426,
                                    0.33679355,  1.0])]
    # make the [0.8, 1.0] band red
    my_cmap.colors[4] = np.array((1.0, 0.0, 0.0))
    return my_cmap, my_norm


def fix_longitudes(lons):
    """Convert longitudes in [0, 360] to [-180, 180]
    """
    idx = lons > 180.0
    lon[idx] = -1.0 * ((-1.0 * lons[idx]) % 180.0)
    return lon


if __name__ == "__main__":
    infl_rats, lon, lat = parse_inflection_points('hyperbola_fits.nc')
    print "parsed inflection points"
    lon = fix_longitudes(lon)
    print "longitudes fixed"
    wcm = WorldCalMap()
    print "initialized maps"
    my_cmap, my_norm = define_colormap()
    print "defined colormap"
    wcm.plot(infl_rats, lon, lat,
             vmin=0.0,
             vmax=2.0,
             midpoint=1.0,
             bands_above=6,
             bands_below=6,
             # \u2013 is unicode en dash
             cbar_tstr=(u'(annual pcp at GPP\u2013pcp slope change) / '
                        u'(1948\u20132004 annual mean pcp)'),
             cmap_arg=my_cmap)
    print "drew maps"
    wcm.label_crop_save("PCP_GPP_ratio_map.png")
    print "overlaid labels, etc."
