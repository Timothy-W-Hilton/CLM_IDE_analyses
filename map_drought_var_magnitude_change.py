"""draw a map of (control simultation to drought simulation) change in
(1) magnitude, (2) percent magnitude, and (3) day-of-year timing for a
variable's maximum (minimum)
"""

import os
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from world_cal_maps import WorldCalMap
from IDE_locations import CLMf05g16_get_spatial_info


class DOYVar(object):
    """read a netCDF file of gridded mean annual cycle data
    """
    def __init__(self, fname, varname):
        self.fname = fname
        self.varname = varname
        self.lon = None
        self.lat = None

    def parse(self):
        nc = netCDF4.Dataset(self.fname)
        self.data = nc.variables[self.varname][...]
        # read longitude *v*ectors
        self.vlon = nc.variables['lon'][...]
        self.vlat = nc.variables['lat'][...]
        self.lon, self.lat = np.meshgrid(self.vlon, self.vlat)
        nc.close()

    def get_annual_min(self):
        return np.amin(self.data, axis=0)

    def get_annual_max(self):
        return np.amax(self.data, axis=0)

    def get_annual_min_DOY(self):
        return np.argmin(self.data, axis=0)

    def get_annual_max_DOY(self):
        return np.argmax(self.data, axis=0)


if __name__ == "__main__":

    sp_info = CLMf05g16_get_spatial_info()
    domain = sp_info[0]
    locs = sp_info[1:]

    fnames = {}
    varobjs = {}
    runs = ['ctl', 'redpcp']
    for this_run in runs:
        fnames[this_run] = os.path.join(
            '/', 'global', 'cscratch1', 'sd', 'twhilton',
            'daily_CLM_output', 'output',
            'IDE_{}_daily_FPSN_doymean.nc'.format(this_run))
        varobjs[this_run] = DOYVar(fnames[this_run], 'FPSN')
        varobjs[this_run].parse()

    d_mag = (varobjs['redpcp'].get_annual_max() -
             varobjs['ctl'].get_annual_max())
    d_mag_pct = ((varobjs['redpcp'].get_annual_max() -
                 varobjs['ctl'].get_annual_max()) /
                 varobjs['ctl'].get_annual_max()) * 100.0
    d_doy = (varobjs['redpcp'].get_annual_max_DOY() -
             varobjs['ctl'].get_annual_max_DOY())

    wcm_mag = WorldCalMap()
    wcm_mag.plot(data=d_mag,
                 lon=varobjs['ctl'].lon,
                 lat=varobjs['ctl'].lat,
                 vmin=d_mag.min(), vmax=0.2,
                 midpoint=0.0,
                 bands_above=3,
                 bands_below=8,
                 cmap_arg=plt.get_cmap('PuOr'),
                 cbar_tstr=(r'$\Delta$GPP annual max (control to drought)'
                            ' [$\mu$mol m$^{-2}$ s$^{-1}$]'),
                 locations=locs,
                 site_labels="names")
    wcm_mag.crop_save('FPSN_Delta_ann_max.png')

    wcm_mag_pct = WorldCalMap()
    wcm_mag_pct.plot(data=d_mag_pct,
                     lon=varobjs['ctl'].lon,
                     lat=varobjs['ctl'].lat,
                     vmin=-100.0, vmax=25.0,
                     midpoint=0.0,
                     bands_above=3,
                     bands_below=10,
                     cmap_arg=plt.get_cmap('PuOr'),
                     cbar_tstr=(r'$\Delta$GPP annual max'
                                ' (controll to drought) [percent]'),
                     locations=locs,
                     site_labels="names")
    wcm_mag_pct.crop_save('FPSN_PctDelta_ann_max.png')

    wcm_doy = WorldCalMap()
    wcm_doy.plot(data=d_doy,
                 lon=varobjs['ctl'].lon,
                 lat=varobjs['ctl'].lat,
                 vmin=-50, vmax=50,
                 bands_above=7,
                 bands_below=7,
                 cmap_arg=plt.get_cmap('PuOr'),
                 cbar_tstr=(r'change in timing of GPP annual max'
                            '(control to drought) [days]'),
                 locations=locs,
                 site_labels="names")
    wcm_doy.crop_save('FPSN_DeltaDOY_ann_max.png')
