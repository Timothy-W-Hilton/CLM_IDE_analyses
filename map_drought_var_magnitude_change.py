"""draw a map of control--drought change in magnitude in a variable's max (min)
"""
import os
import netCDF4
import numpy as np
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
                 varobjs['ctl'].get_annual_max())
    d_doy = (varobjs['redpcp'].get_annual_max_DOY() -
             varobjs['ctl'].get_annual_max_DOY())

    wcm_mag = WorldCalMap()
    wcm_mag.plot(data=d_mag,
                 lon=varobjs['ctl'].lon,
                 lat=varobjs['ctl'].lat,
                 vmin=d_mag.min(), vmax=d_mag.max(),
                 cbar_tstr=r'$\Delta$ annual max',
                 locations=locs)

    wcm_mag_pct = WorldCalMap()
    wcm_mag_pct.plot(data=d_mag_pct,
                     lon=varobjs['ctl'].lon,
                     lat=varobjs['ctl'].lat,
                     vmin=-2.0, vmax=2.0,
                     cbar_tstr=r'$\% \Delta$ annual max',
                     locations=locs)

    wcm_doy = WorldCalMap()
    wcm_doy.plot(data=d_doy,
                 lon=varobjs['ctl'].lon,
                 lat=varobjs['ctl'].lat,
                 vmin=-20, vmax=20,
                 cbar_tstr=r'$\Delta$ DOY',
                 locations=locs)
