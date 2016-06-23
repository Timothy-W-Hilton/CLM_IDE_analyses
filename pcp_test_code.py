import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import os
import numpy as np
from geopy.geocoders import Nominatim
import netCDF4
from datetime import datetime
from mpl_toolkits.basemap import Basemap

if __name__ == "__main__":
    nc = netCDF4.Dataset(os.path.join(os.getenv('SCRATCH'),
                                     'qian_pcp_monthly_totals.nc'))
    lat = nc.variables['LATIXY'][...]
    lon = nc.variables['LONGXY'][...]


    lon2 = lon.copy()
    idx = lon2 > 180
    lon2[idx] = lon2[idx] - 360.0
    # result = rg.search(('37.38605', '-122.08385'))
    # for this_lat, this_lon in zip(lat.flatten()[:5], lon.flatten()[:5]):
    #     print this_lat, this_lon

    results = rg.search(zip(lat.flatten(), lon2.flatten()))

    isusa = np.array([this_cell['cc'] == 'US' for this_cell in results])
    iscal = np.array([this_cell['admin1'] == 'California' for this_cell in results])

    usamask = np.reshape(isusa, lat.shape)
    calmask = np.reshape(iscal, lat.shape)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,4.5))
    m = Basemap(llcrnrlon=-180, llcrnrlat=-80,
                urcrnrlon=180, urcrnrlat=80,
                projection='mill',
                ax=ax)
    m.drawcoastlines(linewidth=1.25)
    m.fillcontinents(color='0.8', zorder=0)
    m.drawparallels(np.arange(-80,81,20),labels=[1,1,0,0])
    m.drawmeridians(np.arange(0,360,60),labels=[0,0,0,1])
    m.scatter(lon2, lat, latlon=True, marker='x', color='k', s=2)
    m.scatter(lon2[usamask], lat[usamask], latlon=True, marker='x', color='r', s=2)
    fig.savefig(os.path.join(os.getenv('HOME'), 'plots', 'maptest', 'usamap.png'))
    plt.close(fig)
