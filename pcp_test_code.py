import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import os
import sys
import numpy as np
import netCDF4
from datetime import datetime
from mpl_toolkits.basemap import Basemap
from RegionTester.region_tester import InUSState

if __name__ == "__main__":
    qian_pcp_file = os.path.join(os.getenv('SCRATCH'),
                                 'qian_pcp_monthly_totals.nc')

    nc = netCDF4.Dataset(qian_pcp_file)
    lat = nc.variables['LATIXY'][...]
    lon = nc.variables['LONGXY'][...]

    lon2 = lon.copy()
    idx = lon2 > 180
    lon2[idx] = lon2[idx] - 360.0

    California = InUSState('')
    California.get_state_shape('California')

    iscal = np.empty(lon2.shape, dtype=bool)
    iscal[...] = False
    print datetime.now()
    for (x, y), this_lon in np.ndenumerate(lon2):
        iscal[x, y] = California.point_inside(lon2[x, y], lat[x, y])
    print datetime.now()

    # usamask = np.reshape(isusa, lat.shape)
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
    m.scatter(lon2[calmask], lat[calmask], latlon=True, marker='x', color='r', s=2)
    fig.savefig(os.path.join(os.getenv('HOME'), 'plots', 'maptest', 'usamap.png'))
    plt.close(fig)
