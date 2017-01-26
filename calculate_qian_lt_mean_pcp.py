"""Calculate 1948-2004 annual mean Qian et al (2006) precipitation

calculates from $CSCRATCH/qian_pcp_annual_pcp_totals.nc

qian_pcp_annual_pcp_totals.nc is produced by qian_totaller.QianTotaller

REFERENCES
Qian, T., A. Dai, K. E. Trenberth, and K. W. Oleson (2006)
Simulation of Global Land Surface Conditions from 1948 to 2004.
Part I: Forcing Data and Evaluations, Journal of Hydrometeorology, 7(5), 953'
"""

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import netCDF4
from qian_pcp_manipulator import get_f05g16_pcp
from timutils.io import delete_if_exists


def qian_ltmean_pcp_ncdf(outfile, qd):
    delete_if_exists(outfile)
    nc = netCDF4.Dataset(outfile, mode='w', format="NETCDF3_CLASSIC")
    nc.createDimension("lat", qd.dlat.shape[0])
    nc.createDimension("lon", qd.dlon.shape[1])
    latvar = nc.createVariable("lat", "f8", ("lat", "lon"))
    lonvar = nc.createVariable("lon", "f8", ("lat", "lon"))
    pcpvar = nc.createVariable("pcp", "f8", ("lat", "lon"))
    latvar[:] = qd.dlat[...]
    latvar.units = 'deg N'
    lonvar[:] = qd.dlon[...]
    lonvar.units = 'deg E'
    n_secs_per_year = 60 * 60 * 24 * 365  # seconds in a year
    pcpvar[:] = qd.data.mean(axis=0) * n_secs_per_year
    pcpvar.units = 'mm yr-1'
    pcpvar.description = (
        'Qian et al (2006) 1948-2004 long-term mean precipitation.\n'
        'Qian, T., A. Dai, K. E. Trenberth, and K. W. Oleson (2006),\n'
        'Simulation of Global Land Surface Conditions from 1948 to 2004.\n'
        'Part I: Forcing Data and Evaluations,\n'
        'Journal of Hydrometeorology, 7(5), 953')
    nc.close()

if __name__ == "__main__":
    qd = get_f05g16_pcp(interp_flag=False)
    qdi = get_f05g16_pcp(interp_flag=True)
    qian_ltmean_pcp_ncdf('test.nc', qdi)

    # --- plotting code below ---
    vmax = np.max((qdi.data.max(), qd.data.max()))
    plt.imshow(qd.data[0, ...],
               cmap='Blues',
               vmax=vmax)
    plt.colorbar()
    plt.figure()
    plt.imshow(qdi.data[0, ...], cmap='Blues', vmax=vmax)
    plt.colorbar()

    plt.figure()
    m = Basemap()
    m.pcolormesh(qdi.dlon, qdi.dlat,
                 qdi.data.mean(axis=0),
                 cmap='Blues',
                 vmax=vmax,
                 latlon=True)
    m.drawcoastlines()
    plt.colorbar()

    qds = qd.data.sum(axis=(1, 2))
    qdis = qdi.data.sum(axis=(1, 2))
    plt.figure()
    # plt.pcolormesh(np.stack((qds, qdis)).transpose(), cmap='RdBu')
    plt.pcolormesh((qds / qdis)[:, np.newaxis], cmap='RdBu')
    plt.colorbar()
