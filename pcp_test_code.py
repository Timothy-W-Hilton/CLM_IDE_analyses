import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import os
import numpy as np
import numpy.ma as ma
from mpl_toolkits.basemap import Basemap
from RegionTester.region_tester import InUSState
import spinup_diagnostics
import qian_pcp_manipulator
from IDE_locations import CLMf05g16_get_spatial_info

def get_data():

    qd = qian_pcp_manipulator.get_f05g16_pcp()

    CLM_f05_g16 = spinup_diagnostics.CLM_spinup_analyzer(
        os.path.join('/', 'global',
                     'cscratch1', 'sd',
                     'twhilton',
                     'archive',
                     'CLM_f05_g16',
                     'lnd',
                     'hist'),
        'CLM_f05_g16')
    wt = CLM_f05_g16.parse_var('WT')
    return (qd, wt)

if __name__ == "__main__":
    (domain_f05_g16, santacruz, mclaughlin,
     sierra_foothills, loma_ridge, sedgewick,
     boxsprings, ARM_SGP) = CLMf05g16_get_spatial_info()
    qd, wt = get_data()

    lon = domain_f05_g16.lon
    lon[lon > 180] -= 360
    calmask_maker = qian_pcp_manipulator.CalMask(lon, domain_f05_g16.lat)
    calmask = calmask_maker.mask()

    tidx0 = 0
    tidx1 = 609
    pcparr = ma.masked_where(np.broadcast_to(calmask, [tidx1-tidx0, 384, 576]),
                             qd.pcp[tidx0:tidx1, ...])
    wtarr = ma.masked_where(np.broadcast_to(calmask, [tidx1-tidx0, 384, 576]),
                            wt.data[tidx0:tidx1, ...])
    fig, ax = plt.subplots()
    print "plotting"
    ax.scatter(pcparr.flatten(), wtarr.flatten())
    ax.set_title('spinup month 0 - 609')
    ax.set_xlabel('monthly pcp, mm')
    ax.set_ylabel('monthly water storage (mm)')
    print "saving"
    fig.savefig(os.path.join(os.getenv('HOME'), 'plots', 'pcp_scatter_mon0-609.png'))
    print "done saving"
    plt.close(fig)
