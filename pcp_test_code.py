import matplotlib
# matplotlib.use('AGG')
import matplotlib.pyplot as plt
import os
import numpy as np
import numpy.ma as ma
from mpl_toolkits.basemap import Basemap
from RegionTester.region_tester import InUSState
import spinup_diagnostics
import qian_pcp_manipulator
from IDE_locations import CLMf05g16_get_spatial_info

def get_wt_ann_sum(wt):
    yr_list = np.split(wt.data, np.arange(12, 609, 12))
    yr_sums = np.stack(map(lambda x: np.sum(x, axis=0), yr_list))
    return yr_sums


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
     boxsprings, ARM_SGP, harvard, wlef) = CLMf05g16_get_spatial_info()
    qd, wt = get_data()
    wt_sum = get_wt_ann_sum(wt)
    lon = domain_f05_g16.lon
    lon[lon > 180] -= 360
    calmask_maker = qian_pcp_manipulator.CalMask(lon, domain_f05_g16.lat)
    calmask = calmask_maker.mask()

    tidx0 = 0
    tidx1 = 50
    pcparr = ma.masked_where(np.broadcast_to(calmask, [tidx1-tidx0, 384, 576]),
                             qd.pcp[tidx0:tidx1, ...])
    wtarr = ma.masked_where(np.broadcast_to(calmask, [tidx1-tidx0, 384, 576]),
                            wt_sum[tidx0:tidx1, ...])
    fig, ax = plt.subplots()
    print "plotting"
    ax.scatter(pcparr.flatten(), wtarr.flatten())
    ax.set_title('California cells, spinup year {} - {}'.format(tidx0, tidx1))
    ax.set_xlabel('annual pcp, mm')
    ax.set_ylabel('annual water storage (mm)')
    print "saving"
    fig.savefig(os.path.join(os.getenv('HOME'), 'plots',
                             'pcp_scatter_yr{}_{}.png'.format(tidx0, tidx1)))
    print "done saving"
    plt.close(fig)

    tidx0 = 0
    tidx1 = 50
    loc = sedgewick
    pcparr = qd.pcp[tidx0:tidx1, loc.clm_y, loc.clm_x]
    wtarr =  wt_sum[tidx0:tidx1, loc.clm_y, loc.clm_x]
    fig, ax = plt.subplots()
    print "plotting"
    ax.scatter(pcparr.flatten(), wtarr.flatten())
    ax.set_title('{}, spinup year {} - {}'.format(loc.name, tidx0, tidx1))
    ax.set_xlabel('annual pcp, mm')
    ax.set_ylabel('annual water storage (mm)')
    print "saving"
    fig.savefig(os.path.join(os.getenv('HOME'), 'plots',
                             'pcp_scatter_{}_yr{}_{}.png'.format(
                                 loc.name.replace(' ', ''), tidx0, tidx1)))
    print "done saving"
    plt.close(fig)

    loc = santacruz
    for loc in (santacruz, mclaughlin,
                sierra_foothills, loma_ridge, sedgewick,
                boxsprings, ARM_SGP, harvard, wlef):
        pcparr = qd.pcp[tidx0:tidx1, loc.clm_y, loc.clm_x]
        wtarr =  wt_sum[tidx0:tidx1, loc.clm_y, loc.clm_x]
        fig, ax = plt.subplots()
        ax2 = ax.twinx()
        pcp_lines = ax.plot(pcparr, 'b-', marker='x', label='pcp')
        wt_lines = ax2.plot(wtarr, 'k--', marker='+', label='WT')
        ax.set_ylim(bottom=100, top=1600)
        ax2.set_ylim(bottom=46000, top=69000)
        ax.set_xlim((-5, 50))
        ax2.set_xlim((-5, 50))
        ax.set_xlabel('spinup year')
        ax.set_ylabel('pcp (mm)')
        ax2.set_ylabel('water storage (mm)')
        ax.legend(loc='upper left')
        ax2.legend(loc='upper right')
        ax.set_title(loc.name)
        fig.savefig(os.path.join(os.getenv('HOME'), 'plots',
                                 'pcp_WT_timeseries_{}.pdf'.format(
                                     loc.name.replace(' ', ''))))
        plt.close(fig)
