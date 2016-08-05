import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import os
import numpy as np
import numpy.ma as ma
import spinup_diagnostics
import qian_pcp_manipulator
from IDE_locations import CLMf05g16_get_spatial_info


def get_LE_ann_sum(LE):
    yr_list = np.split(LE.data, np.arange(12, 609, 12))
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
    qsoil = CLM_f05_g16.parse_var('QSOIL')
    qvege = CLM_f05_g16.parse_var('QVEGE')
    qvegt = CLM_f05_g16.parse_var('QVEGT')

    le = spinup_diagnostics.CLM_var(
        varname='LE',
        varname_long='latent heat flux',
        data=spinup_diagnostics.mm_s_2_w_m2_s(qsoil.data +
                                              qvege.data +
                                              qvegt.data),
        time=qsoil.time,
        units='W/m^2',
        missing_value=qsoil.missing_value)
    return (qd, wt, le)


def draw_pcp_scatter(pcp, var, mask, tidx0, tidx1=50,
                     xlim=(-100, 2000),
                     ylim=(-500, 3000)):
    """draw a scatter plot of precipitation vs. a water-driven land variable
    """
    pcparr = ma.masked_where(np.broadcast_to(mask, [tidx1-tidx0, 384, 576]),
                             qd.pcp[tidx0:tidx1, ...])
    varsum = get_LE_ann_sum(var)
    arr = ma.masked_where(np.broadcast_to(mask, [tidx1-tidx0, 384, 576]),
                          varsum[tidx0:tidx1, ...])
    fig, ax = plt.subplots()
    print "plotting"
    ax.scatter(pcparr.flatten(), arr.flatten())
    ax.set_title('California cells, spinup year {} - {}'.format(tidx0, tidx1))
    ax.set_xlabel('annual pcp, mm')
    ax.set_ylabel('annual total {} ({})'.format(var.varname, var.units))
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    print "saving"
    fig.savefig(os.path.join(os.getenv('HOME'), 'plots', 'new_pcp',
                             'pcp_{}_scatter_yr{}_{}.png'.format(var.varname,
                                                                 tidx0,
                                                                 tidx1)))
    print "done saving"
    plt.close(fig)


def draw_pcp_scatter_loc(pcp, var, loc, tidx0, tidx1=50,
                         xlim=(-100, 2000),
                         ylim=(-500, 3000)):
    """location-specific PCP--LE scatter
    """
    pcparr = qd.pcp[tidx0:tidx1, loc.clm_y, loc.clm_x]
    varsum = get_LE_ann_sum(var)
    arr = varsum[tidx0:tidx1, loc.clm_y, loc.clm_x]
    fig, ax = plt.subplots()
    print "plotting"
    ax.scatter(pcparr.flatten(), arr.flatten())
    ax.set_title('{}, spinup year {} - {}'.format(loc.name, tidx0, tidx1))
    ax.set_xlabel('annual pcp, mm')
    ax.set_ylabel('annual total {} ({})'.format(var.varname, var.units))
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    fname = os.path.join(os.getenv('HOME'), 'plots', 'new_pcp',
                         'pcp_{}_scatter_{}_yr{}_{}.pdf'.format(
                             var.varname,
                             loc.name.replace(' ', ''), tidx0, tidx1))
    print "saving {}".format(fname)
    fig.savefig(fname)
    print "done saving"
    plt.close(fig)


def draw_pcp_timeseries(pcp, var, loc, tidx0, tidx1=50,
                        ylim=(0, 850)):
    """location-specific PCP--LE timeseries
    """
    pcparr = qd.pcp[tidx0:tidx1, loc.clm_y, loc.clm_x]
    varsum = get_LE_ann_sum(var)
    arr = varsum[tidx0:tidx1, loc.clm_y, loc.clm_x]
    fig, ax = plt.subplots()
    ax2 = ax.twinx()
    ax.plot(pcparr, 'b-', marker='x', label='pcp')
    ax2.plot(arr, 'k--', marker='+', label='LE')
    ax.set_xlim((-5, 50))
    ax2.set_xlim((-5, 50))
    ax.set_ylim(bottom=0, top=1600)
    ax2.set_ylim(ylim)
    ax.set_xlabel('spinup year')
    ax.set_ylabel('pcp (mm)')
    ax2.set_ylabel('annual {} ({})'.format(var.varname, var.units))
    ax.legend(loc='upper left')
    ax2.legend(loc='upper right')
    ax.set_title(loc.name)
    fig.savefig(os.path.join(os.getenv('HOME'), 'plots', 'new_pcp',
                             'pcp_{}_timeseries_{}.pdf'.format(
                                 var.varname,
                                 loc.name.replace(' ', ''))))
    plt.close(fig)


def get_range(loclist, vals):
    """get the range of values in CLM variable var across all locations in loclist
    """
    x = [loc.clm_x for loc in loclist]
    y = [loc.clm_y for loc in loclist]
    loc_vals = ma.masked_greater(vals[:, y, x], 1e20)
    return (loc_vals.min(), loc_vals.max())


def get_region_range(mask, vals):
    x, y = np.where(np.logical_not(mask))
    region_vals = ma.masked_greater(vals[:, x, y], 1e20)
    return(region_vals.min() * 0.95, region_vals.max() * 1.05)

if __name__ == "__main__":
    (domain_f05_g16, santacruz, mclaughlin,
     sierra_foothills, loma_ridge, sedgewick,
     boxsprings, ARM_SGP, harvard, wlef) = CLMf05g16_get_spatial_info()
    qd, wt, LE = get_data()
    LE_sum = get_LE_ann_sum(LE)
    lon = domain_f05_g16.lon
    lon[lon > 180] -= 360
    calmask_maker = qian_pcp_manipulator.CalMask(lon, domain_f05_g16.lat)
    calmask = calmask_maker.mask()

    locations = (santacruz, mclaughlin,
                 sierra_foothills, loma_ridge, sedgewick,
                 boxsprings, ARM_SGP, harvard, wlef)
    for this_var in (wt, LE):
        this_sum = get_LE_ann_sum(this_var.data)
        ylim = get_range(locations, this_sum)
        ylim_cal = get_region_range(calmask, this_sum)
        for tidx0 in (0, 40, 45):
            draw_pcp_scatter(qd, this_var, calmask, tidx0, ylim=ylim_cal)
            for loc in locations:
                draw_pcp_scatter_loc(qd, this_var, loc, tidx0, ylim=ylim)
        for loc in locations:
            draw_pcp_timeseries(qd, this_var, loc,
                                tidx0=0, tidx1=50, ylim=ylim)
