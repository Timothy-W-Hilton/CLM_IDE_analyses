import os
import matplotlib.pyplot as plt

def site_pcp_scatter(qd, wt, tidx0, tidx1, site):
    """pcp vs water storage scatter plot for one site

    Parameters
    ----------
    qd : qian_pcp_manipulater.QianMonthlyPCPData object
    	contains the monthly total precipitation time series in field pcp
    wt : CLM_Var object
	contains the water storage time series in field data

    """
    fig, ax = plt.subplots()
    print "plotting"
    pcp_sc = qd.pcp[tidx0:tidx1, site.clm_y, site.clm_x].squeeze()
    wt_sc = wt.data[tidx0:tidx1, site.clm_y, site.clm_x].squeeze()
    ax.scatter(pcp_sc, wt_sc)
    ax.set_xlabel('monthly pcp (mm)')
    ax.set_ylabel('monthly water storage (mm)')
    ax.set_title('{} spinup month {}-{}'.format(site.name, tidx0, tidx1))
    fname = os.path.join(os.getenv('HOME'), 'plots',
                         'pcp_scatter_{}_mon{}_{}.png'.format(
                             site.name.replace(' ', ''),
                             tidx0, tidx1))
    print "saving {}".format(fname)
    fig.savefig(fname)
    print "done saving"
    plt.close(fig)
