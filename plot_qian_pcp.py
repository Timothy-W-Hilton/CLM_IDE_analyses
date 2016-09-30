import os
import matplotlib.pyplot as plt
import itertools

import qian_pcp_manipulator as qpm
from timutils import colormap_nlevs

if __name__ == "__main__":
    qian_dir = os.path.join('/project', 'projectdirs', 'ccsm1',
                            'inputdata', 'atm', 'datm7',
                            'atm_forcing.datm7.Qian.T62.c080727',
                            'Precip6Hrly')
    years = [1972, 1972]
    months = [1, 2]
    fnames = [os.path.join(qian_dir,
                           'clmforc.Qian.c2006.T62.Prec.{}-{:02d}.nc'.format(
                               y, m))
              for y in years for m in months]

    pcp = [qpm.QianMonthlyPCPData(f, 'PRECTmms') for f in fnames]
    for this_pcp in pcp:
        print "parsing {}".format(os.path.basename(this_pcp.fname))
        this_pcp.read_nc()

    fig, ax = plt.subplots(nrows=len(years), ncols=1)
    for y, m, this_pcp, this_ax in zip(itertools.cycle(years), months,
                                       pcp, list(ax)):
        this_m = qpm.setup_worldmap(this_ax)
        cmap, norm = colormap_nlevs.setup_colormap(0.0, 0.0053803, nlevs=11,
                                                   cmap=plt.get_cmap('Blues'),
                                                   extend='neither')
        this_m.pcolormesh(this_pcp.lon,
                          this_pcp.lat,
                          this_pcp.data_all[0, ...],
                          latlon=True)
        this_ax.annotate("{} {} pcp".format(y, m),
                         xy=(0, 1),
                         xycoords='axes fraction')
    fig.savefig(os.path.join(os.getenv('CSCRATCH'), 'qian_pcp.png'))
