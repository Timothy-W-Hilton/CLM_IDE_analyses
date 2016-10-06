import matplotlib
# matplotlib.use('AGG')

import os
import matplotlib.animation as animation
import matplotlib.pyplot as plt
from datetime import datetime

from mpl_toolkits.basemap import Basemap

from timutils import midpt_norm
from timutils import colorbar_from_cmap_norm
import IDE_locations
import spinup_diagnostics as sd


def map_init(data):

    cmap, norm = get_norm_cmap(plt.get_cmap('Blues'))

    fig = plt.figure()
    ax = plt.subplot2grid((2, 12), (0, 0), rowspan=2, colspan=10)
    cbar_ax = plt.subplot2grid((2, 12), (0, 11), rowspan=2, colspan=1)

    clm_map = Basemap(ax=ax)
    clm_map.drawcoastlines()
    clm_map.drawstates()
    clm_map.drawcountries()

    locs = IDE_locations.CLMf05g16_get_spatial_info()
    domain = locs[0]
    locs = locs[1:]
    lat = domain.get_lat()
    lon = domain.get_lon()
    map_data = clm_map.pcolormesh(lon, lat,
                                  data[0, :, :].squeeze(),
                                  cmap=cmap,
                                  norm=norm,
                                  latlon=True)

    t_idx = 0
    colorbar_from_cmap_norm.colorbar_from_cmap_norm(
        cmap, norm, cbar_ax, None, data[t_idx, ...])
    cbar_ax.set_title('data')
    cbar_ax = None
    return (clm_map, fig, ax, cbar_ax, map_data)


def get_norm_cmap(cmap_arg=plt.get_cmap('Blues')):
    cmap, norm = midpt_norm.get_discrete_midpt_cmap_norm(vmin=4799,
                                                         vmax=4807,
                                                         midpoint=4803,
                                                         bands_above_mdpt=7,
                                                         bands_below_mdpt=7,
                                                         extend='both')
    return(cmap, norm)


def map_update(i, m, fig, ax, cbar_ax, t_idx, data, map_data):

    z_sfc = 0  # array index of surface
    map_data.set_array(data[i, :-1, :-1].squeeze().flatten())
    # ax[0].set_xlim([-160.0, -40.0])
    # ax[0].set_ylim([15.0, 85])

    t_str = '{}'.format(t_idx[i])
    ax.set_title(t_str)

    print('plotting for {}'.format(t_str))
    # fig.canvas.draw()


if __name__ == "__main__":
    outfile = 'clm.mp4'

    cases = ['IDE_ctl', 'IDE_redpcp']
    data_dirs = [os.path.join(os.getenv('CSCRATCH'), 'archive',
                              this_case, 'lnd', 'hist')
                 for this_case in cases]
    clm_runs = [sd.CLM_spinup_analyzer(data_dir=this_data,
                                       CASE=this_case)
                for this_data, this_case in zip(data_dirs, cases)]
    this_run = clm_runs[0]
    this_run.gather_filenames(glob_pat='*h1*')
    vardata = this_run.parse_var('WT')

    m, fig, ax, cbar_ax, map_data = map_init(vardata.data)

    im_ani = animation.FuncAnimation(fig,
                                     func=map_update,
                                     frames=10,   # len(cos['t']),
                                     interval=100,
                                     # blit=True,  # only update parts
                                     #             # that changed
                                     fargs=[m, fig, ax, cbar_ax, vardata.time,
                                            vardata.data, map_data])
    # note: not having ffmpeg available (`module load ffmpeg` at
    # NERSC) can create some undecipherable errors - check on that if
    # it's not working
    # im_ani.save(outfile, metadata={'artist': 'CLM'}, bitrate=100000)
    # plt.close(fig)
