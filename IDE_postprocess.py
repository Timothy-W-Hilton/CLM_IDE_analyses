import os
import matplotlib.pyplot as plt
from spinup_diagnostics import CLM_spinup_analyzer as CSA
import IDE_locations


if __name__ == "__main__":
    d_archive = os.path.join(os.getenv('CSCRATCH'), 'archive')
    locs = IDE_locations.CLMf05g16_get_spatial_info()
    d = locs[0]
    locs = locs[1:]
    runs = {'ctl': CSA(os.path.join(d_archive, 'IDE_ctl',
                                    'lnd', 'hist'),
                       'ide_ctl'),
            'ide': CSA(os.path.join(d_archive, 'IDE_redpcp',
                                    'lnd', 'hist'),
                       'ide_redpcp')}

    s_per_tstep = 60 * 60 * 6  # tsteps are 6 hours
    vars = {'RAIN': {}, 'QVEGE': {}, 'QVEGT': {}}
    for run_name, run in runs.iteritems():
        for this_var in vars.keys():
            run.gather_filenames(glob_pat="*h2*.nc")
            vars[this_var][run_name] = run.parse_var(this_var,
                                                     lon_idx=locs[0].clm_x,
                                                     lat_idx=locs[0].clm_y)
            vars[this_var][run_name].data *= s_per_tstep
            vars[this_var][run_name].units = (
                vars[this_var][run_name].units.replace('/s', ''))

    fig, ax = plt.subplots(ncols=1, nrows=3, figsize=(7, 10))
    for this_var, this_ax in zip(vars.values(), ax):
        this_ax.plot(this_var['ctl'].data.squeeze(), 'b-', label='CTL')
        this_ax.plot(this_var['ide'].data.squeeze(), 'k--', label='IDE')
        this_ax.set_ylabel('{} ({})'.format(this_var['ctl'].varname,
                                            this_var['ctl'].units))
    ax[-1].set_xlabel('day of simulation')
    ax[0].set_title('{}'.format(locs[0].name))
    ax[0].legend(loc='best')
    fig.tight_layout()
    fig.savefig(os.path.join(os.getenv('CSCRATCH'),
                             '{}.pdf'.format(locs[0].name)))
    # plt.show()
