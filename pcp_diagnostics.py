from IDE_locations import CLMf05g16_get_spatial_info
from spinup_diagnostics import CLM_spinup_analyzer
from qian_pcp_manipulator import QianMonthlyPCPData

if __name__ == "__main__":
    # get sites
    (domain_f05_g16, santacruz, mclaughlin,
     sierra_foothills, loma_ridge, sedgewick,
     boxsprings, ARM_SGP) = CLMf05g16_get_spatial_info()

    # get pcp
    pcp_ncfile = os.path.join(os.getenv('SCRATCH'),
                              'qian_pcp_annual_totals.nc')
    qd = QianMonthlyPCPData(pcp_ncfile)
    qd.read_nc()
    d = CLM_Domain(fname=os.path.join('/', 'global',
                                      'cscratch1', 'sd',
                                      'twhilton',
                                      'archive',
                                      'CLM_f05_g16',
                                      'lnd',
                                      'hist',
                                      'CLM_f05_g16.clm2.h0.0050-01.nc'))
    qd.interpolate(d.get_lon(), d.get_lat())

    # get total water storage
    CLM_f05_g16 = CLM_spinup_analyzer(os.path.join('/', 'global',
                                                   'cscratch1', 'sd',
                                                   'twhilton',
                                                   'archive',
                                                   'CLM_f05_g16',
                                                   'lnd',
                                                   'hist'),
                                      'CLM_f05_g16')
    wt = CLM_f05_g16.parse_var('WT', verbose=True)
