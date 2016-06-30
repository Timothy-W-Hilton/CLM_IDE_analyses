import os
from clm_domain import CLM_Domain, Location

def CLMf05g16_get_spatial_info():
    domain_f05_g16 = CLM_Domain(os.path.join(os.getenv('SCRATCH'),
                                             'archive',
                                             'CLM_f05_g16',
                                             'lnd',
                                             'hist',
                                             'CLM_f05_g16.clm2.h0.0051-09.nc'))
    santacruz = Location((-122.03089741, ), (36.9741, ), 'Santa Cruz')
    santacruz.get_clm_xy(domain_f05_g16)
    mclaughlin = Location((-122.431667, ), (38.873889, ), 'McLaughlin NRS')
    mclaughlin.get_clm_xy(domain_f05_g16)
    sierra_foothills = Location((-121.311623, ), (39.249407, ),
                                'Sierra Foothill Research Extension Center')
    sierra_foothills.get_clm_xy(domain_f05_g16)
    loma_ridge = Location((-117.7, ), (33.74, ),
                          'Loma Ridge Global Change Experiment')
    loma_ridge.get_clm_xy(domain_f05_g16)
    ARM_SGP = Location((-97.4888, ), (36.6058, ),
                       'ARM Southern Great Plains')
    ARM_SGP.get_clm_xy(domain_f05_g16)
    # coordinates from https://en.wikipedia.org/wiki/Sedgwick_Reserve
    sedgewick = Location((-120.016667, ), (34.7, ), 'Sedgewick NRS')
    sedgewick.get_clm_xy(domain_f05_g16)
    # Box Springs coordinates from http://www.gps-coordinates.net
    # using the description on the NRS website: "Riverside County,
    # adjacent to Riverside campus on Box Springs Mountain; 4 km (2.5
    # mi.) east of the city of Riverside"
    # (http://www.ucnrs.org/reserves/box-springs-reserve.html)
    boxsprings = Location((-117.303321, ), (33.98472, ), 'Box Springs')
    boxsprings.get_clm_xy(domain_f05_g16)
    return (domain_f05_g16, santacruz, mclaughlin,
            sierra_foothills, loma_ridge, sedgewick, boxsprings, ARM_SGP)
