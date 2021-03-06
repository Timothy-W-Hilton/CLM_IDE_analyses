import os
import numpy as np
import pandas as pd
from clm_domain import CLM_Domain, Location


def CLMf05g16_get_spatial_info():
    domain_f05_g16 = CLM_Domain(os.path.join(os.getenv('CSCRATCH'),
                                             'archive',
                                             'CLM_f05_g16',
                                             'lnd',
                                             'hist',
                                             'CLM_f05_g16.clm2.h0.0051-09.nc'))
    santacruz = Location((-122.03089741, ), (36.9741, ), 'Younger Lagoon')
    santacruz.get_clm_xy(domain_f05_g16)
    mclaughlin = Location((-122.431667, ), (38.873889, ), 'McLaughlin')
    mclaughlin.get_clm_xy(domain_f05_g16)
    sierra_foothills = Location((-121.311623, ), (39.249407, ),
                                'Sierra Foothill')
    sierra_foothills.get_clm_xy(domain_f05_g16)
    loma_ridge = Location((-117.7, ), (33.74, ),
                          'Loma Ridge')
    loma_ridge.get_clm_xy(domain_f05_g16)
    ARM_SGP = Location((-97.4888, ), (36.6058, ),
                       'ARM Southern Great Plains')
    ARM_SGP.get_clm_xy(domain_f05_g16)
    # coordinates from https://en.wikipedia.org/wiki/Sedgwick_Reserve
    sedgewick = Location((-120.016667, ), (34.7, ), 'Sedgwick')
    sedgewick.get_clm_xy(domain_f05_g16)
    # Box Springs coordinates from http://www.gps-coordinates.net
    # using the description on the NRS website: "Riverside County,
    # adjacent to Riverside campus on Box Springs Mountain; 4 km (2.5
    # mi.) east of the city of Riverside"
    # (http://www.ucnrs.org/reserves/box-springs-reserve.html)
    boxsprings = Location((-117.303321, ), (33.98472, ), 'Box Springs')
    boxsprings.get_clm_xy(domain_f05_g16)
    #http://www.findlatitudeandlongitude.com/?loc=+mammoth+lakes%2C+california#.WBN35neZMo9
    mammoth_lakes = Location((-118.972079, ), (37.648546, ), 'SNARL')
    mammoth_lakes.get_clm_xy(domain_f05_g16)
    # coordinates from http://atmos.seas.harvard.edu/lab/hf/hfsite.html
    harvard = Location((-72.171478, ), (42.537755, ), 'Harvard Forest')
    harvard.get_clm_xy(domain_f05_g16)
    # coordinates from http://fluxnet.ornl.gov/site/1036
    wlef = Location((-90.2723, ), (45.9459, ), 'WLEF')
    wlef.get_clm_xy(domain_f05_g16)
    carrizo_plain = Location((-119.8633, ), (35.1899, ), 'Carrizo Plain')
    carrizo_plain.get_clm_xy(domain_f05_g16)
    return (domain_f05_g16, santacruz, mclaughlin,
            sierra_foothills, loma_ridge, sedgewick, boxsprings, ARM_SGP,
            harvard, wlef, mammoth_lakes, carrizo_plain)


def sites_to_csv(fname='IDE_sites.csv'):
    """create a CSV file of locations specified in CLMf05g16_get_spatial_info
    """
    sites = CLMf05g16_get_spatial_info()
    sites = sites[1:]  # first element is the domain
    sites_df = pd.DataFrame(
        {'name': np.array([s.name for s in sites]),
         'lat': np.array([s.lat[0] for s in sites]),
         'lon': np.array([s.lon[0] for s in sites]),
         'clm_x': np.array([s.clm_x for s in sites]),
         'clm_y': np.array([s.clm_y for s in sites])})
    sites_df.to_csv(fname)
