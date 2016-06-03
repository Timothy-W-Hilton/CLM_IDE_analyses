"""a set of tools useful for parsing and returning descriptions of the
STEM domain (grid cell longitude and latitude coordinates, grid cell
vertical coordinates), and also for some domain-oriented calculations
(is a point inside or outside of the domain, etc.)
"""

import os
import os.path
import matplotlib
import numpy as np
import warnings
import netCDF4
import stem_pytools.STEM_parsers as sp
from mpl_toolkits import basemap
from scipy.spatial import cKDTree


class Location(object):
    """class to contain the location of a point.  Contains user-accessible
    fields lon and lat.
    """
    def __init__(self, lon, lat, name):
        self.lon = lon
        self.lat = lat
        self.name = name

    def get_clm_xy(self, clm_domain):
        """find the CLM grid cellccontaining location

        place the nearest CLM grid cell x and y coordinates in
        fields clm_x, clm_y
        ARGS:
        clm_domain (CLM_Domain object): description of the domain
        """
        self.clm_y, self.clm_x = clm_domain.find_nearest_xy(
            self.lon,
            self.lat)

class CLM_Domain(object):
    """class to contain CLM domain attributes.
    """
    def __init__(self, fname):
        """Populates fields: fname, lon, lat, topo, area
        ARGS:
        fname: CLM output file from which to read lon, lat
        """
        self.fname = fname
        nc = netCDF4.Dataset(fname)
        self.lon1D = nc.variables['lon'][:]
        self.lat1D = nc.variables['lat'][:]
        self.topo =  nc.variables['topo'][:]
        self.area = nc.variables['area'][:]
        self.lon, self.lat = np.meshgrid(self.lon1D, self.lat1D)
        return(None)

    def get_lat(self):
        """return the CLM latitude grid"""
        return(self.STEM_lat)

    def get_lon(self):
        """return the CLM longitude grid"""
        return(self.STEM_lon)

    def get_topo(self):
        """return the CLM topo grid"""
        return(self.topo)

    def find_nearest_xy(self, lon, lat):
        """find indices of the nearest grid point for a set of lon, lat points.

        Given a set of arbitrary (lon, lat) positions, find the horizontal
        (x, y) CLM grid indices of the nearest CLM grid cell to each
        position.

        PARAMETERS
        ----------
        lon, lat: ndarray; of arbitrary longitudes and latitudes.  Must
           contain the same number of elements.

        RETURNS:
        an N-element tuple of X and Y indices, one index per observation.
            The indices are the closest point on the CLM grid to each
            point in (lon, lat).  N is therefore equal to the number of
            dimensions in lon and lat.

        """
        # convert spherical lon, lat coordinates to cartesian coords. Note
        # that these x,y,z are 3-dimensional cartesian coordinates of
        # positions on a sphere, and are thus different from the x,y,z
        # *indices* of the STEM grid.
        x, y, z = lon_lat_to_cartesian(lon, lat)
        xs, ys, zs = lon_lat_to_cartesian(self.lon, self.lat)

        # use a K-dimensional tree to find the nearest neighbor to (x,y,z)
        # from the points within (xs, ys, zs).  A KD tree is a data
        # structure that allows efficient queries of K-dimensional space (K
        # here is 3).
        tree = cKDTree(np.dstack((xs.flatten(),
                                  ys.flatten(),
                                  zs.flatten())).squeeze())

        d, inds = tree.query(
            np.dstack((x, y, z)).squeeze(), k=1)

        return(np.unravel_index(inds, self.lon.shape))



def lon_lat_to_cartesian(lon, lat, R=1):
    """
    Convert spherical coordinates to three-dimensional Cartesian
    coordinates.

    calculates three dimensional cartesian coordinates (x, y, z) for
    specified longitude, latitude coordinates on a sphere with radius
    R.  Written and posted at earthpy.org by Oleksandr Huziy.
    http://earthpy.org/interpolation_between_grids_with_ckdtree.html
    accessed 19 Mar 2014 by Timothy W. Hilton.

    PARAMETERS
    ==========
    lon; np.ndarray: longitude values
    lat; np.ndarray: latitude values
    R: scalar; radius of the sphere.

    RETURNS
    =======
    three element tuple containing X, Y, and Z, one element per
    lon,lat pair.
    """
    lon_r = np.radians(lon)
    lat_r = np.radians(lat)

    x = R * np.cos(lat_r) * np.cos(lon_r)
    y = R * np.cos(lat_r) * np.sin(lon_r)
    z = R * np.sin(lat_r)
    return x, y, z


def demo():
    """Demonstrates using the domain module to find the indices of a
    couple of key locations on a 1 degree by 1 degree Northern
    Hemisphere grid that are inside the area specifed by a STEM topo
    file.

    """
    domain = CLM_Domain(fname=os.path.join('/', 'global',
                                           'cscratch1', 'sd',
                                           'twhilton',
                                           'archive',
                                           'CLM_f05_g16',
                                           'lnd',
                                           'hist',
                                           'CLM_f05_g16.clm2.h0.0050-01.nc'))
    santacruz = Location((122.03089741, ), (36.9741, ))
    idx = domain.find_nearest_xy(np.array((santacruz.lon)),
                                 np.array((santacruz.lat)))
    return idx
    # mymap = domain.get_mapobj()
    # mymap.drawcoastlines()
    # mymap.plot(domain.bnd_lon, domain.bnd_lat, latlon=True)
    # mymap.scatter(sib_lon,
    #               sib_lat,
    #               latlon=True,
    #               c='gray',
    #               facecolors='none',
    #               marker='o')
    # mymap.scatter(sib_lon[in_STEM_domain],
    #               sib_lat[in_STEM_domain],
    #               latlon=True,
    #               c='red',
    #               marker='x')
