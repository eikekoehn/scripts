#!/usr/bin/env python

'''
Regridding observational data to ROMS grid using fast triangulatin routine.

Author: Martin Frischknecht
Date:   Jan 2018
'''

#########################################################
# Load requirements
#########################################################

# Import Python modules
import numpy as np
import netCDF4

# Import my modules
from myromstools import getGrid

# ROMS setup
setup = 'pactcs30'

# Set grid/clm file
grdfile = '/net/kryo/work/martinfr/Roms/Inputs/{}/{}_grd.nc'.format(setup,setup)
# Read grid
romsGrd = getGrid(grdfile)
# Add lat lon
romsGrd.getLatLon()

# Input File (data to be regridded)
infile = '/net/kryo/work/updata/POCexport/deVries2017/SIMPLE_TRIM_output.nc'
varin = 'NPP'
varout = 'NPP'

# Output File
outfile = '/net/kryo/work/martinfr/Data/Obs/POCexport/deVries2017_NPP_interp2{}.nc'.format(setup)

#########################################################
# Helper functions
#########################################################

# This is the fast interpolation routine based on some triangulation...
# It is superfast, and I think for all of our purposes just perfectly well suited.

import scipy.spatial.qhull as qhull

def interp_tri(xy):
    tri = qhull.Delaunay(xy)
    return tri

def interpolate(values, tri,uv,d=2):
    simplex = tri.find_simplex(uv)
    vertices = np.take(tri.simplices, simplex, axis=0)
    temp = np.take(tri.transform, simplex, axis=0)
    delta = uv- temp[:, d]
    bary = np.einsum('njk,nk->nj', temp[:, :d, :], delta)
    return np.einsum('nj,nj->n', np.take(values, vertices),  np.hstack((bary, 1.0 - bary.sum(axis=1, keepdims=True))))


#########################################################
# Regrid data and write output
#########################################################

# Read latlon of observational data (might be different every time!)
ncf = netCDF4.Dataset(infile, 'r')
lonOBS, latOBS = ncf.variables['LON'][0],ncf.variables['LAT'][0]

# Obs Lon Lat
lonOBS_flat = lonOBS.flatten()
latOBS_flat = latOBS.flatten()
xy = np.zeros([len(latOBS_flat),2])
xy[:,0] = lonOBS_flat
xy[:,1] = latOBS_flat

# ROMS Lon Lat
lonROMS_flat = romsGrd.lon_rho.flatten()
latROMS_flat = romsGrd.lat_rho.flatten()
xyROMS = np.zeros([len(latROMS_flat),2])
xyROMS[:,0] = lonROMS_flat
xyROMS[:,1] = latROMS_flat

# Compute triangulation of grid points
tri = interp_tri(xy)

# Third dimension (mostly time, except for deVries NPP where it is actually different model versions)
tidx = range(len(ncf.variables[varin]))

    
#########################################################
# Write Output
#########################################################

# Open new NetCDF file to write to
ncfile = netCDF4.Dataset(outfile, 'w', format='NETCDF4')

# Create dimensions
ncfile.createDimension('time', 0)
ncfile.createDimension('xi_rho', romsGrd.SX)
ncfile.createDimension('eta_rho', romsGrd.SY)

# Create NetCDF Variable
varDims = ('time','eta_rho','xi_rho')
outdata = ncfile.createVariable(varout,'f4',varDims,fill_value=romsGrd.FillValue)

# Loop through third dimension
for t in tidx:
    print('Interpolating {}'.format(t))
    data_tmp = ncf.variables[varin][t]

    # Interpolate horizontally
    tmp = interpolate(data_tmp.ravel(),tri,xyROMS)
    tmp = tmp.reshape(romsGrd.SY,romsGrd.SX)
    # Mask invalid data, extrapolation and land points
    tmp = np.ma.masked_invalid(tmp)
    tmp = np.ma.masked_greater(tmp,data_tmp.max())
    tmp = np.ma.masked_less(tmp,data_tmp.min())
    tmp = np.ma.masked_where(romsGrd.mask_rho==0,tmp)

    # Write to NetCDF
    outdata[t] = tmp

# Add attribtue (if necessary or available...)
outdata.description = ncf.variables[varin].description

# Close the file and end the program.
ncfile.close()
