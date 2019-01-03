# -*- coding: utf-8 -*-
"""
Created on Tue Nov  6 15:45:23 2018

@author: martinfr
"""

#########################################################
# Load requirements
#########################################################

# Import Python modules
import numpy as np
import netCDF4

# Import my modules
from myromstools import getGrid
from myplotlib import easyplot

# ROMS setup
setup = 'pactcs30'

# Set grid/clm file
grdfile = '/net/kryo/work/martinfr/Roms/Inputs/{}/{}_grd.nc'.format(setup,setup)
# Read grid
romsGrd = getGrid(grdfile)
# Add lat lon
romsGrd.getLatLon()

# Load mask for the CalCS
f = '/net/kryo/work/martinfr/Data/Masks/pactcs30/pactcs30_masks_C_CalCS_only_BudgetAnalysis_final.nc'
ncfmasks = netCDF4.Dataset(f, 'r')
maskCalCS = ncfmasks.variables['mask_CalCS'][:]
easyplot(maskCalCS)

# Load dist2coast file
f = '/net/kryo/work/martinfr/Data/pactcs30/pactcs30_dist2coast.nc'
ncfmasks = netCDF4.Dataset(f, 'r')
dcoast = ncfmasks.variables['dcoast'][:]

#########################################################
# Load some data
#########################################################

# ROMS setup
setup = 'pactcs30'
runtag = 'newsrc_hc05_monthly_pactcs30'

# Infile
f = '/net/kryo/work/martinfr/Roms/Output/{}/{}/avg/z_pactcs30_avg_clim.nc'.format(setup,runtag)
ncf = netCDF4.Dataset(f,'r')
# Load temperature data and compute annual mean
data = np.ma.mean(ncf.variables['temp'][:],axis=0)
NZclim = data.shape[0]

#########################################################
# Possibility 1: Create offshore mean section across the 
# California Current System (central coast according to 
# the mask maskCalCS)
#########################################################

# Offshore distance array
dist_offshore = np.arange(0,105,5).tolist()
dist_offshore.extend(np.arange(110,510,10))
dist_offshore.extend(np.arange(520,1020,20))
dist_offshore.extend(np.arange(1025,1525,25))
dist_offshore = np.asarray(dist_offshore)
# Define center of offshore distance bins for plots and axis labelling
dist_offshore_center = dist_offshore[:-1]+(dist_offshore[1:]-dist_offshore[:-1])/2.

# Apply mask CalCS
dcoast_masked = np.ma.masked_where(maskCalCS==0.0,dcoast)
# Apply mask CalCS
dcoast_masked = np.ma.masked_where(maskCalCS==0.0,dcoast)

# Compute mean cross shore section
data_offshore =  np.ones([NZclim,len(dist_offshore)-1])*np.nan
print('Computing offshore section...')
# Loop through offshore distance array and average selected data
for d in range(1,len(dist_offshore)):
    # Select data within a distance bin
    tmp_i, tmp_j =np.ma.where(np.logical_and(dcoast_masked>=dist_offshore[d-1],dcoast_masked<=dist_offshore[d]))
    tmp_data =  data[:,tmp_i,tmp_j]
    # Compute masked array mean
    data_offshore[:,d-1] = np.ma.mean(tmp_data,axis=1)
# Mask invalid data (and np.nan from creating the array)
data_offshore = np.ma.masked_invalid(data_offshore)

# Rotate matrix and plot it
easyplot(np.fliplr(np.flipud(data_offshore)))


#########################################################
# Possibility 2: Regrid to a certain line...
# Same as regridding (either use triangulation routine or
# griddata)
#########################################################

## USING TRIANGULATION ROUTINE

import scipy.spatial.qhull as qhull

def interp_tri(xy):
    tri = qhull.Delaunay(xy)
    return tri

def interpolate(values, tri,uv,d=2):
    simplex = tri.find_simplex(uv)
    vertices = np.take(tri.simplices, simplex, axis=0)
    temp = np.take(tri.transform, simplex, axis=0)
    delta = uv-temp[:, d]
    bary = np.einsum('njk,nk->nj', temp[:, :d, :], delta)
    return np.einsum('nj,nj->n', np.take(values, vertices),  np.hstack((bary, 1.0 - bary.sum(axis=1, keepdims=True))))

# Line coordinates (e.g. at latitude 35N)
LatLine = np.ones([50])*35.
LonLine = np.linspace(220,240,50)
LonLat_Line = np.zeros([len(LatLine),2])
LonLat_Line[:,0] = LonLine
LonLat_Line[:,1] = LatLine

# Obs Lon Lat
lonROMS_flat = romsGrd.lon_rho.flatten()
latROMS_flat = romsGrd.lat_rho.flatten()
xy = np.zeros([len(latROMS_flat),2])
xy[:,0] = lonROMS_flat
xy[:,1] = latROMS_flat

# Compute triangulation of grid points
tri = interp_tri(xy)

# Interpolate
data_offshore = np.ones([NZclim,len(LatLine)])*np.nan
for d in range(NZclim):
    data_offshore[d] = interpolate(data[d].ravel(),tri,LonLat_Line,2)

# Rotate matrix and plot it
easyplot(np.flipud(data_offshore))



## USING GRIDDATA

from scipy.interpolate import griddata
# Regrid to line (super slow with griddata)
data_offshore = np.ones([NZclim,len(LatLine)])*np.nan
for d in range(NZclim):
    print d
    data_offshore[d] = griddata((romsGrd.lon_rho.ravel(),romsGrd.lat_rho.ravel()), data[d].ravel(), (LonLine,LatLine),fill_value=1e33)

# Rotate matrix and plot it
easyplot(np.flipud(data_offshore))


# YOU WILL SEE THAT BOTH RESULTS ARE BASICALLY THE SAME. JUST ANOTHER ARGUMENT TO ALWAYS GO WITH THE FAST TRIANGULATION METHOD !!!
