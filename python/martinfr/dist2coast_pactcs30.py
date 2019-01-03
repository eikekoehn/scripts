#!/usr/bin/env python
'''
Computes distance to the coastline and saves to netCDF file.

Martin Frischknecht, Mar 2015
'''

#########################################################
# Load requirements
#########################################################

# Import modules
import numpy as np
import netCDF4
# My stuff
from myfunc import points2dist
from myromstools import getGrid,writeNetCDF
from myplotlib import easyplot

# Load ROMS grid data
grdfile =  '/net/kryo/work/martinfr/Roms/Inputs/pactcs30/pactcs30_grd.nc'
romsGrd = getGrid(grdfile)
# Load LatLon
romsGrd.getLatLon()

#########################################################
# Find and correct coastline
#########################################################

# Coast points in xi and eta direction
is_coastx= (romsGrd.mask_rho[:-1,:]+romsGrd.mask_rho[1:,:]==1)
is_coasty= (romsGrd.mask_rho[:,:-1]+romsGrd.mask_rho[:,1:]==1)

# Control plots
#easyplot(is_coastx)
#easyplot(is_coasty)

# Make MANUAL corrections to islands etc (if needed... use easyplot)

# Correct Santa Barbara Channel Islands
is_coastx[470:489,354:395] = 0
is_coastx[410:420,410:420] = 0
is_coastx[470:473,414:415] = 0
is_coasty[470:489,354:395] = 0
is_coasty[410:420,410:420] = 0
is_coasty[471:473,413:415] = 0

# lon/lat of grid interfaces at coast
lonu=0.5*(romsGrd.lon_rho[:-1,:]+romsGrd.lon_rho[1:,:])
latu=0.5*(romsGrd.lat_rho[:-1,:]+romsGrd.lat_rho[1:,:])
lonv=0.5*(romsGrd.lon_rho[:,:-1]+romsGrd.lon_rho[:,1:])
latv=0.5*(romsGrd.lat_rho[:,:-1]+romsGrd.lat_rho[:,1:])
lon_cx=lonu[is_coastx]
lat_cx=latu[is_coastx]
lon_cy=lonv[is_coasty]
lat_cy=latv[is_coasty]

# All coastal lon/lat
lonc=np.hstack((lon_cx[:],lon_cy[:]))
latc=np.hstack((lat_cx[:],lat_cy[:]))
ncoast=len(lonc)

# Compute closest distance to coast (this takes some time, but you only have to do it ONCE, hopefully)
dcoast=np.ma.ones([romsGrd.SY,romsGrd.SX])*1e33
for i in xrange(romsGrd.SY):
    print('Eta: {} out of {}').format(i,romsGrd.SY)
    for j in xrange(romsGrd.SX):
        tmp = []
        for k in xrange(ncoast):
            tmp.append(points2dist((romsGrd.lat_rho[i,j],romsGrd.lon_rho[i,j]),(latc[k],lonc[k])))
        # Save closest distance to array
        tmp = np.array(tmp)
        dcoast[i,j] = np.min(tmp)
        
dcoast = np.ma.masked_greater(dcoast,1e30)
dcoast = dcoast/1000. # in km
dcoast = np.ma.masked_where(romsGrd.mask_rho==0,dcoast)


#########################################################
# Write NetCDF output
#########################################################

# Output file
outfile = '/net/kryo/work/martinfr/Data/pactcs30/pactcs30_dist2coast_TESTTESTTEST.nc'

# Create dictionary and provide it to the writing function
dData = {'dcoast': dcoast}
# Create dictionary (for each variable that you want to save, specify dimension names and size)
dDims = {'dcoast':{'dims':['eta','xi'],'dimlens':[romsGrd.SY,romsGrd.SX]}}

# Write dictionary and all the variables that it contains to netCDF file
writeNetCDF(romsGrd,dData,outfile,dDims)

