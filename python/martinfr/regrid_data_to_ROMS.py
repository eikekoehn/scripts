#!/usr/bin/env python

'''
Regridding observational data to ROMS grid using numpy griddata

Author: Martin Frischknecht
Date:   11.03.2014
'''

#########################################################
# Load requirements
#########################################################

# Import Python modules
import numpy as np
import netCDF4
from scipy.interpolate import griddata

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


# ROMS setup
setup = 'pactcs30'

# Input File
varin = 'nitrate' # 'nitrate', 'oxygen' 'phosphate' or 'salinity'
varout = 'NO3'
infile = '/net/kryo/work/updata/CARS/{}_cars2009.nc'.format(varin)

# Output File
outfile = '/net/kryo/work/martinfr/Data/Obs/CARS/{}_cars2009_{}.nc'.format(varin,setup)


#########################################################
# Regrid OBS to ROMS
#########################################################

# OBS file
ncf = netCDF4.Dataset(infile, 'r')
# Read latlon (might be different every time)
lonOBS = ncf.variables['lon'][:]
latOBS = ncf.variables['lat'][:]
lonOBS, latOBS = np.meshgrid(lonOBS,latOBS)

# Obs mask
maskOBS = ncf.variables['mean'][0].mask
# Mask LatLon on land
lonOBS = np.ma.masked_array(lonOBS,maskOBS)
latOBS = np.ma.masked_array(latOBS,maskOBS)
# Obs data
data = ncf.variables['mean'][0]

# Prepare ROMS latlon
lonROMS = np.ma.masked_where(romsGrd.mask_rho==0,romsGrd.lon_rho)
latROMS = np.ma.masked_where(romsGrd.mask_rho==0,romsGrd.lat_rho)
    
# Regrid surface values (this is only one time step here!! Annual mean values.)
# If you need to regrid monthly values, or even daily fields, use the fast interpolation routine
tmp = griddata((lonOBS.ravel(),latOBS.ravel()), data.ravel(), (lonROMS,latROMS),fill_value=1e33)
# Mask land, fill values, extrapolation
tmp = np.ma.masked_where(romsGrd.mask_rho==0,tmp)
tmp = np.ma.masked_greater(tmp,data.max())
tmp = np.ma.masked_less(tmp,data.min())

# Open new NetCDF file to write to
ncfile = netCDF4.Dataset(outfile, 'w', format='NETCDF4')


#########################################################
# Write Output
#########################################################

# Create dimensions
ncfile.createDimension('time', 0)
ncfile.createDimension('xi_rho', romsGrd.SX)
ncfile.createDimension('eta_rho', romsGrd.SY)

# Create NetCDF Variable
varDims = ('eta_rho','xi_rho')
outdata = ncfile.createVariable(varout,'f4',varDims,fill_value=1e33)

# Write to NetCDF
outdata[:] = tmp

# Close the file and end the program.
ncfile.close()
