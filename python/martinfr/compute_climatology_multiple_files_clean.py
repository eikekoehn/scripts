#!/usr/bin/python

# Computes a monthly climatology from multiple ROMS output files.
# Martin Frischknecht, Mar 2015 (modified from Damian Loher)

import netCDF4
import numpy as np
import os
import sys

# Where is the data stored?
host = 'kryo'
# ROMS setup
setup = 'pactcs30'
# Run
run_name = 'newsrc_hc05_monthly_{}'.format(setup)


# Input files (example for hindcast simulation from 1979 to 2016):
dir_in = '/net/{}/work/martinfr/Roms/Output/{}/{}/avg'.format(host,setup,run_name)
filelist = ['{}_{}_avg.nc'.format(setup,yr) for yr in range(1979,2017)]

# Output directory
dir_out = dir_in

# Add path to files
filelist = [dir_in+'/'+f for f in filelist]

# Output file (clm):
clm_file = '{}/{}_avg_clim.nc'.format(dir_out,setup)

# Find out number of time records:
trec=0
for f in filelist:
    nc_in = netCDF4.Dataset(f,'r')
    trec += len(nc_in.dimensions['time'])

# Read attributes from first file in list
nc_in = netCDF4.Dataset(filelist[0],'r')
navg = nc_in.getncattr('navg')
dt = float(nc_in.getncattr('dt'))
print('Number of time records used: {}'.format(trec))
print('dt = {}'.format(dt))
print('navg = {}'.format(navg))

# Number of days between two consecutive time records:
days_per_trec = float(dt*navg)/86400
# Number of days per month:
days_per_month = 365./12. # cycle_length/12.
# Number or time records per month and year:
trec_per_month = 1.
trec_per_year = 12*trec_per_month
print('Number of time records per month: {}'.format(trec_per_month))
# How many months/years are covered by the avg file:
months = trec / trec_per_month
years = trec / (12*trec_per_month)
print('Months covered by avg file: {}'.format(months))
print('Years covered by avg file: {}'.format(years))
if years == 0:
    # Time period is too short:
    print('Need at least one year. Exiting')
    sys.exit(1)
if years-np.floor(years) != 0.0:
    print('WARNING: final year seems to be incomplete')
    print('Will use {} time records out of {}.'.format(np.floor(years)*12*trec_per_month,trec))
    # Adapt years:
    years = np.floor(years)
print('Number of complete years: {}'.format(years))
print('Dividing fields by: {}'.format(years*trec_per_month))

# Number of time records used for computing climatology:
nr_trec = years*12*trec_per_month
nr_trec_per_file = nr_trec/len(filelist)
# Number of time records in clim (output) file:
nr_trec_clim = 12
print('number of time records in clim: {}'.format(nr_trec_clim))

# Create clim file:
print('Creating file {}'.format(clm_file))
nc_out = netCDF4.Dataset(clm_file,'w', format='NETCDF4')
# Create dimensions:
for dim, the_dim in nc_in.dimensions.iteritems():
    # Found this code on the net... try it next time! MF: 30.7.2015
    nc_out.createDimension(dim, len(the_dim) if not the_dim.isunlimited() else None)
    
# Create global attributes:
att_dict = dict()
for att in nc_in.ncattrs():
    att_dict[att] = getattr(nc_in,att)
att_dict['remark1'] = 'Monthly climatology for run {}'.format(run_name)
att_dict['remark2'] = 'Number of time records used in original output files: {} out of {}'.format(nr_trec,trec)
att_dict['remark3'] = 'Script used for computing monthly clim: {}'.format(os.path.abspath(__file__))
# Write global atts:
nc_out.setncatts(att_dict)

# Now loop over variables:
for var, v_obj in nc_in.variables.items():
    # Skip variables time, time_step:
    if var == 'time_step' or var == 'time' or var == 'ocean_time':
        continue
    sys.stdout.write('{}\n'.format(var))
    sys.stdout.flush()
    if 'time' in v_obj.dimensions:
        dtuple = (nr_trec_clim,) + v_obj.shape[1:]
        # Allocate buffer for output (i.e. clim) variable:
        buf = np.zeros(list(dtuple), dtype=float)
        # Loop through filelist
        tctrl = 0 # Make sure only full years enter climatology
        for f in filelist:
            nc_in = netCDF4.Dataset(f,'r')
            print('Processing file: {}'.format(f))
            # Add up values for individual months:
            for t in range(len(nc_in.variables['time'])):
                if tctrl<nr_trec:
                    # Read time record t of variable var:
                    try:
                        v = nc_in.variables[var][t,:] # shape: (eta_rho,xi_rho) or (s_rho,eta_rho,xi_rho)
                    except KeyError: # Exception for TOT_PROD (pactcs304_avg does NOT contain summed productions yet!!!)
                        if var=='TOT_PROD':
                            tmp = nc_in.variables['PHOTOC_DIAT'][t,:]
                            tmp[tmp<0] = 0.0
                            v = tmp
                            tmp = nc_in.variables['PHOTOC_SP'][t,:]
                            tmp[tmp<0] = 0.0
                            v += tmp    
                            tmp = nc_in.variables['PHOTOC_DIAZ'][t,:]
                            tmp[tmp<0] = 0.0
                            v += tmp
                    # Number of month:
                    m = t//int(trec_per_month)
                    # Which month in year this time record belongs to (0=Jan, 1=Feb,...,11=Dec):
                    m_in_y = m%12
                    print('#trec = {}, #month = {}, month_in_yr = {}'.format(t,m,m_in_y))
                    # Add to corresponding month
                    buf[m_in_y,:] += v[:]
                    # t control (to make sure to only include full years)
                    tctrl +=1
        # Divide by number of time records for each month:
        buf = np.divide(buf,years*trec_per_month)
    else:
        # Time independent variable: just read values
        buf = nc_in.variables[var][:]

    # Define variable in output file:
    v_new = nc_out.createVariable(var, v_obj.datatype, v_obj.dimensions)
    # Copy variable attributes to output file:
    v_new.setncatts({k: v_obj.getncattr(k) for k in v_obj.ncattrs()})
    # Fill land points with fill value
#    if '_FillValue' in v_obj.ncattrs():
#        buf[buf==0.0] = getattr(v_obj,'_FillValue') 
    # Write values to output file:
    v_new[:] = buf
print('')

# Close NetCDF file
nc_out.close()
