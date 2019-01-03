#!/usr/bin/python

# Computes a monthly climatology over one single netCDF file (e.g. 1 file containing 10 years of output > 120 months).
# Martin Frischknecht, Mar 2015 (modified from Damian Loher)

import netCDF4
import numpy as np
import sys
import os

# Where is the data stored?
host = 'kryo'
# ROMS setup
setup = 'pactcs30'
# Run
run_name = 'newsrc_monthly_{}'.format(setup)
job = '4'
dir_in = '/net/{}/work/martinfr/Roms/Output/{}/{}/avg'.format(host,setup,run_name)
dir_out = dir_in


# Input avg file:
avg_file = '{}/{}{}_avg.nc'.format(dir_in,setup,job)
print('Using avg file {}'.format(avg_file))

# Output file (clm):
clm_file = '{}/{}{}_avg_clim.nc'.format(dir_out,setup,job)

# Find out number of time records:
nc_in = netCDF4.Dataset(avg_file, 'r')
trec = len(nc_in.dimensions['time'])
navg = nc_in.getncattr('navg')
dt = float(nc_in.getncattr('dt'))
print('Number of time records used: {}'.format(trec))
print('dt = {}'.format(dt))
print('navg = {}'.format(navg))

# Number of days between two consecutive time records:
days_per_trec = float(dt*navg)/86400
# Number of days per month:
days_per_month = days_per_trec
# Number or time records per month and year:
trec_per_month = days_per_month/days_per_trec
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

# Number of time records used for computing climatology:
nr_trec = years*12*trec_per_month
# Number of time records in clim (output) file:
nr_trec_clim = 12
print('number of time records in clim: {}'.format(nr_trec_clim))

# Create clim file:
print('Creating file {}'.format(clm_file))
nc_out = netCDF4.Dataset(clm_file,'w', format='NETCDF4')
# Create dimensions:
for dim, the_dim in nc_in.dimensions.iteritems():
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
        # Add up values for individual months:
        for t in range(int(nr_trec)):
            # Read time record t of variable var:
            v = nc_in.variables[var][t,:] # shape: (eta_rho,xi_rho) or (s_rho,eta_rho,xi_rho)
            # Number of month:
            m = t//int(trec_per_month)
            # Which month in year this time record belongs to (0=Jan, 1=Feb,...,11=Dec):
            m_in_y = m%12
            print('#trec = {}, #month = {}, month_in_yr = {}'.format(t,m,m_in_y))
            # Add to corresponding month
            buf[m_in_y,:] += v[:]
        # Divide by number of time records for each month:
        buf = np.divide(buf,years*trec_per_month)
    else:
        # Time independent variable: just read values
        buf = nc_in.variables[var][:]

    # Define variable in output file:
    if '_FillValue' in v_obj.ncattrs():
        v_new = nc_out.createVariable(var,v_obj.datatype,v_obj.dimensions,fill_value=getattr(v_obj,'_FillValue'))
        buf[buf==0.0] = getattr(v_obj,'_FillValue') 
    else:
        v_new = nc_out.createVariable(var,v_obj.datatype,v_obj.dimensions)
    
    # Copy variable attributes to output file:
    v_att_dict = {}
    for att in v_obj.ncattrs():
        if att=='_FillValue':
            pass
        else:
            v_att_dict[att] = getattr(v_obj,att)
    v_new.setncatts(v_att_dict)
    # Write values to output file:
    v_new[:] = buf
print('')

# Close NetCDF file
nc_out.close()
