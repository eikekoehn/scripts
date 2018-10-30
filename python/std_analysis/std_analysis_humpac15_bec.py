#!/usr/bin/env python

# Script from Domitille (09.04.2018)

# Make sure all modules are found:
import sys
import os

os.environ["PYNGL_NCARG"] = "/usr/local/up/lib/python2.7/site-packages/PyNGL/ncarg"
sys.path.append('/home/loher/python/GenPlots')
sys.path.append('/home/loher/python/modules')

import gen_plots
import numpy as np

case1 = 'humpac15_bec_1979_2016_glodap_eq_grd'
#case1='sa_amz_full_dda_51'
#infile_bgc_flux1 = '/net/hydro/work/up/turig/comp_frc_files/frc_test_all/calcs_15km_comp_frc_files_era2005_monthly_Y6_Y10_bgc_flux_avg.nc'
#infile1='/net/kryo/work2/dlouchar/ROMS/outputs/amacan50/t_amacan50_2/avg/t_amacan50_2_00004.nc'
#infile1='/net/kryo/work2/fana/hindcast_bec/hindcast/output/avg/avg_1980.nc'
infile1='/net/kryo/work2/fana/hindcast_bec/analysis/humpac15_avg_1979-2016_clim.nc'


#ROMS_clim_1='/net/kryo/work2/dlouchar/ROMS/inputs/large/amaze50/amaze50_clm2.nc'
#ROMS_clim_1 = '/net/kryo/work2/dlouchar/ROMS/inputs/amaze50_DDAsg/amaze50sg/amaze50sg_clm2_stan.nc'
#ROMS_clim_1='/net/kryo/work2/dlouchar/ROMS/inputs/large/amacan50/amacan50_clm.nc'
#ROMS_clim_1='/net/kryo/work/fana/Roms_Tools_Outputs_humpac/humpac15_clim_BEC_normal_yr_365/humpac15_clm.nc'
#ROMS_clim_2='/net/kryo/work/fana/Roms_Tools_Outputs_humpac/humpac15_clim_BEC_normal_yr_365/humpac15_clm2.nc'
ROMS_clim_1 = '/net/kryo/work2/fana/romstools_outputs_humpac15/BEC_files_glodapv2_rf015/humpac15_clm2.nc'

#grdfile='/net/kryo/work2/dlouchar/ROMS/inputs/large/amaze50/inputs_14-03-18/amaze50_grd.nc'
#grdfile='/net/kryo/work2/dlouchar/ROMS/inputs/large/amacan50/amacan50_grd.nc'
#grdfile='/net/kryo/work2/fana/romstools_outputs_humpac15/humpac15_grd.nc'
grdfile='/net/nardus/work/Ana/masks_and_maps/humpac_grids/humpac15_sm5rf015_grd.nc' # glodap
#grdfile='/net/nardus/work/Ana/masks_and_maps/humpac_grids/humpac15_grd.nc'

# Comments about the simulations (compilers and compilations options, ...):
#comment_str = 'US West Coast, 15km resolution<br>Case1: run with ERA interim forcing<br>'\
#                'Case2: COADS QuikSCAT forcing<br>'
#comment_str = 'US West Coast, 15km resolution<br>Case1: run with ERA interim forcing<br>'\
#                 'Case2: control run'
#comment_str = 'US West Coast, 15km resolution<br>Case1: control run<br>'\
#                'Case2: COADS QuikSCAT forcing<br>'
#comment_str = 'US West Coast, 15km resolution<br>Case1: control run (old COADS QuikScat forcing<br>'\
#              'Case2: ERA interim 2005 monthly forcing<br>'
#comment_str = 'US West Coast, 15km resolution<br>Case1: ERA forcing from Ivy<br>'\
#              'Case2: ERA interim 2005 monthly forcing<br>'
#comment_str = 'Grid testing<br>'\
comment_str = 'Comparison of Year 1979-2016 of Humpac15-BEC run with climatology'

# Construct class:
x = gen_plots.Gen_Plots(case1,grdfile)

# Set label to used at the bottom of the plots:
#x.set_bottom_plotlabel('New grid (Amacan50)')
x.set_bottom_plotlabel('New grid (Humpac15)')

# Simulation years contained in infile:
x.yrstr = [ 'DJF','MAM','JJA','SON','AnnMean']
#x.yrstr = ['AnnMean']
# Set input file for time intervals:
infile_dict1 = {'1979-2016': infile1 }

# Directory where plots + webpage will be put:
#x.outdir = '/net/malva/work/loher/plots/%s' %(case1)
x.outdir = '/home/koehne/Documents/scripts/python/%s' %(case1)

x.startrec = [ -1, 2, 5, 8, 0]  # first time record (first record = 0)
x.endrec =   [ 2, 5, 8, 11, 12]  # 12 last time record + 1(!)
#x.startrec = [0]
#x.endrec = [12]

#x.startrec_obs = [ -1, 2, 5, 8, 0]  # first time record (first record = 0)
#x.endrec_obs = [ 2, 5, 8, 11, 12]   # last time record + 1
x.set_default_start_rec_obs_3d([ -1, 2, 5, 8, 0 ])
x.set_default_end_rec_obs_3d([ 2, 5, 8, 11, 12 ]) #12
x.set_default_start_rec_obs_2d([ -1, 2, 5, 8, 0 ])
x.set_default_end_rec_obs_2d([ 2, 5, 8, 11, 12 ]) #12
x.set_infile(infile1)

# Title to be used in Web page:
x.sim_title = 'Comparison %s vs. Observations' %(case1)
x.sim_comment = comment_str

# Generate one table for each time period:
x.tvert = True
# Type of sigma levels used in ROMS:
x.siglev=3;

# Add some tracers: use add_tracer for 3D tracers
#x.add_tracer_top('pCO2','pCO2',infile_bgc_flux1,infile_bgc_flux2,'pCO2','pCO2','ppm',1,'rho')
x.add_tracer('Nitrate','NO3',infile_dict1,'NO3','mmol m-3',1,'rho')
#x.add_tracer('Phosphate','PO4',infile_dict1,'PO4','mmol m-3',1,'rho')
#x.add_tracer('Iron','Fe',infile_dict1,'Fe','mmol m-3',1,'rho')
#x.add_tracer('SiO3','SiO3',infile_dict1,'SiO3','mmol m-3',1,'rho')
#x.add_tracer_top('MLD','MLD',infile_dict1,'hbls','mmol m-3',1,'rho')
#x.add_tracer('Total_Chlorophyll','Tot_CHL',infile_dict1,'','mg Chl m-3',1,'rho')
#x.add_tracer_top('pH','pH',infile_bgc_flux1,infile_bgc_flux2,'pH','pH','',1,'rho')

# Comparison with observation: the following arguments are required:
#
#  - name of tracer
#  - Netcdf file of observational data
#  - name of variable to use in Netcdf file with obs. data
#  - conversion factor (needed for converting obs. data to units of model results)
#

#x.add_obs('Temperature','/net/kryo/work/fana/Roms_Tools_Outputs_humpac/humpac15_clim_BEC_normal_yr_365/humpac15_clm.nc',\
#          'temp','WOA 2013',1)
#x.add_obs('Salinity','/net/kryo/work/fana/Roms_Tools_Outputs_humpac/humpac15_clim_BEC_normal_yr_365/humpac15_clm.nc',\
#          'salt','WOA 2013',1)
#x.add_obs('Temperature','/net/hydro/work/loher/lunaria_data/usw15km/obs/comp_dat/woa05_usw15km.nc',\          
#          'TEMP_WOA05','WOA 2005',1)
#x.add_obs('Salinity','/net/hydro/work/loher/lunaria_data/usw15km/obs/comp_dat/woa05_usw15km.nc',\
#          'SALT_WOA05','WOA 2005',1)
#x.add_obs('Temperature','/net/kryo/work2/dlouchar/ROMS/inputs/large/amaze50/amaze50_clm.nc',\
#          'temp','WOA 2013',1)
#x.add_obs('Salinity','/net/kryo/work2/dlouchar/ROMS/inputs/large/amaze50/amaze50_clm.nc',\
#          'salt','WOA 2013',1)
#x.add_obs('Sea Surface Height','/net/kryo/work2/dlouchar/ROMS/inputs/large/amaze50/amaze50_clm.nc',\
#          'zeta','WOA 2013',1) 
#x.add_obs('MLD','/net/kryo/work2/dlouchar/ROMS/inputs/large/amaze50/amaze50_clm.nc',\
#          'mld','MLD Argo',1)  
#x.add_obs('Temperature','/net/kryo/work2/dlouchar/ROMS/inputs/amaze50_DDAsg/amaze50sg/originalfiles/amaze50sg_clm.nc',\
#          'temp','WOA 2013',1)
#x.add_obs('Temperature','/net/kryo/work2/dlouchar/ROMS/inputs/amaze50_DDAsg/amaze50sg/originalfiles/amaze50sg_clm.nc',\
#          'temp','WOA 2013',1)
#x.add_obs('Salinity','/net/kryo/work2/dlouchar/ROMS/inputs/amaze50_DDAsg/amaze50sg/originalfiles/amaze50sg_clm.nc',\
#          'salt','WOA 2013',1)
#x.add_obs('Sea Surface Height','/net/kryo/work2/dlouchar/ROMS/inputs/amaze50_DDAsg/amaze50sg/originalfiles/amaze50sg_clm.nc',\
#          'zeta','WOA 2013',1) 
#x.add_obs('MLD','/net/kryo/work2/dlouchar/ROMS/inputs/amaze50_DDAsg/amaze50sg/originalfiles/amaze50sg_clm.nc',\
#          'mld','MLD Argo',1) 
        
#x.add_obs('Oxygen','/net/hydro/work/loher/lunaria_data/usw15km/obs/comp_dat/woa05_usw15km.nc',\
#          'O2_WOA05','WOA 2005',1000.0 / 22.4)
#x.add_obs('Chlorophyll_A','/net/hydro/work/loher/lunaria_data/usw15km/obs/chl_USWC_00_06.nc',\
#          'CHL_00_06','SeaWiFS 2000-2006',1)

# Add ROMS climatology for comparison with observations:
x.add_obs_from_ROMS_clim(ROMS_clim_1)
#x.add_obs_from_ROMS_clim(ROMS_clim_2)
# Flag for switching comparisons with observations on/off
x.comp_obs = True

#
# Set plot resources:
#
hres = x.get_hmres()
hres.cnLinesOn = False # no contour lines
#hres.mpGridAndLimbOn = False # no grid lines
# Higher resolution for coast line:
hres.mpDataBaseVersion = 'Ncarg4_1'
hres.mpDataSetName = 'Earth..1'
# Draw tick marks on x and y axis:
hres.tmXBOn = True
hres.tmYLOn = True
hres.tmXBMajorLengthF = 0.01 # length of tick marks
##
## Zoom in for horizontal sections: uncomment to plot the whole domain
##
## Set domain to be covered by plot: specify lower
## left and upper right corner:
#hres.mpLimitMode    = "LatLon"  # Limit portion of map that is viewed.
#hres.mpMinLatF      = -20
#hres.mpMaxLatF      = 5
# Note that with the SAWC grid, longitudes have to be given as deg. east:
#hres.mpMinLonF      = 265
#hres.mpMaxLonF      = 290
x.set_hmres(hres)

x.fill_land = True

# Scaling factor for conversion to PNG:
x.png_resize = 22
#x.add_comp_obs_atdepth('Chlorophyll_A',[0]) # compare to observations only at surface
#x._hplot_lev['Chlorophyll_A'] = np.arange(0,2,0.1)
#x._hplot_lev['Sea_Surface_Height'] = np.arange(-30,30,2)
#x._hplot_lev['Temperature'] = np.arange(10,20,0.25)
#x._hplot_lev_diff['Temperature'] = np.arange(-4,4,0.25)
#x._hplot_lev_diff_obs['Temperature'] =np.arange(-12,12,0.25)
#x._hplot_lev_diff_obs['Salinity'] =np.arange(-4,4,0.25)
#x._hplot_lev['Salinity'] =np.array([5,10,20,30,32,34,34.5,35,35.5,36,36.5,37,37.5,38])
#x._hplot_lev['Iron'] =np.arange(0,0.004,0.0005)
#x._hplot_lev_diff_obs['Iron'] =np.arange(-0.004,0.004,0.00025)
#x._hplot_lev['SiO3'] =np.array([0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5.5,6,10,20,30,40])

#x.set_hdepths([ 0, 100 ])  # depths in meters for horiz. sections
x.set_hdepths([ 0, 100, 200 ])
x.set_longitudes([ '100W', '80W' ],[[],[]]) # longitudes for sections along y axis
x.set_latitudes([ '5S', '12S' ],[[],[]])
x.set_vdepths([ 200]) # generate additionally section of top 200m


# Do plots:
x.horiz_sections()
x.vert_sections_x()
x.vert_sections_y()

# Convert plots to PNG and create webpage:
x.convert_plots_png()
x.gen_webpage()
