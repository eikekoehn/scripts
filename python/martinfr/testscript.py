# -*- coding: utf-8 -*-
"""
Testfile to see how the scripts written by Martin Frischknecht work.

Created on Wed Nov  7 11:42:51 2018

@author: koehne
"""

from myromstools import getGrid
from myplotlib import *

setup = 'pactcs60'

# Set grid/clm file
grdfile = '/net/kryo/work/martinfr/Roms/Inputs/{}/{}_grd.nc'.format(setup,setup)
clmfile = '/net/kryo/work/martinfr/Roms/Inputs/{}/N64ts10tb4hc250_grd_merged_SiO3_PO4_fix/{}_clm.nc'.format(setup,setup)

# Read grid
romsGrd = getGrid(grdfile)
# Add necessary attributes
romsGrd.getAttrs(clmfile)
# Add lat lon
romsGrd.getLatLon()
# Add grid area
romsGrd.getArea()
# Add angle
romsGrd.getAngle()
# Add topography
romsGrd.getTopo()

easyplot(romsGrd.h)