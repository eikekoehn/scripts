#!/usr/bin/env python

'''
Module containing important functions for several purposes

Martin Frischknecht, 27.01.2014
'''
import numpy as np

def loadNetCDF(datafile,variable):
    '''
    Loads data from NetCDF file and returns it.
    
    Arguments:
    datafile: string indicating path and name of NetCDF file
    variable: string indicating variable within NetCDF to load
  
    Martin Frischknecht, 29.01.2014
    '''
    # Import python modules
    import netCDF4
    
    # Reading data from file
    f= netCDF4.Dataset(datafile,'r')
    data = f.variables[variable][:]
    
    return data

	
def points2dist(start_latlon,stop_latlon):
    '''
    Calculates distance on a sphere given start-point (lat/lon)
    and end-point (lat/lon). Returns distance between points in (m).
 
    Arguments:
    start_latlon: tuple of (lat, lon) of start point in degrees
    stop_latlon: tuple of (lat, lon) of end point in degrees
    
    Martin Frischknecht, 28.01.2014
    '''
    # Import python modules
    import math
    
    # Converting to radians and saving into single components
    start_lat = math.radians(start_latlon[0])
    start_lon = math.radians(start_latlon[1])
    stop_lat = math.radians(stop_latlon[0])
    stop_lon = math.radians(stop_latlon[1])
    
    # Calculating and defining needed quantities
    dLat = stop_lat - start_lat
    dLon = stop_lon - start_lon
    Re = 6371009 # in m[] from http://en.wikipedia.org/wiki/Earth_radius
    
    # Calculating Haversine formula from http://en.wikipedia.org/wiki/Haversine_formula
    a = math.sin(dLat/2)**2 + math.cos(start_lat) * math.cos(stop_lat) * math.sin(dLon/2)**2
    c = 2 * math.asin(math.sqrt(a))
    dist = Re * c

    return dist


def dist2points(start_latlon,bearing,distance):
    '''
    Calculates point on a sphere given a start-point (lat/lon),
    a bearing (due north), and a distance. Returns tuple of (lat, lon)
    for the calculated end-point in degrees.
  
    Arguments:
    start_latlon: tuple of (lat, lon) of start point in degrees
    bearing: angle due north in radians
    distance: distance to go from start-point in km
  
    Martin Frischknecht, 28.01.2014
    '''
    # Import python modules  
    import math
    
    stop_latlon = np.zeros([2])
    R = 6378.1 #Radius of the Earth
    brng = bearing
    d = distance
    
    lat1 = math.radians(start_latlon[0]) #Current lat point converted to radians
    lon1 = math.radians(start_latlon[1]) #Current long point converted to radians
    lat2 = math.asin( math.sin(lat1)*math.cos(d/R) + math.cos(lat1)*math.sin(d/R)*math.cos(brng))
    lon2 = lon1 + math.atan2(math.sin(brng)*math.sin(d/R)*math.cos(lat1),math.cos(d/R)-math.sin(lat1)*math.sin(lat2))

    stop_latlon[0] = math.degrees(lat2)
    stop_latlon[1] = math.degrees(lon2)

    return stop_latlon
	
 
def calcbearing(origin, destination):
    '''
    Calculates bearing (due north) of two points given by 
    tuples of (lat, lon). Returns bearing positive counter-clockwise.
    
    Arguments:
    origin: tuple of (lat, lon) of start points in degrees
    destination: tuple of (lat, lon) of end point in degrees
    
    Martin Frischknecht, 28.01.2014
    '''
    # Import python modules
    import math
    
    lat1, lon1 = origin
    lat2, lon2 = destination
    
    # Convert to radians
    rlat1 = math.radians(lat1)
    rlat2 = math.radians(lat2)
    dlon = math.radians(lon2-lon1)
    
    # Calculate bearing (return in radians and degrees)
    b = math.atan2(math.sin(dlon)*math.cos(rlat2),math.cos(rlat1)*math.sin(rlat2)-math.sin(rlat1)*math.cos(rlat2)*math.cos(dlon))
    bd = math.degrees(b)
    br,bn = divmod(bd+360,360) # the bearing remainder and final bearing
    
    return -b,-bd

	
def movingaverage(values, window):
    import numpy as np
    weights = np.repeat(1.0, window)/window
    return np.convolve(values, weights, 'valid')


def runningMean(x, N, mode='full'):
    import numpy as np
    #return np.convolve(x, np.ones((N,))/N, mode=mode)[(N-1):] # Old return: with default mode='full'
    return np.convolve(x, np.ones((N,))/N, mode=mode)


def uv2rho(u,v):
    '''
    Interpolates ROMS fields from u, v grids to
    rho points and returns the interpolated fields.
    
    Martin Frischknecht, Apr 2015
    '''
    import numpy as np
    # Read/Define shapes
    sz = list(u.shape)
    rho_shape = sz
    rho_shape[-1] = sz[-1] + 1
      
    # Initialize arrays
    u_rho = np.zeros(rho_shape)
    v_rho = np.zeros(rho_shape)
    
    # Regrid u/v on rho grid
    if len(sz)==4:
        # input fields have temporal/spatial dimension
        for t in range(sz[0]):
            print('Processing time step {}...'.format(t))
            for d in range(sz[1]):
                 # Interpolate u
                tmp = np.ma.masked_array(0.5*(u[t,d,:,0:-1] + u[t,d,:,1:]))
                u_rho[t,d,:,:] = np.ma.concatenate((u[t,d,:,0:1],tmp,u[t,d,:,-1:]),axis=1)
                # Interpolate v
                tmp = np.ma.masked_array(0.5*(v[t,d,0:-1,:] + v[t,d,1:,:]))
                v_rho[t,d,:,:] = np.ma.concatenate((v[t,d,0:1,:],tmp,v[t,d,-1:,:]),axis=0)               
      
    elif len(sz)==3:
        # input fields have temporal dimension
        for t in range(sz[0]):
            # Interpolate u
            tmp = np.ma.masked_array(0.5*(u[t,:,0:-1] + u[t,:,1:]))
            u_rho[t,:,:] = np.ma.concatenate((u[t,:,0:1],tmp,u[t,:,-1:]),axis=1)
            # Interpolate v
            tmp = np.ma.masked_array(0.5*(v[t,0:-1,:] + v[t,1:,:]))
            v_rho[t,:,:] = np.ma.concatenate((v[t,0:1,:],tmp,v[t,-1:,:]),axis=0)
    else:
        # Interpolate u
        tmp = np.ma.masked_array(0.5*(u[:,0:-1] + u[:,1:]))
        u_rho = np.ma.concatenate((u[:,0:1],tmp,u[:,-1:]),axis=1)
        # Interpolate v
        tmp = np.ma.masked_array(0.5*(v[0:-1,:] + v[1:,:]))
        v_rho = np.ma.concatenate((v[0:1,:],tmp,v[-1:,:]),axis=0)
    
    return u_rho, v_rho
    
def u2rho(u):
    '''
    Interpolates ROMS fields from u, v grids to
    rho points and returns the interpolated fields.
    
    Martin Frischknecht, Apr 2015
    '''
    import numpy as np
    
    # Initialize arrays
    u_shape = np.array(u.shape)
    u_shape[-1] = u_shape[-1] + 1
    u_rho = np.zeros(u_shape)
    
    # Regrid u/v on rho grid
    if u.ndim==4:
        # if input fields have temporal/vertical dimension
        for i in range(u.shape[0]):
            # Interpolate u
            tmp = np.ma.masked_array(0.5*(u[i,:,:,0:-1] + u[i,:,:,1:]))
            u_rho[i] = np.ma.concatenate((u[i,:,:,0:1],tmp,u[i,:,:,-1:]),axis=-1)
    elif u.ndim==3:
        # if input fields have temporal/vertical dimension
        tmp = np.ma.masked_array(0.5*(u[:,:,0:-1] + u[:,:,1:]))
        u_rho = np.ma.concatenate((u[:,:,0:1],tmp,u[:,:,-1:]),axis=-1)
    else:
        # Interpolate u
        tmp = np.ma.masked_array(0.5*(u[:,0:-1] + u[:,1:]))
        u_rho = np.ma.concatenate((u[:,0:1],tmp,u[:,-1:]),axis=1)
        
    del tmp
    return u_rho

def v2rho(v):
    '''
    Interpolates ROMS fields from u, v grids to
    rho points and returns the interpolated fields.
    
    Martin Frischknecht, Apr 2015
    '''
    import numpy as np
    
    # Initialize arrays
    v_shape = np.array(v.shape)
    v_shape[-2] = v_shape[-2] + 1
    v_rho = np.zeros(v_shape)
    
    # Regrid u/v on rho grid
    if v.ndim==4:
        # if input fields have temporal/vertical dimension
        for i in range(v.shape[0]):
            # Interpolate v
            tmp = np.ma.masked_array(0.5*(v[i,:,0:-1,:] + v[i,:,1:,:]))
            v_rho[i] = np.ma.concatenate((v[i,:,0:1,:],tmp,v[i,:,-1:,:]),axis=-2)
    elif v.ndim==3:
        # Interpolate v
        tmp = np.ma.masked_array(0.5*(v[:,0:-1,:] + v[:,1:,:]))
        v_rho = np.ma.concatenate((v[:,0:1,:],tmp,v[:,-1:,:]),axis=-2)
    else:
        # Interpolate v
        tmp = np.ma.masked_array(0.5*(v[0:-1,:] + v[1:,:]))
        v_rho = np.ma.concatenate((v[0:1,:],tmp,v[-1:,:]),axis=-2)
        
    del tmp
    return v_rho


def rotate_uv(u,v,angle):
    '''
    Rotates u,v components to zonal and meridional currents.
    Interpolate to rho grid first!
    
    Martin Frischknecht, Jan 2016
    '''
    import numpy as np

    if u.ndim==4 or v.ndim==4:
        # Prepare angle array  
        angle = np.array([angle]*u.shape[1])
        zonal = np.zeros(u.shape)
        meridional = np.zeros(u.shape)
        for t in range(u.shape[0]):
            # Rotate currents
            zonal[t] = u[t]*np.cos(angle) - v[t]*np.sin(angle)
            meridional[t] = u[t]*np.sin(angle) + v[t]*np.cos(angle)
    else:
        # Prepare angle array    
        if u.ndim==3 or v.ndim==3:
            angle = np.array([angle]*u.shape[0])
        # Rotate currents
        zonal = u*np.cos(angle) - v*np.sin(angle)
        meridional = u*np.sin(angle) + v*np.cos(angle)
    
    # Mask array
    if isinstance(u,np.ma.MaskedArray):
        zonal = np.ma.masked_where(u.mask==True,zonal)
        meridional = np.ma.masked_where(u.mask==True,meridional)
    elif isinstance(v,np.ma.MaskedArray):
        zonal = np.ma.masked_where(v.mask==True,zonal)
        meridional = np.ma.masked_where(v.mask==True,meridional)

    return zonal, meridional
    

def smooth(im, n=15):
    from scipy import signal
    """
    Smooth a 2D array im by convolving with a Gaussian kernel of size n
    Input:
    im (2D array): Array of values to be smoothed
    n (int) : number of points to include in the smoothing
    Output:
    improc(2D array): smoothed array (same dimensions as the input array)
    """
    g = gaussKern(n)
    improc = signal.convolve2d(im, g, mode='same', boundary='symm')
    return(improc)


def gaussKern(size):
    """
    Calculate a normalised Gaussian kernel to apply as a smoothing
    function.
    Input:
    size (int): the size of the kernel to use (how many points will be
                used in the smoothing operation).
    Output:
    g (array(size,size)): Normalised 2D kernel array for use in
                          convolutions
    """
    size = int(size)
    x,y = np.mgrid[-size:size+1,-size:size+1]
    g = np.exp(-(x**2/float(size)+y**2/float(size)))
    return g / g.sum()


def get_special_tracer(var,ncf):
    '''
    Get special tracer from ROMS output. Define your own if you need it!
    '''
    import netCDF4
    Q = 0.137
    
    if var=='TOT_PROD' or var=='NPP':
        tmp = ncf.variables['PHOTOC_SP'][:]
        tmp += ncf.variables['PHOTOC_DIAT'][:]
        tmp += ncf.variables['PHOTOC_DIAZ'][:]
        return tmp
    elif var=='DueBiol' or var=='NCP':
        Q = ncf.variables['Q'][:]
        DONrefract = ncf.variables['DONrefract'][:]
        bgc_fluxes = ['DENITRIF', 'SP_NO3_UPTAKE', 'DIAT_NO3_UPTAKE', 'DIAZ_NO3_UPTAKE',
                      'SP_NH4_UPTAKE', 'DIAT_NH4_UPTAKE', 'DIAZ_NH4_UPTAKE', 'DON_REMIN', 
                      'DONr_REMIN', 'ZOO_LOSS_DIC', 'SP_LOSS_DIC', 'DIAT_LOSS_DIC', 
                      'DIAZ_LOSS_DIC', 'SP_GRAZE_DIC', 'DIAT_GRAZE_DIC', 'DIAZ_GRAZE_DIC',
                      'POC_REMIN']
        factors = {'DENITRIF':-1., 'SP_NO3_UPTAKE':-1., 'DIAT_NO3_UPTAKE':-1., 
                   'DIAZ_NO3_UPTAKE':-1., 'SP_NH4_UPTAKE':-1., 'DIAT_NH4_UPTAKE':-1., 
                   'DIAZ_NH4_UPTAKE':-1., 'DON_REMIN':1., 'DONr_REMIN':1., 'ZOO_LOSS_DIC':Q, 
                   'SP_LOSS_DIC':Q, 'DIAT_LOSS_DIC':Q, 'DIAZ_LOSS_DIC':Q, 'SP_GRAZE_DIC':Q,
                   'DIAT_GRAZE_DIC':Q, 'DIAZ_GRAZE_DIC':Q, 'POC_REMIN':Q*(1-DONrefract)}
        # If NCP, dont consider NO3 loss through DENITRIF
        if var=='NCP':
            del bgc_fluxes[0]
            del factors['DENITRIF']
        # Loop through components and read them
        tmp = None
        for key in bgc_fluxes:
            print('Reading {} to compute {}...'.format(key,var))
            if tmp==None:
                tmp = ncf.variables[key][:]*factors[key]
            else:
                tmp += ncf.variables[key][:]*factors[key]
        return tmp
    elif var=='TIN':
        tmp = ncf.variables['NO3'][:]
        tmp += ncf.variables['NH4'][:]
        return tmp
    elif var=='TON':
        tmp = ncf.variables['SPC'][:]*Q
        tmp += ncf.variables['DIATC'][:]*Q
        tmp += ncf.variables['DIAZC'][:]*Q
        tmp += ncf.variables['ZOOC'][:]*Q
        tmp += ncf.variables['DON'][:]
        tmp += ncf.variables['DONR'][:]
        return tmp
    elif var=='TOT_PHYTO':
        tmp = ncf.variables['SPC'][:]
        tmp += ncf.variables['DIATC'][:]
        tmp += ncf.variables['DIAZC'][:]
        return tmp
    elif var=='Nstar':
        tmp = ncf.variables['NO3'][:]
        tmp += ncf.variables['NH4'][:]
        tmp -= 16*ncf.variables['PO4'][:]
        tmp += 2.9
        return tmp
    elif var=='tau':
#        f = '/net/kryo/work/martinfr/Roms/Output/pactcs30/hcBEC_T002_pactcs30/avg/z_pactcs30_avg_1979-2015_clim_phys_residence_time.nc'
        f = '/net/kryo/work/martinfr/Roms/Output/pactcs30/hcBEC_T002_pactcs30/avg/z_pactcs30_avg_1979-2015_clim_phys_residence_time_newnew.nc'
        ncf = netCDF4.Dataset(f,'r')
        tmp = ncf.variables['tau'][:]
        return tmp
    elif var=='TN':
        tmp = ncf.variables['NO3'][:]
        tmp += ncf.variables['NH4'][:]
        tmp += ncf.variables['SPC'][:]*Q
        tmp += ncf.variables['DIATC'][:]*Q
        tmp += ncf.variables['DIAZC'][:]*Q
        tmp += ncf.variables['ZOOC'][:]*Q
        tmp += ncf.variables['DON'][:]
        tmp += ncf.variables['DONR'][:]
        return tmp       
    elif var=='NO3_UPTAKE':
        tmp = ncf.variables['SP_NO3_UPTAKE'][:]
        tmp += ncf.variables['DIAT_NO3_UPTAKE'][:]
        tmp += ncf.variables['DIAZ_NO3_UPTAKE'][:]
        return tmp   
    elif var=='NH4_UPTAKE':
        tmp = ncf.variables['SP_NH4_UPTAKE'][:]
        tmp += ncf.variables['DIAT_NH4_UPTAKE'][:]
        tmp += ncf.variables['DIAZ_NH4_UPTAKE'][:]
        return tmp   

  
