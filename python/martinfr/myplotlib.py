#!/usr/bin/env python

'''
Module containing the my main plotting tools.

Martin Frischknecht, 15.07.2014
'''

#########################################################
# Import modules shared by most functions
#########################################################

import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from pylab import *

#########################################################
# General methods for plot design/cosmetics
#########################################################

def mycolmap(ncolr=16,name='jet'):
    return plt.get_cmap(name, ncolr)

def mycolmap_diff(ncolr=16,name='RdBu_r'):
    return plt.get_cmap(name, ncolr)

def monthnames(option=0):
    '''
    Return month names either
    1: abbrev. not capital
    2: abbrev. capital
    3: long names
    '''
    if option==0:
        monthnames = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    elif option==1:
        monthnames = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
    elif option==2:
        monthnames = ['January','February','March','April','May','June','July','August','September',
                      'October','November','December']
    elif option==3:
        monthnames = ['J','F','M','A','M','J','J','A','S','O','N','D']   
    # Return list
    return monthnames

def pdf_settings():
    from IPython.display import set_matplotlib_formats
    set_matplotlib_formats('png', 'pdf')
    matplotlib.rcParams['pdf.fonttype'] = 42

def rcParams(fontsize=10,elw=1.0):
    '''
    Set defaults rcParams for my figures, axes
    and plots.
    '''
    # Change matplotlib rcParams
    plt.rcParams['figure.facecolor'] = 'w'
    plt.rcParams['xtick.color'] = '0.15'
    plt.rcParams['xtick.labelsize'] = fontsize
    plt.rcParams['ytick.color'] = '0.15'
    plt.rcParams['ytick.labelsize'] = fontsize
    plt.rcParams['text.color'] = '0.15'
    plt.rcParams['font.size'] = fontsize
    plt.rcParams['grid.color'] = '0.5'
    plt.rcParams['patch.edgecolor'] = 'white'
    plt.rcParams['axes.edgecolor'] = '0.15'
    plt.rcParams['axes.linewidth'] = elw
    plt.rcParams['axes.labelsize'] = fontsize
    plt.rcParams['xtick.direction'] = 'out'
    plt.rcParams['ytick.direction'] = 'out'


def no_spines(ax,top=True,bottom=True,right=True,left=True):
    '''
    Remove unnecessary spines from figure axes.
    '''
    # Plot cosmetics
    if top==True:
        ax.spines["top"].set_visible(False)
    if bottom==True:
        ax.spines["bottom"].set_visible(False)
    if right==True:
        ax.spines["right"].set_visible(False)
    if left==True:
        ax.spines["left"].set_visible(False)


def no_ticks(ax,axis='both',which='both',bottom='on',top='on',labelbottom='off',labelleft='off'):
    plt.tick_params(
    axis=axis,          # changes apply to the x-axis
    which=which,      # both major and minor ticks are affected
    bottom=bottom,      # ticks along the bottom edge are off
    top=top,         # ticks along the top edge are off
    labelbottom=labelbottom,
    labelleft=labelleft)


def add_cbar(im=None,cax=None,units='',fontsize=10,elw=1.0,orientation='vertical',**kwargs):
    '''
    Add colorbar to plot. 
    
    cax: colorbar axes instance (default None, cax is added)
    im: mappable object
    units: colorbar labelling
    '''
    
    # Make colorbar
    if cax==None:
        if im==None:
            cbar = plt.colorbar(orientation=orientation,**kwargs)
        else:
            cbar = plt.colorbar(mappable=im,orientation=orientation,**kwargs)
    else:
        if im==None:
            cbar = plt.colorbar(cax=cax,orientation=orientation,**kwargs)
        else:
            cbar = plt.colorbar(cax=cax,mappable=im,orientation=orientation,**kwargs)
    
    # Change settings
    rcParams(fontsize=fontsize,elw=elw)
    if orientation=='horizontal':
        cbar.ax.set_xlabel(units,fontsize=fontsize)
        plt.setp(plt.getp(cbar.ax.axes, 'xticklines'), visible=False)
    else:
        cbar.ax.set_ylabel(units)
        plt.setp(plt.getp(cbar.ax.axes, 'yticklines'), visible=False)
    cbar.solids.set_edgecolor("face")
    cbar.outline.set_visible(False)    
    plt.setp(plt.getp(cbar.ax.axes, 'yticklines'), visible=False)
    
    return cbar


def bounds(scalar_fields):
    '''
    Get the bounds of a set of scalar_fields
    :param scalar_fields : the scalar field set
    :return: a set of normalized vector field components
    '''
    max_bound = -np.inf
    min_bound =  np.inf
    
    for scalar_field in scalar_fields:
        scalar_field = np.ma.masked_invalid(scalar_field)
        try:
            max_lim, min_lim = np.nanpercentile(scalar_field.filled(np.nan), [95,5])
        except:
            max_lim, min_lim = -np.inf, np.inf
        if max_lim > max_bound:
            max_bound = max_lim
        if min_lim < min_bound:
            min_bound = min_lim
    
    return min_bound, max_bound


def align_yaxis(ax1, v1, ax2, v2):
    """adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1"""
    _, y1 = ax1.transData.transform((0, v1))
    _, y2 = ax2.transData.transform((0, v2))
    adjust_yaxis(ax2,(y1-y2)/2,v2)
    adjust_yaxis(ax1,(y2-y1)/2,v1)

def adjust_yaxis(ax,ydif,v):
    """shift axis ax by ydiff, maintaining point v at the same location"""
    inv = ax.transData.inverted()
    _, dy = inv.transform((0, 0)) - inv.transform((0, ydif))
    miny, maxy = ax.get_ylim()
    miny, maxy = miny - v, maxy - v
    if -miny>maxy or (-miny==maxy and dy > 0):
        nminy = miny
        nmaxy = miny*(maxy+dy)/(miny+dy)
    else:
        nmaxy = maxy
        nminy = maxy*(miny+dy)/(maxy+dy)
    ax.set_ylim(nminy+v, nmaxy+v)


class set_cbar_zero(Normalize):
    """
    set_cbar_zero(midpoint = float)       default: midpoint = 0.

    Normalizes and sets the center of any colormap to the desired value which 
    is set using midpoint. 
    """
    
    def __init__(self, vmin=None, vmax=None, midpoint=0., clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        if isinstance(value, ma.MaskedArray):
            mask = value.mask
            value = value.filled(self.vmax)
        else:
            mask = False
    
            x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
            return np.ma.masked_array(np.interp(value, x, y),mask=mask)

#########################################################
# Plotting tools
#########################################################

def easyplot(data,ax=None,title='',xlab='',ylab='',units='',save=0,output='',**kwargs):
    '''
    Plots data convienient and fast.
    Return plot for a fast and first 
    look at data
    
    Arguments:
    data: data array of no specific length
          (either 1d or 2d)
     
    Martin Frischknecht, 29.01.2014
    '''

    # Check data structure (list,array) and dimensions
    try:
        if data.ndim>2:
            print('Warning: easyplot only accepts 1D/2D data structures...')
            return
    except:
        data = np.array(data)
        print('Converted data structure into an array...')
        if data.ndim>2:
            print('Warning: easyplot only accepts 1D/2D data structures...')
            return
    
    # If no axes given, field is plotted to new figure
    # Figure will be returned or saved then.
    if not ax:
        # Initiate new figure
        plt.figure(facecolor='w')
        ax = plt.subplot(111)
        # return or save figure
        returnfig = 1
    else:
        # plot to existing axes, do not return figure
        returnfig = 0 
    
    # General axes params
    rcParams()
    no_spines(ax)
  
    # Check for data dimensions
    sz = data.shape
  
    # Define color map
    colmap = plt.get_cmap('RdBu_r', 16)
      
    # Plot data
    if len(sz) == 1:
        plt.plot(data,**kwargs)
        plt.grid()
        plt.xlim = (0,sz[0])
        try:
            plt.ylim = (data.min()-0.05*(data.max()-data.min()),
                        data.max()+0.05*(data.max()-data.min()))
        except:
            pass
    elif len(sz) == 2:
        plt.pcolormesh(data,cmap=colmap,**kwargs)
        plt.grid()
        ax.set_xlim(0,sz[1])
        ax.set_ylim(0,sz[0])
        add_cbar(units=units)

    # Axis labelling
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(title)
    
    # If new figure was plotted
    if returnfig:
        # Save or show figure
        if save:
            plt.savefig(output,dpi=220)
            plt.close()
        else:
            plt.show(block=False)
        

def plot_pactcs30(data,ax=None,title='',zoom=0,save=0,output='',cmap='RdBu_r',units='',**kwargs):
    '''
    Plots data on pactcs30 grid.
    
    Arguments:
    data: array of the size [518,604], one layer of a certain data set
    title: optional, default is 'Grid pacsg'
     
    Martin Frischknecht, Oct 2015
    '''
    from mpl_toolkits.basemap import Basemap
      
    # Load ROMS pacsg grid information
    fgr =  '/net/kryo/work/martinfr/Roms/Inputs/pactcs30/pactcs30_grd.nc'
    #fgr =  'Data/pactcs30_grd.nc'
    ncfgr = netCDF4.Dataset(fgr, 'r')
    lon = ncfgr.variables['lon_rho'][:]
    lat = ncfgr.variables['lat_rho'][:]
    
    # Define subset of map (zoom or not):
    if not zoom:
        # Lon/Lat limits
        lonmin, lonmax = 100, 305
        latmin, latmax = -80, 73
        res='c'
            
        xticks = np.arange(120.,310.,30.)
        xtick_labels_W = [ '{:g}W'.format(t) for t in np.arange(150.,55.,-30.)]
        xtick_labels_E = [ '{:g}E'.format(t) for t in np.arange(120.,190.,30.)]
        xtick_labels = xtick_labels_E + xtick_labels_W
    
        yticks = np.arange(-60.,70.,20.)
        ytick_labels_N = [ '{:g}N'.format(t) for t in np.arange(0.,75.,20.)]
        ytick_labels_S = [ '{:g}S'.format(t) for t in np.arange(60.,0.,-20.)]
        ytick_labels = ytick_labels_S + ytick_labels_N
      
    if zoom==1:
        # Lon/Lat limits
        lonmin, lonmax = 223, 248
        latmin, latmax = 27, 52
        res='l'
        
        # Compute ticks on x and y axis:
        xticks = np.arange(225.,250.,5.)
        xtick_labels = [ '{:g}W'.format(t) for t in np.arange(135.,110.,-5.)]
        yticks = np.arange(30.,52.,5.)
        ytick_labels = [ '{:g}N'.format(t) for t in np.arange(30.,52.,5.)]

    if zoom==2:
        # Lon/Lat limits
        lonmin, lonmax = 213, 255
        latmin, latmax = 20, 60
        res='l'
        
        # Compute ticks on x and y axis:
        xticks = np.arange(215.,255.,5.)
        xtick_labels = [ '{:g}W'.format(t) for t in np.arange(145.,105.,-5.)]
        yticks = np.arange(25.,62.,5.)
        ytick_labels = [ '{:g}N'.format(t) for t in np.arange(25.,62.,5.)]
        
    if zoom==3:
        # Lon/Lat limits
        lonmin, lonmax = 225, 255
        latmin, latmax = 23, 57
        res='l'
        
        # Compute ticks on x and y axis:
        xticks = np.arange(230.,255.,5.)
        xtick_labels = [ '{:g}W'.format(t) for t in np.arange(145.,105.,-5.)]
        yticks = np.arange(25.,62.,5.)
        ytick_labels = [ '{:g}N'.format(t) for t in np.arange(25.,62.,5.)]
      
        
    # Set up map
    m = Basemap(llcrnrlon=lonmin,llcrnrlat=latmin,urcrnrlon=lonmax,urcrnrlat=latmax,
            resolution=res,projection='cyl',lon_0=200,lat_0=0.)
    
    x, y = m(lon, lat)
    
    # Convert ticks to map coordinates:
    xticks_map, _ = m( xticks,latmin*np.ones(len(xticks)) )
    _, yticks_map = m( lonmin*np.ones(len(yticks)),yticks )
      
    # Define color map
    if isinstance(cmap, basestring):
        colmap = plt.get_cmap(cmap, 16)
    else:
        colmap = cmap
         
    # If no axes given, field is plotted to new figure
    # Figure will be returned or saved then.
    if not ax:
        # Initiate new figure
        plt.figure(facecolor='w')
        ax = plt.subplot(111)
        # return or save figure
        returnfig = 1
    else:
        # plot to existing axes, do not return figure
        returnfig = 0 
      
    # General axes params    
    rcParams()
    
    # Plot data
    m.pcolormesh(x,y,data,cmap=colmap,**kwargs)
    add_cbar(units=units)
    
    # Basemap bakcground methods
    m.fillcontinents(color='0.65')
    m.drawmapboundary(color='w',linewidth=0)
    m.drawcoastlines()
    
    # Set ticks
    plt.xticks(xticks_map[::2],xtick_labels[::2])
    plt.yticks(yticks_map[::2], ytick_labels[::2])  
    
    # Title
    plt.title(title, fontsize=12)

    # If new figure was plotted
    if returnfig: 
        # Save or show figure
        if save:
            plt.savefig(output,dpi=220)
            plt.close()
        else:
            plt.show(block=False)


  
def plot_vector_pacsg(xi_component,eta_component):
  '''
  Plots wind fields as vector fields on pacsg using a 
  Mercator map projection.
  
  Arguments:
  xi_component: xi_component of vector (shape [348,418])
  eta_component: eta_component of vector (shape [348,418])

  IMPORTANT:
  Basemap.quiver assumes vectors to be in East-, North-, directions.
  Vectors from pacsg oriented along xi, eta axis of grid must therefore
  be transformed before plotting results.

  Vector field needs to be interpolated to rho-points!

  Martin Frischknecht, 27.01.2014
  '''
  #########################################################
  # Load requirements
  #########################################################

  # Import python modules
  import numpy as np
  from mpl_toolkits.basemap import Basemap
  import matplotlib.pyplot as plt
  import netCDF4

  # Load ROMS pacsg grid data
  grdfile = '/net/kryo/work/munnich/Roms/Inputs/gridfiles/pacsg_grd.nc'
  fgr = netCDF4.Dataset(grdfile,'r')
  mask = fgr.variables['mask_rho']
  angle = fgr.variables['angle'][:] # angle of xi (due east)
  lon_rho = fgr.variables['lon_rho'][:]
  lat_rho = fgr.variables['lat_rho'][:]
  h =  fgr.variables['h'][:]
  SY, SX = mask.shape

  #########################################################
  # Plot data
  #########################################################

  # IMPORTANT: Shift xi,eta vectors to EW,NS vectors
  map_EW = xi_component*np.cos(angle) - eta_component*np.sin(angle)
  map_NS = xi_component*np.sin(angle) + eta_component*np.cos(angle)

  # Prepare map projection
  latmin = np.min(lat_rho)-0.1
  latmax = np.max(lat_rho)+0.1
  lonmin = np.min(lon_rho)-0.1
  lonmax = np.max(lon_rho)+0.1
  latmed = (latmin+latmax)/2 # latitude of true scale

  # Create Mercator map:
  m = Basemap(projection='merc',lat_ts=latmed,\
          llcrnrlat=latmin,urcrnrlat=latmax,\
          llcrnrlon=lonmin,urcrnrlon=lonmax,\
          resolution='c')

  x, y = m(lon_rho, lat_rho)

  # Compute ticks on x and y axis:
  xticks = np.arange(110.,310.,10.)
  xtick_labels_W = [ '{:g}W'.format(t) for t in np.arange(170.,55.,-10.)]
  xtick_labels_E = [ '{:g}E'.format(t) for t in np.arange(110.,190.,10.)]
  xtick_labels = xtick_labels_E + xtick_labels_W

  yticks = np.arange(-40.,70.,10.)
  ytick_labels_N = [ '{:g}N'.format(t) for t in np.arange(0.,75.,10.)]
  ytick_labels_S = [ '{:g}S'.format(t) for t in np.arange(40.,0.,-10.)]
  ytick_labels = ytick_labels_S + ytick_labels_N
  # Convert ticks to map coordinates:
  xticks_map, _ = m( xticks,latmin*np.ones(len(xticks)) )
  _, yticks_map = m( lonmin*np.ones(len(yticks)),yticks )

  # Initiate Plot
  fig = plt.figure(facecolor='w')
  ax = fig.add_axes([0.05,0.05,0.9,0.9])

  # Plot some underlying data
  im = m.pcolormesh(x,y,h,shading='flat',cmap='RdBu_r')


  # Set vector field properties
  every = 10 # Plot every-th point on grid
  scaler = 2.# Scale length of vectors (scaler large = vector small)

  # Plot vector field
  mq = m.quiver(x[::every,::every], y[::every,::every],
          map_EW[::every,::every], 
          map_NS[::every,::every],
          scale=scaler, color='k', alpha=0.7,width=0.001)

  # Define vector key ('Scale for vector length')
  qk = plt.quiverkey(mq, 0.2, 0.05, 1, '1 N/m2', labelpos='W', color='k')

  m.drawcoastlines()
  m.drawcountries()
  m.fillcontinents()
  m.drawparallels(np.arange(-50,70.,10.),dashes=[1,4])
  m.drawmeridians(np.arange(100.,310.,10.),dashes=[1,4])

  plt.xticks(xticks_map,xtick_labels, fontsize=12)
  plt.yticks(yticks_map, ytick_labels, fontsize=12)

  plt.title('vector field pacsg', fontsize=20)

  cbar = plt.colorbar()
  cbar.ax.set_ylabel('Bottom Topography (m)', fontsize=16)
  cl = plt.getp(cbar.ax, 'ymajorticklabels')
  plt.setp(cl, fontsize=12)

  plt.show(block=False)
    

      
      
def get_vel_cmap():
    import numpy as np
    from matplotlib.colors import LinearSegmentedColormap

    palette = [np.linspace(0,1,100), np.linspace(0.1,0.85,100),np.linspace(0.3,0.85,100)]

    b3 = palette[2] # value of blue at sample n
    b2 = palette[2] # value of blue at sample n
    b1 = np.linspace(0,1,len(b2)) # position of sample n - ranges from 0 to 1

    # setting up columns for list
    g3 = palette[1]
    g2 = palette[1]
    g1 = np.linspace(0,1,len(g2))

    r3 = palette[0]
    r2 = palette[0]
    r1 = np.linspace(0,1,len(r2))

    # creating list
    R = zip(r1,r2,r3)
    G = zip(g1,g2,g3)
    B = zip(b1,b2,b3)

    # transposing list
    RGB = zip(R,G,B)
    rgb = zip(*RGB)
    # print rgb

    # creating dictionary
    k = ['red', 'green', 'blue']
    LinearL = dict(zip(k,rgb)) # makes a dictionary from 2 lists
    my_cmap = LinearSegmentedColormap('my_colormap',LinearL)

    return my_cmap


# Helper function for nonlinear cmap
class nlcmap(LinearSegmentedColormap):
    """A nonlinear colormap"""

    name = 'nlcmap'

    def __init__(self, cmap, levels):
        self.cmap = cmap
        self.monochrome = self.cmap.monochrome
        self.levels = np.asarray(levels, dtype='float64')
        self.levmax = self.levels.max()
        self.levmin = self.levels.min()
        self._x = (self.levels - self.levmin) / (self.levmax - self.levmin)
        self._y = np.linspace(0, 1, len(self.levels))

    def __call__(self, xi, alpha=1.0, **kw):
        yi = np.interp(xi, self._x, self._y)
        return self.cmap(yi, alpha)
