#!/usr/bin/env python
import numpy as np
from myromstools import compute_zlev, compute_dz

def bisect_sigma_SingleDepth(romsGrd,depth):
    '''
    Find neighbouring indices around a certain depth (2D array or scalar) supplied 
    to the function. If 2D, depth can be varying from point to point (e.g. MLD).
    
    Used to calculate the mean above a certain depth (see also mean_sigma) or
    values interpolated to certain depth (see also atDepth)

    Inputs
    depth:   either 2D array of shape (eta,xi) or scalar value of constant depth
    romsGrd: ROMS grid class (from myromstools getGrid) including grdfile and
             clmfile attributes
    '''
    
    # Check depth is 2D array or scalar
    if type(depth)==float:
        # Make it a 2D array of ROMS grid dimensions
        depth = np.ones([romsGrd.SY,romsGrd.SX])*depth
    # Make sure depths are positively defined
    depth[depth<0.0] *= -1.
    
    # Compute ROMS depths (rho levels):
    print('Bisecting sigma layers...')
    depth_w = -compute_zlev(romsGrd.ncclm,romsGrd.ncgrd,romsGrd.NZ,type='w',stype=3)
    
    # Find indices above and below specified depth (defined on sigma w-points)
    sortidx = np.argpartition(np.abs(depth_w-depth),2,axis=0)
    i1 = np.minimum(sortidx[0],sortidx[1])
    i2 = np.maximum(sortidx[0],sortidx[1])

    # Load land-sea mask
    mask = np.copy(romsGrd.mask_rho)
    # Apply land-sea mask
    i1 = i1*mask
    i2 = i2*mask
    # If depth is below deepest level (update mask)
    mask[i1==0] = 0
    # If depth is above highest level
    i2[i1==romsGrd.NZ-1] = romsGrd.NZ-1

    # Return
    return i1,i2,depth_w,mask



# Find bisection indices in z direction
def bisect_sigma_MultipleDepths(romsGrd,depths):
    
    import bisect
    
    if depths is None: # Use standard depths
        Z = -np.array([0,5,10,15,20,25,30,35,40,45,50,60,70,80,90,100,125,150,175,200,225,250,275,300,350,400,500])
    else:
        if isinstance(depths,list):
            # Make sure values are negative
            Z = [-1*abs(i) for i in depths]
        else:
            Z = depths
            # Make sure values are negative
            Z[Z>0] *= -1.
            
    print('Bisecting sigma levels to depths: {}'.format(Z))
   
    # Read mask_rho from grid file:
    print('Read masks')
    mask_grd = np.copy(romsGrd.mask_rho)

    # Compute ROMS depths (rho levels):
    print('Compute sigma depths')
    z0 = compute_zlev(romsGrd.ncclm,romsGrd.ncgrd,romsGrd.NZ,type='rho',stype=3)

    # Loop over Z to find nearest sigma levels: only for grid type rho
    # we compute idx1, idx2 (indices of sigma levels enclosing each Z level)
    # and a mask: mask[k,i,j] = 1 if ocean at cell i,j has depth Z[k] or larger.
    print('Compute idx1 and idx2 for rho grid')
    mask = np.zeros([len(Z),romsGrd.SY,romsGrd.SX])
    i1 = np.zeros(mask.shape,np.int)
    i2 = np.zeros(mask.shape,np.int)
    for k in range(len(Z)):
        depth = Z[k]
        mask[k,:,:] = mask_grd
        for j in range(romsGrd.SY):
            for i in range(romsGrd.SX):
                if mask[k,j,i] == 1:
                    tmp = bisect.bisect_right(z0[:,j,i],depth)
                    # z0[:,j,i] is ascending (i.e. k=1 is at bottom) for each
                    # i and j:
                    if tmp == 0:
                        # At this place, the water column is not deep enough so
                        # that no value will be available in the horizontal
                        # section. Therefore we adapt mask:
                        mask[k,j,i] = 0
                    elif tmp == romsGrd.NZ:
                        # Here we are at or above the highest sigma level: take
                        # uppermost value
                        i1[k,j,i] = romsGrd.NZ-1; i2[k,j,i] = romsGrd.NZ-1
                    else:
                        i1[k,j,i] = tmp - 1; i2[k,j,i] = tmp
    # Return indices and mask
    return i1,i2,mask


def compute_weights(romsGrd,depth1,depth2=0.,depth_w=None,depth1_i1=None,depth1_i2=None,depth2_i1=None,depth2_i2=None):

    '''
    Inputs
    depth1: lower bisection depth, scalar
    depth2: upper bisection depth, scalar (default=0.)
    depth_w: depth of sigma level w-points, 3D array
    i1,i2: bisection indices in vertical direction (from bisect_sigma_SingleDepth)    
    '''
    
    print('Computing averaging weights...')
    
    # If depth_w, i1, i2 are None, bisect sigma layers first
    if (depth_w is None) or (depth1_i1 is None) or (depth1_i2 is None):
        depth1_i1,depth1_i2,depth_w,mask = bisect_sigma_SingleDepth(romsGrd,depth1)
    
    if depth2==0.0:
        depth2_i1 = np.ones([romsGrd.SY,romsGrd.SX])*(romsGrd.NZ-1)
        depth2_i2 = np.ones([romsGrd.SY,romsGrd.SX])*(romsGrd.NZ-1)
    else:
        if (depth_w is None) or (depth2_i1 is None) or (depth2_i2 is None):
            depth2_i1,depth2_i2,depth_w,mask = bisect_sigma_SingleDepth(romsGrd,depth2)
                
    if type(depth1)==float:
        depth1 = np.ones([romsGrd.SY,romsGrd.SX])*depth1
    if type(depth2)==float:
        depth2 = np.ones([romsGrd.SY,romsGrd.SX])*depth2

    # Make sure depth is positive
    depth1[depth1<0.0] *= -1
    depth2[depth2<0.0] *= -1

    # dz sigma layers (between sigma w-points)
    weights = depth_w[:-1]-depth_w[1:]
   
    # Adjust mask if chosen depth is below topography
    mask = np.copy(romsGrd.mask_rho)
    
    # Apply mask
    weights = weights*mask
    
    # Compute weights    
    for i in xrange(weights.shape[1]):
        for j in xrange(weights.shape[2]):
            if mask[i,j]==1:
                # Only interpolate weights where topography below depth1
                # Otherwise, the whole water columnn is considered
                if depth1[0,0]<depth_w[0,i,j]:
                    for k in range(romsGrd.NZ):
                        if k<depth1_i1[i,j]:
                            weights[k,i,j] = 0.0
                        elif k==depth1_i1[i,j]:
                            weights[k,i,j] = weights[k,i,j]-(depth_w[k,i,j]-depth1[i,j])
                        elif (k>depth1_i2[i,j]) and (k<depth2_i1[i,j]):
                            weights[k,i,j] = weights[k,i,j]
                        elif k==depth2_i1[i,j]:
                            weights[k,i,j] = (depth_w[k,i,j]-depth2[i,j])
                        elif k>depth2_i1[i,j]:
                            weights[k,i,j] = 0.0
    
    # Make sure values are all positive
    weights[weights<0.0] = 0.0
    
    # Return
    return weights,mask



def atDepth(romsGrd,data_dict,depth,idx1=None,idx2=None,masks=None):
    '''
    Some fancy function help text    
    '''
    # Make sure depth is negatively defined in this case
    if depth>0.0:
        depth*=-1
    # Transform input if necessary
    returnKey = False
    if not isinstance(data_dict,dict):
        returnKey = True
	data_dict = {'data':data_dict}
    
    data_atDepth = intROMS_zlevs(romsGrd,data_dict,depths=depth,idx1=idx1,idx2=idx2,masks=masks)

    if returnKey:
	return data_atDepth['data']
    else:
        return data_atDepth

 

def integrate_sigma(romsGrd,data,depth1=None,depth2=None,weights=None,mask=None):
    '''
    Integrate data on a masked array on sigma layers.
    
    data: 3D or 4D (masked) array (e.g. masked depths>100m)
    dz_sigma: dz between sigma layers (computed from compute_zlev.compute_dz)
    mask: optional mask_rho to mask returned array
    
    martinfr, Nov 2015
    '''
    
    data_mean = mean_sigma(romsGrd,data,depth1,depth2,weights)
    data_integrated = np.zeros(data_mean.shape)
    dz = depth1-depth2
    
    if data_mean.ndim==3: # Time dimension present in data
        data_integrated = np.zeros((data_mean.shape[0],)+data_mean.shape[1:])
    elif data_mean.ndim==2: # Time mean field already
        data_integrated = np.zeros((1,)+data_mean.shape[:])
        data_mean = np.ma.expand_dims(data_mean,axis=0)
        
    for t in range(data_mean.shape[0]):
        data_integrated[t,:,:] = dz*data_mean[t]
                
    # Prepare Output
    data_integrated = np.ma.squeeze(data_integrated)
    if mask is not None:
	if data_integrated.ndim==3:
	    n = data_integrated.shape[0]
            data_integrated = np.ma.masked_where(np.array([mask==0]*n),data_integrated) 
        else:
            data_integrated = np.ma.masked_where(mask==0,data_integrated) 
        
    return data_integrated


def integrate_sigma_WaterColumn(romsGrd,data,weights=None,mask=None):
    '''
    Integrate data on a masked array on sigma layers.
    
    data: 3D or 4D (masked) array (e.g. masked depths>100m)
    dz_sigma: dz between sigma layers (computed from compute_zlev.compute_dz)
    mask: optional mask_rho to mask returned array
    
    martinfr, Nov 2015
    '''
       
    if weights is None:
        weights = compute_dz(romsGrd.ncclm,romsGrd.ncgrd,romsGrd.NZ,stype=3)
    
    if data.ndim==4: # Time dimension present in data
        data_integrated = np.zeros((data.shape[0],)+data.shape[2:])
    elif data.ndim==3: # Time mean field already
        data_integrated = np.zeros((1,)+data.shape[1:])
        data = np.ma.expand_dims(data,axis=0)

    for t in range(data.shape[0]):
        data_integrated[t,:,:] = np.ma.sum(weights*data[t],axis=0)
                
    # Prepare Output
    data_integrated = np.ma.squeeze(data_integrated)
    if mask is not None:
        data_integrated = np.ma.masked_where(mask==0,data_integrated)
    else:
        try:
            data_integrated = np.ma.masked_where(romsGrd.mask_rho==0,data_integrated)
        except:
            pass
        
    return data_integrated

  
def mean_sigma(romsGrd,data,depth1=100.,depth2=0.0,weights=None,mask=None,verbose=False):
    '''
    Computes the mean over a certain depth on sigma layers using weights computed
    in compute_weiths(). If e.g. the mean is computed several times in a loop or so,
    make sure to pass the weights (and possibly the mask) to speed up calcualtions.
    
    data:    3D array of certain property to be averaged
    depth:   either 2D array of spatially varying depths (e.g. MLD) or scalar 
             of one specific depth.
    weights: computed from compute_weights()
    mask:    masked areas depending on depth of averaging (output from compute_weights())
    verbose: True/False for additional feedback from function
           
    martinfr, Nov 2015
    '''
    
    if weights is None:
        weights,mask = compute_weights(romsGrd,depth1,depth2)   
    
    if data.ndim==4: # Time dimension present in data
        data_mean = np.ma.empty((data.shape[0],)+data.shape[2:])
    elif data.ndim==3: # Time mean field already
        data_mean = np.ma.empty((1,)+data.shape[1:])
        data = np.ma.expand_dims(data,axis=0)

    for t in range(data.shape[0]):
        if verbose:
            print('Computing vertical mean for t {}...'.format(t) )
        data_mean[t,:,:] = np.ma.average(data[t],axis=0,weights=weights)  
        # Mask data
	try:
            data_mean[t] = np.ma.masked_where(msk==0,data_mean[t])
        except:
            pass
     
    # Prepare Output
    data_mean = np.ma.squeeze(data_mean)
        
    return data_mean


# Interpolate ROMS output to fixed depths.
def intROMS_zlevs(romsGrd,data_dict,type='rho',depths=None,idx1=None,idx2=None,masks=None):
    '''
    Some fancy function help text
    '''

    # Comput ROMS depths (rho levels):
    z0 = compute_zlev(romsGrd.ncclm,romsGrd.ncgrd,romsGrd.NZ,type='rho',stype=3)

    # Set depth array
    if depths is None: # Use standard depths
        Z = -np.array([0,5,10,15,20,25,30,35,40,45,50,60,70,80,90,100,125,150,175,200,225,250,275,300,350,400,500])
    else:
        if isinstance(depths,list):
            Z = depths
            Z = [-1*i for i in Z if i>0.]
        else:
            if isinstance(depths,int):
		Z = [depths]
	    elif isinstance(depths,float):
	        Z = [depths]
	    else:
		Z = depths
            if len(Z)>1:
                Z = [-1*i for i in Z if i>0.]
            elif len(Z)==1:
                if Z[0]>0.:
                    Z[0] *= -1.
                
    # If no bisection indices are passed as arguments, bisect vertical grid
    if (idx1 is None) or (idx2 is None) or (masks is None):
        if  len(Z) > 1:
            idx1,idx2,masks = bisect_sigma_MultipleDepths(romsGrd,depths)
        elif len(Z)==1:
            idx1,idx2,depth_w,masks = bisect_sigma_SingleDepth(romsGrd,depths)
    
    # Create output file
    data_out = dict()
    # Loop over variables in clim file and interpolate to Z depths:
    for var in data_dict.keys():
        buf = data_dict[var][:]
        # Assign type specific entities:
        if  len(Z)==1:
            i1 = np.expand_dims(idx1, axis=0)
            i2 = np.expand_dims(idx2, axis=0)
            mask = np.expand_dims(masks, axis=0)        
        else:
            i1 = idx1
            i2 = idx2
            mask = masks
        # Prepare output buffer:
        if buf.ndim==4:
            data_out[var] = np.ma.empty([buf.shape[0],len(Z),romsGrd.SY,romsGrd.SX])
        else:
            data_out[var] = np.ma.empty([len(Z),romsGrd.SY,romsGrd.SX])
        print('Interpolating {}'.format(var))
        for k in range(len(Z)):
            depth = Z[k]
            for i in range(romsGrd.SX):
                for j in range(romsGrd.SY):
                    if mask[k,j,i] == 1:
                        if i1[k,j,i] == i2[k,j,i]:
                            if buf.ndim==4:
                                data_out[var][:,k,j,i] = buf[:,i1[k,j,i],j,i]
                            else:
                                data_out[var][k,j,i] = buf[i1[k,j,i],j,i]
                        else:
                            # Here we do the linear interpolation:
                            z1 = z0[i1[k,j,i],j,i]
                            z2 = z0[i2[k,j,i],j,i]
                            if buf.ndim==4:
                                t1 = buf[:,i1[k,j,i],j,i]
                                t2 = buf[:,i2[k,j,i],j,i]
                                data_out[var][:,k,j,i] = (depth-z2) / (z1-z2) * (t1-t2) + t2
                            else:
                                t1 = buf[i1[k,j,i],j,i]
                                t2 = buf[i2[k,j,i],j,i]
                                data_out[var][k,j,i] = (depth-z2) / (z1-z2) * (t1-t2) + t2
                    else:
                        if buf.ndim==4:
                            data_out[var][:,k,j,i] = np.NaN
                        else:
                            data_out[var][k,j,i] = np.NaN
        # Squeeze array if len(Z)==1
        if len(Z)==1:
            data_out[var] = np.squeeze(data_out[var])
            data_out[var] = np.ma.masked_where(masks==0,data_out[var])
        # Mask invalid
	data_out[var] = np.ma.masked_invalid(data_out[var])

    # Return interpolated fields
    return data_out
