
# coding: utf-8

# This file is to produce a Taylor diagram of a timeseries at a single location.
# The datasets are thus NOT weighted. 
# The input datasets "reference" and "model" have to be of identical length.
# 
# author: Eike Köhn
# date of file creation: 29.10.2018
# 
# Change Log:
# Date:         Action:
# 29.10.2018    Write fundamental script

# In[7]:


def plot_taylor_diagram(reference,model):
    import matplotlib.pyplot as plt
    import numpy as np
    
    a = reference
    b = model
    Na = np.prod(a.shape)
    Nb = np.prod(b.shape)
    if Na != Nb:
        print('Datasets do not have identical length')
        return
    N = Na
    n = np.arange(N)
    
    # caluculate the statistical characteristics
    a_mean = np.mean(a)
    b_mean = np.mean(b)
    a_std = np.std(a)
    b_std = np.std(b)

    R = 1./N * np.sum((a[i]-a_mean)*(b[i]-b_mean) for i in n) / (a_std*b_std) # Correlation is korrekt
    E_dash = (1./N * np.sum(((b[i]-b_mean)-(a[i]-a_mean))**2 for i in n))**0.5

    # calculate the E_dash contour field
    sigr = a_std
    sigf_dum = np.linspace(0,2*sigr,500)
    Rf_dum = np.linspace(0,1,500)
    sigf,Rf = np.meshgrid(sigf_dum,Rf_dum)
    E_dash_field = (sigf**2+sigr**2-2*sigf*sigr*Rf)**0.5
    
    # set plotting settings 
    xtickval = np.array([1.0,0.99,0.95,0.9,0.8,0.6,0.4,0.2,0])
    xtickang = np.arccos(xtickval)
    xtickstr = [str(xtickval[i]) for i in range(len(xtickval))]

    # Plot the Taylor Diagram
    plt.rcParams['font.size'] = 20
    fig = plt.figure(figsize=(10,10)) 
    ax = plt.subplot(111, projection='polar')
    ax.set_thetamin(0)
    ax.set_thetamax(90)
    ax.set_ylabel('Std. Deviation')
    ax.set_ylim(0,2*a_std)
    c0 = ax.contour(np.arccos(Rf),sigf,E_dash_field,np.array([0.5,1,1.5])*sigr,linestyles='dashed',colors='#999999')
    ax.clabel(c0,fmt='%1.2f')
    ax.scatter(np.arccos(1),a_std,300,marker='*',clip_on=False,label='Reference/Obs.',zorder=3)
    ax.scatter(np.arccos(R),b_std,200,'r',clip_on=False,label='Model',zorder=3)
    ax.set_xticks(xtickang)
    ax.set_xticklabels(xtickstr)
    ax.tick_params(pad=10)
    plt.text(x=0.8,y=2*a_std+0.2,s='Correlation',rotation=-50)
    ax.legend()
    ax.set_title('Taylor Diagram')
    return fig,ax

