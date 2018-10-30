
# coding: utf-8

# This file contains a number of different statistical functions.
# author: Eike Köhn
# date of file creation: 29. Oct. 2018
# 
# change log: 
# Date:         Action:
# 29.10.2018    Added functions m, cov, and corr

# In[2]:


import numpy as np

def m(x, w):
    """Weighted Mean"""
    return np.sum(x * w) / np.sum(w)

def cov(x, y, w):
    """Weighted Covariance"""
    return np.sum(w * (x - m(x, w)) * (y - m(y, w))) / np.sum(w)

def corr(x, y, w):
    """Weighted Correlation"""
    return cov(x, y, w) / np.sqrt(cov(x, x, w) * cov(y, y, w))

