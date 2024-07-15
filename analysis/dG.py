#%% Import Stuff
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d

#%% Let's get dG!
def get_dG(ffrg, kBT, rmin, rmax):
    rcut=0.65
    ffrg = ffrg.replace([np.inf, -np.inf], np.nan) #drop infinite values
    ffrg = ffrg.dropna(subset=["ffrg"], how="all")
    ffrg = ffrg.set_index(np.arange(0, len(ffrg.rg), 1))
    rmin = np.min(ffrg.rg) #find the minimum rg value
    rmax = np.max(ffrg.rg) #find the maximum rg value
    smooth = gaussian_filter1d(ffrg.ffrg, 1)
    diff1 = np.gradient(smooth, axis=0)
    diff2 = np.gradient(diff1, axis=0)
    inflections = np.where(np.diff(np.sign(diff1)))[0]
    
    infl=np.abs(ffrg.rg[inflections] - rcut).argmin()
    
    rcut = ffrg.rg[inflections[infl]]
    
    get_rcut = (ffrg.rg - rcut).abs().argsort() #find the index associated with cutoff
    rcut_index = get_rcut.iloc[0]
    get_rmax = (ffrg.rg - rmax).abs().argsort() #find the index associated with max
    rmax_index = get_rmax.iloc[0]
    get_rmin = (ffrg.rg - rmin).abs().argsort() #find the index associated with min
    rmin_index = get_rmin.iloc[0]
    
    #Get first dG
    dRg = ffrg.rg[1]-ffrg.rg[0]
    num = np.sum(np.exp(-1 * ffrg.ffrg[rcut_index:rmax_index] / kBT))*dRg #integrate numerator
    den = np.sum(np.exp(-1 * ffrg.ffrg[rmin_index:rcut_index] / kBT))*dRg #integrate denominator
    result = np.log(num / den) * -1 * kBT
    
    wErr = ffrg.stdv # PMF error, stdv from 3 runs
    eErr = np.abs(np.exp(-1 * (ffrg.ffrg)/kBT))*np.abs((1/kBT)*wErr) # exponential error
    iErr = dRg * np.sqrt(np.sum(np.square(eErr))) # integration error
    lnErr = iErr / num # ln error, numerator
    ldErr = iErr / den # ln error, denominator
    
    dgErr = kBT * np.sqrt(lnErr**2 + ldErr**2) # total dG error
                     
    return result, dgErr, rcut
