#%% Import Some Stuff
import sys
sys.path.append('C:/Research/.Projects/Universal-Files')
import blockaveragingstdv
import numpy as np
import MDAnalysis as md
from MDAnalysis.analysis import contacts
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.ticker as ticker
from numpy import savetxt

#%% User Inputs

k=47
z="arg"
z1="arg"
title = "{}{}".format(k, z1)
path="C:/Research/.Projects/Virus-Stabilization/Hydrophobic-Polymer/REUS/arginine/{}".format(title)
pmf_path = "C:/Research/.Projects/Virus-Stabilization/Hydrophobic-Polymer/REUS/arginine/avg-pmfs"
opath="C:/Research/.Projects/Virus-Stabilization/Hydrophobic-Polymer/REUS/arginine/analysis/solvation-new/output"

runs = [1, 2, 3]

# Figure Parameters
font = "Helvetica"
hfont = {'fontname':'Helvetica'}
big = 22
medium = 18
small = 14
xsmall = 10
alpha=0.2

#%% Extract Simulation Data
def open_files():
    solv = pd.read_table("{}/pref-int-data/local-{}-rcut0.55-{}-reus-{}.xvg".format(path, z, title, run), index_col=False, skipinitialspace=True, skiprows=24, sep=" ", names=['time','solv'])
    solv = solv.set_index(np.arange(0, len(solv.solv), 1))
    
    rg = pd.read_table("{}/gyrate-data/gyrate-{}-reus-{}.xvg".format(path, title, run), index_col=False, skipinitialspace=True, skiprows=27, sep=" ", names=['time','rg'])
    rg = rg.set_index(np.arange(0, len(rg.rg), 1))
    return solv, rg

def get_data():
    data = pd.read_table("{}/average_pmf_{}.dat".format(pmf_path, title), skiprows=1, sep=" ", names=['rg', 'ffrg', 'stdv'])
    data = data.set_index(np.arange(0, len(data.rg), 1))
    data.replace([np.inf, -np.inf], np.nan, inplace=True) #drop infinite values
    data.dropna(subset=["ffrg"], how="all", inplace=True)
    pdf = np.exp(-1 * data.ffrg) / np.sum(np.exp(-1 * data.ffrg))
    pdf = pd.concat((data.rg, pdf), axis=1)
    pdf = pdf.set_index(np.arange(0, len(pdf.rg), 1))
    
    return pdf, data.ffrg

#%% Plotting functions
def plot(x, y, err):
    fig, ax = plt.subplots()

    plt.plot(x, y, color="black", alpha=1, lw=1)
    plt.fill_between(x, (y - err), (y + err), color="black", alpha=0.4, lw=0)

    plt.autoscale()
    plt.xlim(0.35, 0.9)
    # plt.ylim(0.85, 1.05)
    ax.xaxis.set_major_locator(ticker.AutoLocator())
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    # ax.yaxis.set_major_locator(ticker.MultipleLocator(0.05))
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    plt.xticks(fontsize=small, **hfont)
    plt.yticks(fontsize=small, **hfont)
    plt.xlabel("R$_g$ [nm]", **hfont, fontsize=big)
    plt.ylabel(r"Y", **hfont, fontsize=big)

    plt.tight_layout()
    # plt.savefig('{}/desolv-{}-norm-ens-avg.svg'.format(opath, title), transparent=True)
    plt.show()
    
#%% Unbiasing Function
walls=[]
centers=[]
offset=0.01
n_centers = 12
iters = n_centers * (0.05/offset) + 1
for i in range(int(iters)):
    val = round(((0.35-(offset/2)) + (i*offset)), 3)
    cen = round((val+(offset/2)), 3)
    walls.append(val)
    centers.append(cen)
centers.remove(centers[-1])

def bin_data():
    output = np.vstack((rg.rg, solv.solv)).T
    # output = np.vstack((rg.rg, fluct)).T
    output = pd.DataFrame(output[output[:,0].argsort()], columns=["rg", "s"])
    
    i=0
    all_bins = []
    for i in range(len(walls)-1):
        lower_wall = walls[i]
        upper_wall = walls[i+1]
        current_bin = output[(output['rg'] >= lower_wall) & (output['rg'] < upper_wall)]
        all_bins.append(current_bin)
        
    bin_means = []
    bin_err = []
    for i in range(len(all_bins)):
        current_bin = all_bins[i]
        bin_means.append(np.mean(current_bin.s))
        bin_err.append(np.std(current_bin.s))
    
    bins = np.column_stack((centers, bin_means, bin_err))
    bins = pd.DataFrame(bins, columns=["rg", "s", "err"])
    
    return bins, all_bins

def unbias_data(output):
    i=0
    all_bins = []
    mean_pdfs = []
    mean_rgs = []
    for i in range(len(pdf)-1):
        lower_wall = pdf.rg[i]
        upper_wall = pdf.rg[i+1]
        current_bin = output[(output['rg'] >= lower_wall) & (output['rg'] < upper_wall)]
        mean_rg = np.mean(pdf.rg[i] + pdf.rg[i+1]) / 2
        mean_pdf = np.mean(pdf.ffrg[i] + pdf.ffrg[i+1]) / 2
        all_bins.append(current_bin)
        mean_rgs.append(mean_rg)
        mean_pdfs.append(mean_pdf)
    
    bin_means = []
    bin_err = []
    for i in range(len(all_bins)):
        current_bin = all_bins[i]
        bin_means.append(np.mean(current_bin.s))
        bin_err.append(np.std(current_bin.s))
    
    bins = np.column_stack((mean_rgs, bin_means, bin_err, mean_pdfs))
    
    return bins

def unbias_nw():
    output = np.vstack((rg.rg, solv.solv)).T
    output = pd.DataFrame(output[output[:,0].argsort()], columns=["rg", "s"])
    bins = unbias_data(output)
    bins = pd.DataFrame(bins, columns=["rg", "s", "err", "pdf"])
    u_solv = np.sum(bins.s * bins.pdf) / np.sum(bins.pdf)
    u_solv_sq = np.sum((bins.s**2) * bins.pdf) / np.sum(bins.pdf)
    
    return u_solv, u_solv_sq

def set_min(s):
    get_diff = (bins.rg - 0.8).abs().argsort()
    zero = get_diff[0]
    val = s[zero]
    s = s/val
    return s

#%%
bins_list = []
all_bins_list = []
all_fluct = []
all_s = []

for run in runs:
    pdf, ffrg = get_data()
    solv, rg = open_files()
    
    # Put into bins along RC
    bins, all_bins = bin_data()
    
    fluct_list = []
    for window in all_bins:
        s = window.s
        sq = window.s**2
        s = s.mean(axis=0)
        sq = sq.mean(axis=0)
        sig_nw = (sq) - (s ** 2)
        fluct = sig_nw / s
        
        fluct_list.append(fluct)
    
    # # Option to normalize by ensemble average
    # u_solv, u_solv_sq = unbias_nw()
    # bins.s = bins.s / u_solv
    
    # # Option to normalize to elongated state
    bins.s = set_min(bins.s)

    # Append to lists
    all_fluct.append(fluct_list)
    all_s.append(bins.s)
    bins_list.append(bins)
    all_bins_list.append(all_bins)

fluct_list = np.vstack((all_fluct)).T
fluct_err = np.std(fluct_list, axis=1)
fluct = np.mean(fluct_list, axis=1)

s_list = np.vstack((all_s)).T
s_err = np.std(s_list, axis=1)
s = np.mean(s_list, axis=1)

plot(bins.rg, fluct, fluct_err)
plot(bins.rg, s, s_err)

# output1 = np.vstack((bins.rg, s, s_err)).T
# header = "Rg Nw Err"
# savetxt('{}/desolv-{}-rc-ens-avg.dat'.format(opath, title), output1, header=header, delimiter=' ', fmt='%s') 

output1 = np.vstack((bins.rg, s, s_err)).T
header = "Rg Nw Err"
savetxt('{}/desolv-{}-{}-rc.dat'.format(opath, title, z), output1, header=header, delimiter=' ', fmt='%s') 

output2 = np.vstack((bins.rg, fluct, fluct_err)).T
header = "Rg Fluct Err"
savetxt('{}/fluct-{}-{}-rc.dat'.format(opath, title, z), output2, header=header, delimiter=' ', fmt='%s') 
