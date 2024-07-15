
#%% Import Stuff
import numpy as np
import pandas as pd
from numpy import savetxt
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import matplotlib.ticker as ticker
import matplotlib.pylab as pl
import tqdm as tqdm

#%% Define Some Variables
T = 300
kBT=T*8.314462618*0.001

font = "Helvetica"
hfont = {'fontname':'Helvetica'}
big = 22
medium = 18
small = 14
xsmall = 10

title="47arg"
runs=[1,2,3]

#%% Open File Functions

def open_file(path):
    data = pd.read_table("{}".format(path), skiprows=1, sep=" ", names=['rg', 'ffrg', 'stdv'])
    data = data.set_index(np.arange(0, len(data.rg), 1))
    return data

def open_ie(path):
    data = pd.read_table("{}".format(path), sep=" ", skipinitialspace=True, skiprows=27, names=['time', 'csr', 'ljsr', 'c14', 'lj14'])
    return data

def get_rg(rg_path):
    file = open('{}'.format(rg_path))
    lines = file.readlines()
    
    del lines[0:27]
    n = len(lines)
    t_list = []
    v_list = []
    
    for j in range(n):
        [time, value, skip, skip, skip] = lines[j].split()
        t_list.append(time)
        v_list.append(value)
    
    time = np.array(t_list, dtype=float)
    time = time/1000 #Change to ns
    value = np.array(v_list, dtype=float)
    
    rg = np.vstack((value)).T
    return rg

def set_min_barrier(ffrg):
    global get_diff, zero
    get_diff = (ffrg.rg - 0.4).abs().argsort()
    zero = get_diff[0]
    # print('Found', zero)
    min = ffrg.ffrg[zero]
    ffrg.ffrg = ffrg.ffrg - min
    ffrg.ffrg = ffrg.ffrg
    return ffrg

def bin_data(pmf, ie):
    centers=pmf.rg
    step=pmf.rg[1] - pmf.rg[0]
    all_bins = []
    for i in range(len(centers)):
        lower_wall = centers[i]-step
        upper_wall = centers[i]+step
        current_bin = ie[(ie['rg'] >= lower_wall) & (ie['rg'] < upper_wall)]
        all_bins.append(current_bin)

    pwie_means = []
    pwie_err = []
    
    paie_means = []
    paie_err = []
    
    pcie_means = []
    pcie_err = []
    
    for i in range(len(all_bins)):
        current_bin = all_bins[i]
        pwie_means.append(np.mean(current_bin.pwie))
        pwie_err.append(np.std(current_bin.pwie))
        
        paie_means.append(np.mean(current_bin.paie))
        paie_err.append(np.std(current_bin.paie))
        
        pcie_means.append(np.mean(current_bin.pcie))
        pcie_err.append(np.std(current_bin.pcie))
    
    output1 = pd.DataFrame(np.vstack((centers, pwie_means, pwie_err)).T, columns=['rg', 'ffrg', 'stdv'])
    output2 = pd.DataFrame(np.vstack((centers, paie_means, paie_err)).T, columns=['rg', 'ffrg', 'stdv'])
    output3 = pd.DataFrame(np.vstack((centers, pcie_means, pcie_err)).T, columns=['rg', 'ffrg', 'stdv'])

    return output1, output2, output3

def bin_data_vac(pmf, ie):
    centers=pmf.rg
    step=pmf.rg[1] - pmf.rg[0]
    all_bins = []
    for i in range(len(centers)):
        lower_wall = centers[i]-step
        upper_wall = centers[i]+step
        current_bin = ie[(ie['rg'] >= lower_wall) & (ie['rg'] < upper_wall)]
        all_bins.append(current_bin)

    means = []
    err = []
    for i in range(len(all_bins)):
        current_bin = all_bins[i]
        means.append(np.mean(current_bin.ffrg))
        err.append(np.std(current_bin.stdv))
    output = pd.DataFrame(np.vstack((centers, means, err)).T, columns=['rg', 'ffrg', 'stdv'])

    return output

#%% Analysis!

output_list = []
for run in runs:
    print("Run {}".format(run))
    
    # Define the path
    arg_path="C:/Research/.Projects/Virus-Stabilization/Hydrophobic-Polymer/REUS/arginine/{}/check-time/avg-pmfs/average_pmf_100ns.dat".format(title)
    vac_path="C:/Research/.Projects/Virus-Stabilization/Hydrophobic-Polymer/REUS/arginine/bootstrapping/0arg-vac/bootstrapped_average_pmf_1.dat"
    pol_pol_path="C:/Research/.Projects/Virus-Stabilization/Hydrophobic-Polymer/REUS/arginine/{}/interactions/run{}/{}-reus-pol-pol-ie.xvg".format(title, run, title)
    pol_wat_path="C:/Research/.Projects/Virus-Stabilization/Hydrophobic-Polymer/REUS/arginine/{}/interactions/run{}/{}-reus-pol-wat-ie.xvg".format(title, run, title)
    pol_arg_path="C:/Research/.Projects/Virus-Stabilization/Hydrophobic-Polymer/REUS/arginine/{}/interactions/run{}/{}-reus-pol-arg-ie.xvg".format(title, run, title)
    pol_cl_path="C:/Research/.Projects/Virus-Stabilization/Hydrophobic-Polymer/REUS/arginine/{}/interactions/run{}/{}-reus-pol-cl-ie.xvg".format(title, run, title)
    rg_path="C:/Research/.Projects/Virus-Stabilization/Hydrophobic-Polymer/REUS/arginine/{}/gyrate-data/gyrate-{}-reus-{}.xvg".format(title, title, run)

    # Get the data
    pmf = open_file(arg_path)
    vac = open_file(vac_path)
    pol_pol_ie = open_ie(pol_pol_path)
    pol_wat_ie = open_ie(pol_wat_path)
    pol_arg_ie = open_ie(pol_arg_path)
    pol_cl_ie = open_ie(pol_cl_path)
    rg = np.squeeze(get_rg(rg_path))
    
    pmf_err = pmf.stdv
    vac_err = vac.stdv
    
    # Convert to kT units
    ppie = (pol_pol_ie.ljsr + pol_pol_ie.lj14 / kBT)
    pwie = (pol_wat_ie.ljsr / kBT)
    paie = (pol_arg_ie.ljsr / kBT)
    pcie = (pol_cl_ie.ljsr / kBT)
    
    # Stack interaction energy terms
    ie = pd.DataFrame(np.vstack((rg, ppie, pwie, paie, pcie)).T, columns=['rg', 'ppie', 'pwie', 'paie', 'pcie'])
    
    # Bin interaction energy terms and vacuum terms to match original PMF
    pwie, paie, pcie = bin_data(pmf, ie)
    vac = bin_data_vac(pmf, vac)
    
    # Normalize everything!
    pmf = set_min_barrier(pmf)
    vac = set_min_barrier(vac)
    pwie, paie, pcie = set_min_barrier(pwie), set_min_barrier(paie), set_min_barrier(pcie)
    
    # Compute cavitation contribution
    cav = pmf.ffrg - vac.ffrg - pwie.ffrg - paie.ffrg - pcie.ffrg
    cav_err = np.sqrt(pmf.stdv**2 + vac.stdv**2 + pwie.stdv**2 + paie.stdv**2 + pcie.stdv**2)
    cav = pd.DataFrame(np.vstack((pmf.rg, cav, cav_err)).T, columns=['rg', 'ffrg', 'stdv'])

    output_list.append((pmf, vac, cav, pwie, paie, pcie))
    
    xlabels = ["$R_{g}$ [nm]"]
    ylabels = [r"$W(R_{g})$ [kT]", r"$W_{vac}(R_{g})$ [kT]", 
                r"$W_{cav}(R_{g})$ [kT]", r"$\langle$$U_{pw}(R_{g})$$\rangle$ [kT]",
                r"$\langle$$U_{pa}(R_{g})$$\rangle$ [kT]", r"$\langle$$U_{pc}(R_{g})$$\rangle$ [kT]"]

    xlims = [(0.35, 0.9)]
    ylims = [(-1, 3), (-1,3), (-7,17), (-20,10), (-5,5), (-1,1)]

    x = [pmf.rg]
    y = [pmf.ffrg, vac.ffrg, cav.ffrg, pwie.ffrg, paie.ffrg, pcie.ffrg]
    err = [pmf.stdv, vac.stdv, cav.stdv, pwie.stdv, paie.stdv, pcie.stdv]

    fig, ax = plt.subplots(2,3, figsize=(18,10))
    plt.subplots_adjust(wspace=0.3, hspace=0.3)

    color="black"
    m=0 # Row
    k=0 # Column
    for i in range(len(ylabels)):
        ax[m,k].errorbar(x[0], y[i], err[i], color=color, alpha=0.8, ls='-', lw=2, elinewidth=0.1, capsize=0.1)
        # plt.fill_between(x[i], (y[i] + err[i]), (y[i] - err[i]), alpha=0.1, lw=0, color=color)
        ax[m,k].set_xlabel(xlabel=("{}".format(xlabels[0])), fontsize=big, **hfont)
        ax[m,k].set_ylabel(ylabel=("{}".format(ylabels[i])), fontsize=big, **hfont)
        ax[m,k].tick_params(axis='both', which='both', direction='in', labelsize=medium)
        ax[m,k].xaxis.set_major_locator(ticker.MultipleLocator(0.1))
        ax[m,k].xaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax[m,k].yaxis.set_major_locator(ticker.AutoLocator())
        ax[m,k].yaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax[m,k].set_xlim(xlims[0])
        ax[m,k].set_ylim(ylims[i])
        ax[m,k].set_ylabel(ylabel=("{}".format(ylabels[i])), fontsize=big, **hfont)
        # ax[m,k].tick_params(axis='both', which='major', length=12, width=2)
        # ax[m,k].tick_params(axis='both', which='minor', length=8, width=2)
        labels=ax[m,k].get_xticklabels() + ax[m,k].get_yticklabels()
        [label.set_fontname(font) for label in labels]
        
        k+=1
        if k == 3:
            m+=1
            k=0


# Get averages from replicate runs
means = []
errs = []
for i in range(len(output_list[0])):
    curr_ffrg = []
    for j in range(len(output_list)):
        ffrg = output_list[j][i].ffrg
        curr_ffrg.append(ffrg)
    means.append(np.mean(curr_ffrg, axis=0))
    errs.append(np.std(curr_ffrg, axis=0))

pmf = pd.DataFrame(np.vstack((pmf.rg, means[0], pmf_err)).T, columns=['rg', 'ffrg', 'stdv'])
vac = pd.DataFrame(np.vstack((pmf.rg, means[1], vac_err)).T, columns=['rg', 'ffrg', 'stdv'])   
cav = pd.DataFrame(np.vstack((pmf.rg, means[2], errs[2])).T, columns=['rg', 'ffrg', 'stdv'])        
pwie = pd.DataFrame(np.vstack((pmf.rg, means[3], errs[3])).T, columns=['rg', 'ffrg', 'stdv'])        
paie = pd.DataFrame(np.vstack((pmf.rg, means[4], errs[4])).T, columns=['rg', 'ffrg', 'stdv'])
pcie = pd.DataFrame(np.vstack((pmf.rg, means[5], errs[5])).T, columns=['rg', 'ffrg', 'stdv'])

header = "rg val err"
savetxt('pmf-{}.dat'.format(title), pmf, header=header, delimiter=' ', fmt='%s') 
savetxt('vac-{}.dat'.format(title), vac, header=header, delimiter=' ', fmt='%s')
savetxt('cav-{}.dat'.format(title), cav, header=header, delimiter=' ', fmt='%s') 
savetxt('pwie-{}.dat'.format(title), pwie, header=header, delimiter=' ', fmt='%s') 
savetxt('paie-{}.dat'.format(title), paie, header=header, delimiter=' ', fmt='%s') 
savetxt('pcie-{}.dat'.format(title), pcie, header=header, delimiter=' ', fmt='%s')

#%%Plot
xlabels = ["$R_{g}$ [nm]"]
ylabels = [r"$W(R_{g})$ [kT]", r"$W_{vac}(R_{g})$ [kT]", 
            r"$W_{cav}(R_{g})$ [kT]", r"$\langle$$U_{pw}(R_{g})$$\rangle$ [kT]",
            r"$\langle$$U_{pa}(R_{g})$$\rangle$ [kT]", r"$\langle$$U_{pc}(R_{g})$$\rangle$ [kT]"]

xlims = [(0.35, 0.9)]
ylims = [(-1, 3), (-1,3), (-7,17), (-20,10), (-5,5), (-1,1)]

x = [pmf.rg]
y = [pmf.ffrg, vac.ffrg, cav.ffrg, pwie.ffrg, paie.ffrg, pcie.ffrg]
err = [pmf.stdv, vac.stdv, cav.stdv, pwie.stdv, paie.stdv, pcie.stdv]

fig, ax = plt.subplots(2,3, figsize=(18,10))
plt.subplots_adjust(wspace=0.3, hspace=0.3)

color="black"
m=0 # Row
k=0 # Column
for i in range(len(ylabels)):
    ax[m,k].errorbar(x[0], y[i], err[i], color=color, alpha=0.8, ls='-', lw=2, elinewidth=0.1, capsize=0.1)
    # plt.fill_between(x[i], (y[i] + err[i]), (y[i] - err[i]), alpha=0.1, lw=0, color=color)
    ax[m,k].set_xlabel(xlabel=("{}".format(xlabels[0])), fontsize=big, **hfont)
    ax[m,k].set_ylabel(ylabel=("{}".format(ylabels[i])), fontsize=big, **hfont)
    ax[m,k].tick_params(axis='both', which='both', direction='in', labelsize=medium)
    ax[m,k].xaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ax[m,k].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax[m,k].yaxis.set_major_locator(ticker.AutoLocator())
    ax[m,k].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax[m,k].set_xlim(xlims[0])
    ax[m,k].set_ylim(ylims[i])
    ax[m,k].set_ylabel(ylabel=("{}".format(ylabels[i])), fontsize=big, **hfont)
    # ax[m,k].tick_params(axis='both', which='major', length=12, width=2)
    # ax[m,k].tick_params(axis='both', which='minor', length=8, width=2)
    labels=ax[m,k].get_xticklabels() + ax[m,k].get_yticklabels()
    [label.set_fontname(font) for label in labels]
    
    k+=1
    if k == 3:
        m+=1
        k=0

plt.savefig('PMF-decomposition.svg', transparent=True)
plt.show()

# test = cav.ffrg + vac.ffrg + pwie.ffrg + paie.ffrg + pcie.ffrg
# plt.plot(pmf.rg, pmf.ffrg, ls='-', lw=2)
# plt.plot(pmf.rg, test, ls='--', lw=2)    

