
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
kBT = 1

font = "Helvetica"
hfont = {'fontname':'Helvetica'}
big = 22
medium = 18
small = 14
xsmall = 10
step = 25

#User Inputs
conc = 0.25 # 0.25, 0.5, 1.0, "both"
option = 1 # 1 --> single components, 2 --> mixtures, 3 --> individual species, 4 --> single to mixture
mols = 1 # 1 --> arg, 2 --> glu, 3 --> lys, 4 --> arg/glu, 5 --> arg/lys
manual=True

excluded_titles = ["47gly", "93gly", "185gly"] # Exclusions for pcie

#%% Set Params
if conc == 0.25 and option == 1:
    labels = ["Wat", "Arg", "Glu", "Lys"]
    titles = ["0arg", "47arg", "47glu", "47lys"]
    colors = ["black", "firebrick", "darkcyan", "goldenrod"]

if conc == 0.25 and option == 2:
    labels = ["Wat", "Arg/Arg", "Arg/Glu", "Arg/Lys"]
    titles = ["0arg", "47arg", "48arg-glu", "48arg-lys"]
    colors = ["black", "firebrick", "darkgreen", "chocolate"]

if conc == 0.25 and option == 4:
    labels = ["Wat", "Glu", "Arg/Glu", "Lys", "Arg/Lys"]
    titles = ["0arg", "47glu", "48arg-glu", "47lys", "48arg-lys"]
    colors = ["black", "darkcyan", "darkgreen", "goldenrod", "chocolate"]

if conc == 0.5 and option == 1:
    labels = ["Wat", "Arg", "Glu", "Lys"]
    titles = ["0arg", "93arg", "93glu", "93lys"]
    colors = ["black", "firebrick", "darkcyan", "goldenrod"]

if conc == 1.0 and option == 1:
    labels = ["Wat", "Arg", "Glu", "Lys"]
    titles = ["0arg", "185arg", "185glu", "185lys"]
    colors = ["black", "firebrick", "darkcyan", "goldenrod"]

if conc == 1.0 and option == 2:
    labels = ["Wat", "Arg/Arg", "Arg/Glu", "Arg/Lys"]
    titles = ["0arg", "185arg", "186arg-glu", "186arg-lys"]
    colors = ["black", "firebrick", "darkgreen", "chocolate"]

if conc == 1.0 and option == 4:
    labels = ["Wat", "Glu", "Arg/Glu", "Lys", "Arg/Lys"]
    titles = ["0arg", "185glu", "186arg-glu", "185lys", "186arg-lys"]
    colors = ["black", "darkcyan", "darkgreen", "goldenrod", "chocolate"]
    
if conc == "both" and mols == 1:
    labels = ["0.0 M", "0.25 M", "1.0 M"]
    titles = ["0arg", "47arg", "185arg"]
    colors = ["black", "lightcoral", "firebrick"]

if conc == "both" and mols == 2:
    labels = ["0.0 M", "0.25 M", "1.0 M"]
    titles = ["0arg", "47glu", "185glu"]
    colors = ["black", "paleturquoise", "darkcyan"]

if conc == "both" and mols == 3:
    labels = ["0.0 M", "0.25 M", "1.0 M"]
    titles = ["0arg", "47lys", "185lys"]
    colors = ["black", "khaki", "goldenrod"]

if conc == "both" and mols == 4:
    labels = ["0.0 M", "0.25 M", "1.0 M"]
    titles = ["0arg", "48arg-glu", "186arg-glu"]
    colors = ["black", "mediumseagreen", "darkgreen"]

if conc == "both" and mols == 5:
    labels = ["0.0 M", "0.25 M", "1.0 M"]
    titles = ["0arg", "48arg-lys", "186arg-lys"]
    colors = ["black", "sandybrown", "chocolate"]
    
if conc == 0.25 and option == 4 and mols == 2:
    labels = ["Wat", "Glu", "Arg/Glu"]
    titles = ["0arg", "47glu", "48arg-glu"]
    colors = ["black", "darkcyan", "darkgreen"]

if conc == 1.0 and option == 4 and mols == 2:
    labels = ["Wat", "Glu", "Arg/Glu"]
    titles = ["0arg", "185glu", "186arg-glu"]
    colors = ["black", "darkcyan", "darkgreen"]
    
if conc == 0.25 and option == 4 and mols == 3:
    labels = ["Wat", "Lys", "Arg/Lys"]
    titles = ["0arg", "47lys", "48arg-lys"]
    colors = ["black", "goldenrod", "chocolate"]

if conc == 1.0 and option == 4 and mols == 3:
    labels = ["Wat", "Lys", "Arg/Lys"]
    titles = ["0arg", "185lys", "186arg-lys"]
    colors = ["black", "goldenrod", "chocolate"]
    
if conc == 1.0 and option == 4 and mols == 4:
    labels = ["Wat", "Arg", "Arg/Arg", "Arg/Glu"]
    titles = ["0arg", "93arg", "185arg", "186arg-glu"]
    colors = ["black", "lightcoral", "firebrick", "darkgreen"]

if conc == 1.0 and option == 4 and mols == 5:
    labels = ["Wat", "Arg", "Arg/Arg", "Arg/Lys"]
    titles = ["0arg", "93arg", "185arg", "186arg-lys"]
    colors = ["black", "lightcoral", "firebrick", "goldenrod"]

if conc == 1.0 and option == 4 and mols == 6:
    labels = ["Wat", "Arg", "Arg/Arg", "Arg/Glu", "Arg/Lys"]
    titles = ["0arg", "93arg", "185arg", "186arg-glu", "186arg-lys"]
    colors = ["black", "lightcoral", "firebrick", "darkgreen", "chocolate"]

if conc == 1.0 and option == 4 and mols == 7:
    labels = ["Wat", "Glu", "Glu/Glu" "Arg/Glu"]
    titles = ["0arg", "93glu", "185glu", "186arg-glu"]
    colors = ["black", "paleturquoise", "darkcyan", "darkgreen"]

if conc == 1.0 and option == 4 and mols == 8:
    labels = ["Wat", "Lys", "Lys/Lys" "Arg/Lys"]
    titles = ["0arg", "93lys", "185lys", "186arg-lys"]
    colors = ["black", "khaki", "goldenrod", "chocolate"]

#%% Functions

if manual==True:
    labels = ["Wat", "0.25 M Arg", "0.5 M Arg", "1.0 M Arg"]
    titles = ["0arg", "47arg", "93arg", "185arg"]

    # labels = ["Wat", "0.25 M Arg", "0.25 M Gdm", "0.25 M Gly"]
    # titles = ["0arg", "47arg", "47gdm", "47gly"]
    # colors = ["black", "firebrick", "purple", "darkgreen"]

    # labels = ["Wat", "0.25 M Gdm", "0.5 M Gdm", "1.0 M Gdm"]
    # titles = ["0arg", "47gdm", "93gdm", "185gdm"]
    
    # labels = ["Wat", "0.25 M Gly", "0.5 M Gly", "1.0 M Gly"]
    # titles = ["0arg", "47gly", "93gly", "185gly"]
    
    colors = [pl.cm.Reds(np.linspace(0.3,1,5)), pl.cm.Greys(np.linspace(0.5,0.5,5)), 
              pl.cm.Blues(np.linspace(0.3,1,5)), pl.cm.Greens(np.linspace(0.3,1,5)), 
              pl.cm.YlOrBr(np.linspace(0.3,1,5)), pl.cm.Purples(np.linspace(0.3,1,5))]


def open_file(path):
    pmf = pd.read_table("{}/pmf-{}.dat".format(path, title), skiprows=1, sep=" ", names=['rg', 'ffrg', 'stdv'])
    vac = pd.read_table("{}/vac-{}.dat".format(path, title), skiprows=1, sep=" ", names=['rg', 'ffrg', 'stdv'])
    cav = pd.read_table("{}/cav-{}.dat".format(path, title), skiprows=1, sep=" ", names=['rg', 'ffrg', 'stdv'])
    pwie = pd.read_table("{}/pwie-{}.dat".format(path, title), skiprows=1, sep=" ", names=['rg', 'ffrg', 'stdv'])
    paie = 0
    pcie = 0

    if title != "0arg":
        paie = pd.read_table("{}/paie-{}.dat".format(path, title), skiprows=1, sep=" ", names=['rg', 'ffrg', 'stdv'])
        if title not in excluded_titles:
            pcie = pd.read_table("{}/pcie-{}.dat".format(path, title), skiprows=1, sep=" ", names=['rg', 'ffrg', 'stdv'])
    
    return pmf, vac, cav, pwie, paie, pcie

def plot():
    x = [pmf.rg]
    xlabels = ["$R_{g}$ [nm]"]
    xlims = [(0.35, 0.9)]
    ylims = [(-1, 3), (-1,3), (-7,20), (-20,5), (-7.5,4.5), (-1,1)]
    
    if title not in excluded_titles and title != "0arg":
        y = [pmf.ffrg, vac.ffrg, cav.ffrg, pwie.ffrg, paie.ffrg, pcie.ffrg]
        err = [pmf.stdv, vac.stdv, cav.stdv, pwie.stdv, paie.stdv, pcie.stdv]
        ylabels = [r"$W(R_{g})$ [kT]", r"$W_{vac}(R_{g})$ [kT]", 
                    r"$W_{cav}(R_{g})$ [kT]", r"$E_{pw}(R_{g})$ [kT]",
                    r"$E_{pa}(R_{g})$ [kT]", r"$E_{pc}(R_{g})$ [kT]"]
        alphas=[0.8, 0.0, 0.8, 0.8, 0.8, 0.8]
    
    elif title in excluded_titles:
        y = [pmf.ffrg, vac.ffrg, cav.ffrg, pwie.ffrg, paie.ffrg]
        err = [pmf.stdv, vac.stdv, cav.stdv, pwie.stdv, paie.stdv]
        ylabels = [r"$W(R_{g})$ [kT]", r"$W_{vac}(R_{g})$ [kT]", 
                    r"$W_{cav}(R_{g})$ [kT]", r"$E_{pw}(R_{g})$ [kT]",
                    r"$E_{pa}(R_{g})$ [kT]"]
        alphas=[0.8, 0.0, 0.8, 0.8, 0.8]
        
    if title == "0arg":
        y = [pmf.ffrg, vac.ffrg, cav.ffrg, pwie.ffrg]
        err = [pmf.stdv, vac.stdv, cav.stdv, pwie.stdv]
        ylabels = [r"$W(R_{g})$ [kT]", r"$W_{vac}(R_{g})$ [kT]", 
                    r"$W_{cav}(R_{g})$ [kT]", r"$E_{pw}(R_{g})$ [kT]"]
        alphas=[0.8, 1.0, 0.8, 0.8]

    m=0 # Row
    k=0 # Column
    for i in range(len(ylabels)):
        if title =="0arg" and i != 1:
            color = "black"
        else:
            color = colors[i][j]
        ax[m,k].errorbar(x[0], y[i], err[i], color=color, errorevery=step, alpha=alphas[i], ls='-', lw=3, elinewidth=2, capsize=2)
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
    
#%% Run

fig, ax = plt.subplots(2,3, figsize=(18,10))
plt.subplots_adjust(wspace=0.3, hspace=0.3)

j=0
delta_list = []
for title in titles:
    deltas = []
    path="C:/Research/.Projects/Virus-Stabilization/Hydrophobic-Polymer/REUS/arginine/{}/interactions".format(title)
    pmf, vac, cav, pwie, paie, pcie = open_file(path)
    plot()
    j+=1

plt.tight_layout()
plt.savefig('PMF-decomposition-{}.svg'.format('-'.join(titles)), transparent=True)

plt.show()


