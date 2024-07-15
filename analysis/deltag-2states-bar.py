#%% Import Stuff
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import matplotlib.font_manager as font_manager
import matplotlib.ticker as ticker
import math

T = 300
kBT = 1
#Defaults
rcut = 0.65 #Cutoff to determine folded/unfolded regimes
rmin = 0.35
rmax = 0.90

#Figure Parameters
font = "Helvetica"
hfont = {'fontname':'Helvetica'}
big = 36
medium = 28
small = 28
xsmall = 10
alpha=1.0

path = "C:/Research/.Projects/Virus-Stabilization/Hydrophobic-Polymer/REUS/arginine/avg-pmfs"

#User Inputs
conc = "both" # 0.25, 1.0, "both"
option = 4 # 1 --> single components, 2 --> mixtures, 3 --> individual species, 4 --> single to mixture
mols = 2 # 1 --> arg, 2 --> glu, 3 --> lys, 4 --> arg/glu, 5 --> arg/lys
manual = True # Option to override everything

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
    labels = ["0.0 M", "0.25 M", "0.5 M", "1.0 M"]
    titles = ["0arg", "48arg-glu", "94arg-glu", "186arg-glu"]
    colors = ["black", "plum", "darkorchid", "indigo"]

if conc == "both" and mols == 6:
    labels = ["0.0 M", "0.25 M", "0.5 M", "1.0 M"]
    titles = ["0arg", "48lys-glu", "94lys-glu", "186lys-glu"]
    colors = ["black", "mediumseagreen", "forestgreen", "darkgreen"]

if conc == "both" and mols == 5:
    labels = ["0.0 M", "0.25 M", "1.0 M"]
    titles = ["0arg", "48arg-lys", "94arg-lys", "186arg-lys"]
    colors = ["black", "peachpuff", "sandybrown", "darkorange"]
    
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
    labels = ["Wat", "Arg", "Arg/Glu"]
    titles = ["0arg", "93arg", "186arg-glu"]
    colors = ["black", "lightcoral", "darkgreen"]
    
#%%Manual
if manual == True:
    watColor = ["black"]
    argColors = ["lightcoral", "firebrick", "maroon"]
    gluColors = ["paleturquoise", "darkturquoise", "darkcyan"]
    lysColors = ["khaki", "goldenrod", "darkgoldenrod"]
    lysgluColors = ["palegreen", "mediumseagreen", "darkgreen"]
    arggluColors = ["thistle", "mediumorchid", "darkviolet"]
    arglysColors = ["bisque", "orange", "darkorange"]
    
    watLabel = ["0.0 M"]
    argLabels = ["0.25 M Arg", "0.5 M Arg", "1.0 M Arg"]
    gluLabels = ["0.25 M Glu", "0.5 M Glu", "1.0 M Glu"]
    lysLabels = ["0.25 M Lys", "0.5 M Lys", "1.0 M Lys"]
    arggluLabels = ["0.25 M Arg/Glu", "0.5 M Arg/Glu", "1.0 M Arg/Glu"]
    arglysLabels = ["0.25 M Arg/Lys", "0.5 M Arg/Lys", "1.0 M Arg/Lys"]
    lysgluLabels = ["0.25 M Lys/Glu", "0.5 M Lys/Glu", "1.0 M Lys/Glu"]
    
    watTitle = ["0arg"]
    argTitles = ["47arg", "93arg", "185arg"]
    gluTitles = ["47glu", "93glu", "185glu"]
    lysTitles = ["47lys", "93lys", "185lys"]
    arggluTitles = ["48arg-glu", "94arg-glu", "186arg-glu"]
    arglysTitles = ["48arg-lys", "94arg-lys", "186arg-lys"]
    lysgluTitles = ["48lys-glu", "94lys-glu", "186lys-glu"]
    
    # argTitles = ["47gly", "93gly", "185gly"]
    
    # colors = watColor + argColors + gluColors + lysColors + arggluColors + arglysColors + lysgluColors
    # labels = watLabel + argLabels + gluLabels + lysLabels + arggluLabels + arglysLabels + lysgluLabels
    # titles = watTitle + argTitles + gluTitles + lysTitles + arggluTitles + arglysTitles + lysgluTitles
    
    # colors = ["firebrick", "firebrick", "firebrick"]
    # hatches = ["", "/", "x"]
    # alphas = [1.0, 0.8, 0.8]
    # labels = [r"1.0 M Arg", r"1.0 M Arg$_{GG}$", r"1.0 M Arg$_{HTT}$"]
    # titles = ["185arg", "185arg-gdm", "185arg-htt"]
    
    colors = ["black", "lightcoral", "maroon"]
    hatches = ["", "", ""]
    alphas = [1.0, 1.0, 1.0]
    labels = [r"0.0 M", r"0.25 M Arg", r"1.0 M Arg"]
    titles = ["0arg-c", "47arg-c", "185arg-c"]
    
#%% Define Some Functions

def open_file():
    data = pd.read_table("{}/average_pmf_{}.dat".format(path, title), skiprows=1, sep=" ", names=['rg', 'ffrg', 'stdv'])
    data = data.set_index(np.arange(0, len(data.rg), 1))
    return data

def addlabels(x,y):
    for i in range(len(x)):
        if math.copysign(1, y[i]) == -1:
            factor = 1.5
        else:
            factor = 1
        plt.text(x[i], y[i]+(0.07*math.copysign(1, y[i])*factor), np.round(y[i], 3), ha = 'center',
                 bbox = dict(facecolor = 'white', alpha =0.7),
                 fontsize=medium, **hfont)

#%% Let's get dG!

def get_dG(ffrg, kBT, rmin, rmax):
    rcut=0.7
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
    print(rcut)
    
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

# more_offsets = [0,3,6,9,12,15]
more_offsets = [0,1]
# more_offsets = []
def plot_dg():
    global x_list
    fig, ax = plt.subplots(1, 1, figsize=(8,10))
    x_list = []
    width=1
    x = 0
    i=0
    count=0
    for i in range(len(o_result)):
        offset=width*count
        if i in more_offsets:
            count += 0.25
        # ax.bar(x + offset, o_result[i], yerr=o_err[i], ecolor="black", capsize=10, color=colors[i], width=width)
        ax.bar(x + offset, o_result[i], yerr=o_err[i], ecolor="black", capsize=10, color=colors[i], width=width, hatch=hatches[i], alpha=alphas[i])
        x_list.append(x+offset)
        count += 1
    
    # for i, v in enumerate(o_result):
    #     plt.text(x_list[i], v + (np.sign(v)*0.01), np.round(v, 2), horizontalalignment="center", verticalalignment="bottom", **hfont, fontsize=24)
    #     print(np.sign(v)*0.1)
    
    addlabels(x_list, o_result)
        
    ax.xaxis.set_major_locator(ticker.AutoLocator())
    ax.yaxis.set_major_locator(ticker.AutoLocator())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.tick_params(axis='both', which='both', direction='in', labelsize=medium)
    ax.tick_params(axis='both', which='major', length=12, width=2)
    ax.tick_params(axis='both', which='minor', length=8, width=2)
    ax.set_ylabel("$\Delta$G$_{u}$ [kT]", **hfont, fontsize=big)
    plt.yticks(fontsize=big, **hfont)
    ax.set_xticks(x_list, labels, fontsize=big, **hfont, rotation=30)
    
    for tick in ax.xaxis.get_majorticklabels():
        tick.set_horizontalalignment("right")
        
    # labels=ax.get_xticklabels() + ax.get_yticklabels()
    # [label.set_fontname(font) for label in labels]

    
i=1

o_result = []
o_err = []
for i in range(len(titles)):
    title = titles[i]
    ffrg = open_file()
    result, yerr, rcut = get_dG(ffrg, kBT, 0.35, 0.90)
    o_result.append(result)
    o_err.append(yerr)


# plt.xlim(0.35, 0.9)
# plt.ylim(-1, 3)
# plt.tight_layout()
# plt.show()

plot_dg()

plt.hlines(0, -1,22, color="black", lw=4)
plt.xlim(-1,3.2)
# plt.ylim(-0.5, 1.3)
plt.ylim(-1.1,0.4)

plt.tight_layout()
plt.savefig('dG-{}.svg'.format('-'.join(titles)), transparent=True)

plt.show()
