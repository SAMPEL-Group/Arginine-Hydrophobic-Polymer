#%% Import Some Stuff!!
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import seaborn as sns
import numpy as np
import pandas as pd
import copy

# Figure Parameters
font = "Helvetica"
hfont = {'fontname':'Helvetica'}
big = 44
medium = 32
small = 14
xsmall = 10
alpha=0.2

rcut=0.6

runs = [1,2,3]
nW = 12

rpath="C:/Research/.Projects/Virus-Stabilization/Hydrophobic-Polymer/REUS/arginine/analysis/pref-int-new/revised/data/gyrate"
ppath="C:/Research/.Projects/Virus-Stabilization/Hydrophobic-Polymer/REUS/arginine/avg-pmfs"

#%% Functions and Grab Data
def openFile(title, window, run):
    data = pd.read_table("data/all-angles-time-3-{}-reus-{}-{}.dat".format(title, window, run), skiprows=1, sep=" ", skipinitialspace=True)
    rg = pd.read_table("../gyrate/gyrate-{}-reus-{}-{}.xvg".format(title, window, run), sep=" ", skiprows=27, skipinitialspace=True, names=['time', 'rg', 'rgx', 'rgy', 'rgz'])
    data = data.mean(axis=1)
    data = pd.concat((data, rg.rg), axis=1)
    data.columns = ['theta', 'rg']
    return data
    
def plot_histograms(histoData, ls):
    # Simple Plotting Function
    ax = sns.histplot(histoData,
                      bins=25,
                      kde=True,
                      stat='probability',
                      color=colors[i],
                      line_kws={'linestyle':ls, 'linewidth':4},
                      alpha=0,
                      lw=0)

colors = pl.cm.Reds(np.linspace(0.3,1,3))
titles=["47arg", "93arg", "185arg"]

# colors = pl.cm.Purples(np.linspace(0.3,1,3))
# titles=["48arg-glu", "94arg-glu", "186arg-glu"]

# colors = pl.cm.Oranges(np.linspace(0.3,1,3))
# titles=["48arg-lys", "94arg-lys", "186arg-lys"]
# titles=["185arg", "186arg-glu", "186arg-lys"]
# titles=["47arg"]

i=0
fig, ax = plt.subplots(1, 1, figsize=(12,8))
for title in titles:
    catData=[]
    for run in runs:
        for window in range(nW):
            data = openFile(title, window, run)
            catData.append(data)
    data = pd.concat(catData, axis=0)
    data.dropna(inplace=True)
    f, u = data[data.rg <= 0.6], data[data.rg > 0.6]
    plot_histograms(f.theta, '-')
    plot_histograms(u.theta, ':')
    i+=1

plt.xticks(fontsize=medium, **hfont)
plt.yticks(fontsize=medium, **hfont)
plt.xlabel(r"$\theta$", **hfont, fontsize=big)
plt.ylabel(r"P($\theta$)", **hfont, fontsize=big)

ax.xaxis.set_major_locator(ticker.MultipleLocator(30))
ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(15))
ax.yaxis.set_major_locator(ticker.AutoLocator())
ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

ax.tick_params(axis='both', which='both', direction='in', labelsize=medium)
ax.tick_params(axis='both', which='major', length=12, width=2)
ax.tick_params(axis='both', which='minor', length=8, width=2)

# plt.legend(loc="best")
plt.xlim(40,180)
plt.ylim(-0.005, 0.2)
plt.tight_layout()
plt.savefig('figures/orientation-uf-{}.svg'.format(title), transparent=True)
