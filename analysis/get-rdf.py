#%% What needs to be done?
# 1. Select an arginine molecule, within a cutoff distance to the polymer
# 2. Identify the nearest polymer bead to arginine
# 3. Compute the vector connecting polymer bead to arginine COM
# 4. Compute the vector connecting the arginine COM to Gdm+ C
# 5. Compute the angle between the two vectors

#%% Import Some Stuff
import sys
sys.path.append('C:/Research/.Projects/Universal-Files')
import blockaveragingstdv
import numpy as np
import MDAnalysis as md
from MDAnalysis.analysis import contacts
from MDAnalysis.analysis.rdf import InterRDF
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import matplotlib.ticker as ticker
import matplotlib.pylab as pl

path="C:/Research/.Projects/Virus-Stabilization/Arginine/Mixtures/"

#%% User Inputs
n = 47 # USER INPUT: Number of excipient molecules
mol2 = "Arginine" # USER INPUT: Select second excipient
sel1 = "CZ"
sel2 = "CZ"

# Figure Parameters
font = "Helvetica"
hfont = {'fontname':'Helvetica'}
big = 22
medium = 18
small = 14
xsmall = 10
alpha=0.8
bins=25
fig, ax = plt.subplots()

#%% Extract Simulation Data
if mol2 == "Arginine":    
    title = "arg-cl-{}"
    res = "ARG"
if mol2 == "Glutamate":
    title = "arg-cl-glu-na-{}"
    res = "GLU"
if mol2 == "Lysine":
    title = "arg-lys-cl-{}"
    res = "LYS"

def open_files(mol1, mol2, title):
    gro = r"{}/{}-{}/gro/tip4p05-26n-box-{}-npt-prod.gro".format(path, mol1, mol2, title)
    traj = r"{}/{}-{}/traj/tip4p05-26n-box-{}-npt-prod-trunc.xtc".format(path, mol1, mol2, title)
    u = md.Universe(gro, traj)
    
    return u

u = open_files("Arginine", "{}".format(mol2), title.format(n))

def get_rdf():
    g1 = u.select_atoms('resname ARG and name {}'.format(sel1))
    g2 = u.select_atoms('resname {} and name {}'.format(res, sel2))
    
    rdf = InterRDF(g1, g2, nbins=75, range=(3,30))
    rdf.run(verbose=True)
    
    plt.plot(rdf.bins/10, rdf.rdf, lw=2, alpha=alpha)
    
    return rdf
    
def plot_finish():
    plt.autoscale()
    ax.xaxis.set_major_locator(ticker.AutoLocator())
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.05))
    ax.yaxis.set_major_locator(ticker.AutoLocator())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    # plt.xlim(0, 1.5)
    # plt.ylim(0, 1.5)
    plt.xticks(fontsize=small)
    plt.yticks(fontsize=small)
    plt.title("Arginine({})-{}({}) RDF".format(sel1, mol2.capitalize(), sel2), **hfont, fontsize=big)
    plt.xlabel("r [nm]", **hfont, fontsize=big)
    plt.ylabel("g(r)", **hfont, fontsize=big)
    # ax.legend(prop=font_manager.FontProperties(family=font), loc="best", fontsize=small)
    # plt.legend(loc="best")
    # ax.legend(bbox_to_anchor=(1.05, 0.5), loc="upper left")
    plt.tight_layout()
    plt.savefig('{}-{}-{}-rdf.svg'.format(mol2, sel1, sel2), transparent=True)
    plt.show()

# colors = ["firebrick", "darkcyan", "indigo"]
# states = ["collapsed", "elongated"]
# rdf_list = []
# for i in range(len(states)):
#     u = open_files("{}".format(states[i]))
#     rdf = get_rdf()
#     rdf_list.append(rdf)

rdf = get_rdf()
plot_finish()


