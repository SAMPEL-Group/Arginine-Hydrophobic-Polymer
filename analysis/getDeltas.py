
#%% Import Stuff
import numpy as np
import pandas as pd
from numpy import savetxt
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import matplotlib.ticker as ticker
import matplotlib.pylab as pl
import tqdm as tqdm
import dG

#%% Define Some Variables
T = 300
kBT = 1

font = "Helvetica"
hfont = {'fontname':'Helvetica'}
big = 22
medium = 18
small = 14
xsmall = 10

opath = "C:/Research/.Projects/Virus-Stabilization/Hydrophobic-Polymer/REUS/arginine/0arg/interactions/deltas"

excluded_titles = ["47gly", "93gly", "185gly"] # Exclusions for pcie
labels=["PMF", "Vac", "Cav", "Epw", "Epa", "Epc"]
spe=False
orig=False

#%% Functions

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

def getPMF():
    #Get probability data
    new_path = "C:/Research/.Projects/Virus-Stabilization/Hydrophobic-Polymer/REUS/arginine/avg-pmfs"
    data = pd.read_table("{}/average_pmf_{}.dat".format(new_path,title), skiprows=1, sep=" ", names=['rg', 'ffrg', 'stdv'])
    data = data.set_index(np.arange(0, len(data.rg), 1))
    pdf = np.exp(-1 * data.ffrg) / np.sum(np.exp(-1 * data.ffrg))
    pdf = pd.concat((data.rg, pdf), axis=1)
    pdf = pdf.set_index(np.arange(0, len(pdf.rg), 1))
    return pdf

### Obtain delta by integration
def delta(energy):
    folded = energy[energy['rg'] < rcut].ffrg
    unfolded = energy[energy['rg'] >= rcut].ffrg
    
    num = np.sum(np.exp(-1 * unfolded / kBT))
    den = np.sum(np.exp(-1 * folded / kBT))
    delta = np.log(num / den) * -1 * kBT
    
    unf_err = energy[energy['rg'] >= rcut].ffrg + energy[energy['rg'] >= rcut].stdv
    fol_err = energy[energy['rg'] < rcut].ffrg + energy[energy['rg'] < rcut].stdv
    en = np.sum(np.exp(-1 * unf_err / kBT))
    ed = np.sum(np.exp(-1 * fol_err / kBT))
    de = np.log(en / ed) * -1 * kBT
    
    err = np.abs(delta - de)

    return delta, err

### Obtain delta by.... other means
def deltaIE(energy):
    folded = energy[energy['rg'] < rcut]
    unfolded = energy[energy['rg'] >= rcut]
    
    foldedProbability = pdf[pdf['rg'] < rcut].ffrg
    unfoldedProbability = pdf[pdf['rg'] >= rcut].ffrg
    
    foldedReweight = np.sum(folded.ffrg*foldedProbability) / np.sum(foldedProbability)
    unfoldedReweight = np.sum(unfolded.ffrg*unfoldedProbability) / np.sum(unfoldedProbability)
    
    # dRg = energy.rg[1]-energy.rg[0] # bin width
    # fErr = folded.stdv # error in bins of folded state
    # ifErr = dRg * np.sqrt(np.sum(np.square(fErr))) # integration error for folded state
    # uErr = folded.stdv # error in bins of unfolded state
    # iuErr = dRg * np.sqrt(np.sum(np.square(uErr))) # integration error for unfolded state
    
    # fErr = np.sum(folded.stdv*foldedProbability) / np.sum(foldedProbability)
    # uErr = np.sum(unfolded.stdv*unfoldedProbability) / np.sum(unfoldedProbability)
    
    # new way 01/17/2024
    fErr = np.sqrt(np.sum((folded.stdv**2)*(foldedProbability**2))) / np.sum(foldedProbability)
    uErr = np.sqrt(np.sum((unfolded.stdv**2)*(unfoldedProbability**2))) / np.sum(unfoldedProbability)

    err = np.sqrt(fErr**2 + uErr**2) # total dG error
    
    delta = unfoldedReweight - foldedReweight
    
    #%% Old methods
    # Original integration method
    if orig == True:
        energy = energy[(energy['rg'] > 0.35) & (energy['rg'] < 0.9)]
        folded = energy[energy['rg'] < rcut].ffrg
        unfolded = energy[energy['rg'] >= rcut].ffrg
        
        num = np.sum(np.exp(-1 * unfolded / kBT))
        den = np.sum(np.exp(-1 * folded / kBT))
        delta = np.log(num/den)*-1*kBT
        
        unf_err = energy[energy['rg'] >= rcut].ffrg + energy[energy['rg'] >= rcut].stdv
        fol_err = energy[energy['rg'] < rcut].ffrg + energy[energy['rg'] < rcut].stdv
        en = np.sum(np.exp(-1 * unf_err / kBT))
        ed = np.sum(np.exp(-1 * fol_err / kBT))
        de = np.log(en / ed) * -1 * kBT
        
        err = np.abs(delta - de)
    
    # Single point
    if spe == True:
        folded = energy.iloc[(energy['rg']-0.4).abs().argsort()[:1]].ffrg
        unfolded = energy.iloc[(energy['rg']-0.8).abs().argsort()[:1]].ffrg
        delta = (np.asarray(unfolded) - np.asarray(folded))[0]
        
        folded = energy.iloc[(energy['rg']-0.4).abs().argsort()[:1]].stdv
        unfolded = energy.iloc[(energy['rg']-0.8).abs().argsort()[:1]].stdv
        err = (np.asarray(unfolded) + np.asarray(folded))[0]
    #%%
    return delta, err

#%% Run

j=0
delta_list = []
errs_list = []
titles=["0arg","47arg","93arg","185arg","47glu","93glu","185glu",
        "47lys","93lys","185lys","48arg-glu","94arg-glu","186arg-glu",
        "48arg-lys","94arg-lys","186arg-lys", "48lys-glu", "94lys-glu",
        "186lys-glu", "47gly","93gly","185gly","47gdm","93gdm","185gdm"]
for title in titles:
    deltas = []
    errs = []
    path="C:/Research/.Projects/Virus-Stabilization/Hydrophobic-Polymer/REUS/arginine/{}/interactions".format(title)
    pmf, vac, cav, pwie, paie, pcie = open_file(path)
    pdf=getPMF()
    
    #%% append things
    dG_val, dG_err, rcut = dG.get_dG(pmf, kBT, 0.35, 0.90)
    deltas.append(dG_val)
    deltas.append(delta(vac)[0])
    deltas.append(deltaIE(cav)[0])
    deltas.append(deltaIE(pwie)[0])
    errs.append(dG_err)
    errs.append(delta(vac)[1])
    errs.append(deltaIE(cav)[1])
    errs.append(deltaIE(pwie)[1])

    if title != "0arg":
        deltas.append(deltaIE(paie)[0])
        errs.append(deltaIE(paie)[1])
        
        if title not in excluded_titles:
            deltas.append(deltaIE(pcie)[0])
            errs.append(deltaIE(pcie)[1])
    
    if title == "0arg":
        deltas.append(paie)
        errs.append(paie)
        deltas.append(pcie)
        errs.append(pcie)
        
    if title in excluded_titles:
        deltas.append(pcie)
        errs.append(pcie)
            
    delta_list.append(deltas)
    errs_list.append(errs)
    #%%
    j+=1

    print("==={}===".format(title))
    # append values to output
    output = []
    for i in range(len(deltas)):
        row = labels[i], np.round(deltas[i], 2), np.round(errs[i], 3)
        print(labels[i], np.round(deltas[i], 2), np.round(errs[i], 3))
        output.append(list(row))
    
    # append sum
    row = "sum", np.round(np.sum(deltas[1:]), 2), np.round(np.sqrt(np.sum(np.square(errs[1:])))/len(errs[1:]), 3)
    output.append(list(row))
    print("sum", np.round(np.sum(deltas[1:]), 2), np.round(np.sqrt(np.sum(np.square(errs[1:])))/len(errs[1:]), 3))
    print("")

    header = "label val err"
    savetxt('output/deltas-{}.dat'.format(title), output, header=header, delimiter=' ', fmt='%s') 