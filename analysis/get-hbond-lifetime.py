#%% What needs to be done?
# 1. Determine hydrogen bonding pairs
# 2. Compute the number of hydrogen bonds per timestep
# 3. Compute the lifetime of hydrogen bonds

#%% Import Some Stuff
import numpy as np
import MDAnalysis as md
from MDAnalysis.analysis import contacts
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
import MDAnalysis.analysis.hydrogenbonds.hbond_autocorrel as HBAC
from MDAnalysis.analysis.waterdynamics import WaterOrientationalRelaxation as WOR
from tqdm import tqdm
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.ticker as ticker
from numpy import savetxt

#%% User Inputs

k=0
title = "Arginine"
path="C:/Research/.Projects/Virus-Stabilization/Arginine/Mixtures/Arginine-{}".format(title)

# Figure Parameters
font = "Helvetica"
hfont = {'fontname':'Helvetica'}
big = 22
medium = 18
small = 14
xsmall = 10
alpha=0.2
bins=25

#%% Extract Simulation Data

def open_files():
    gro = r"{}/gro/tip4p05-26n-box-arg-cl-{}-npt-prod.gro".format(path, k)
    # traj = r"{}/traj/tip4p05-26n-box-arg-cl-{}-npt-prod-trunc.xtc".format(path, k)
    traj = r"{}/traj/tip4p05-26n-box-arg-cl-{}-100ps.xtc".format(path, k)
    top = r"{}/tpr/tip4p05-26n-box-arg-cl-{}-npt-prod.tpr".format(path, k)
    u = md.Universe(top, traj)
    
    return u

#%% Plotting functions
def plot():
    fig, ax = plt.subplots()

    plt.plot(frames, hbond_res, color="black", alpha=0.8, lw=2)

    plt.autoscale()
    ax.xaxis.set_major_locator(ticker.AutoLocator())
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
    ax.yaxis.set_major_locator(ticker.AutoLocator())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    plt.xticks(fontsize=small, **hfont)
    plt.yticks(fontsize=small, **hfont)
    plt.xlabel("Time [ns]", **hfont, fontsize=big)
    plt.ylabel("$N^{WW}_{HB}$", **hfont, fontsize=big)

    plt.tight_layout()
    plt.savefig('ww-hb-{}.svg'.format(title), transparent=True)
    plt.show()

#%% Exponential Fit
def fit_biexponential(tau_timeseries, ac_timeseries):
    from scipy.optimize import curve_fit

    def model(t, A, tau1, B, tau2):
        return A * np.exp(-t / tau1) + B * np.exp(-t / tau2)

    params, params_covariance = curve_fit(model, tau_timeseries, ac_timeseries, [1, 0.5, 1, 2])

    fit_t = np.linspace(tau_timeseries[0], tau_timeseries[-1], 1000)
    fit_ac = model(fit_t, *params)

    return params, fit_t, fit_ac
#%%

u = open_files()

arg_atoms = u.select_atoms('resname ARG') 
wat_atoms = u.select_atoms("resname water")
w = len(wat_atoms)/4
e = len(arg_atoms)/27

sel1 = "name OW"
sel2 = "name OW"

import time
t1 = time.perf_counter()
hbonds = HBA(u, update_selections=False, d_a_cutoff=3.5, d_h_a_angle_cutoff=150,
             acceptors_sel="{}".format(sel1), donors_sel="{}".format(sel2))
hbonds.run()

# hbonds.donors_sel = "name OW"
# hbonds.acceptors_sel = "name OW"

hbonds._acceptors.names
hbonds.acceptors_sel
hbonds.hydrogens_sel

# hbonds.guess_acceptors("resname water")

hbonds_t = hbonds.count_by_time()
hbond_res = 2 * (hbonds_t / w)
results = hbonds.results.hbonds
frames = hbonds.times / 1000
plot()

# plt.clf()
# hblt = hbonds.lifetime().T
# tau_times = hblt[:,0] * u.trajectory.dt
# plt.plot(tau_times, hblt[:,1])

# params, fit_t, fit_ac = fit_biexponential(hblt[:,0], hblt[:,1])
# A, tau1, B, tau2 = params
# tc = A * tau1 + B * tau2
# print(f"Time Constant = {tc:.2f} ps")

t2 = time.perf_counter()
print(f"Took {t2 - t1:0.4f} Seconds to Complete Analysis")
