import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.ticker as ticker
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from scipy.optimize import curve_fit

def initialize(t_total, tau, data, verbose: bool):
    global delta_tau, max_blocks, tau_max, iter_tau, s
    # t_total = 100000 # timesteps
    # tau = 1000 # ps; starting tau value
    delta_tau = tau #change in tau per iteration
    max_blocks = 5
    tau_max = t_total/max_blocks #max tau has max number of blocks
    iter_tau = int(tau_max/tau) #number of iterations of tau
    s = []
    error = do_analysis(tau, data, t_total, verbose)
    return error
    
#Figure parameters
hfont = {'fontname':'Helvetica'}
big = 22
medium = 16
small = 14
xsmall = 12
N = 500
alpha = 1.0

def do_analysis(tau, data, t_total, verbose: bool):
    RMS = get_RMS(data, t_total)
    tau, tau_list = block_average(data, tau, t_total, verbose)
    error = fit_line(tau_list, s, RMS, t_total, verbose)
    return error
    
def statistical_inefficiency(tau, data, verbose: bool):
    interval = tau
    block_avg_list = []
    block_std_list = []
    block1 = []
    num_blocks = round(len(data)/interval)
    # num_blocks = max_blocks
    if verbose == True:
        print('Number of Blocks is: ', num_blocks)
    for i in range(num_blocks):
        block = data[(i * interval):((i + 1) * interval)] 
        block1.append(block)
        block1 = np.array(block1, dtype=float)
        block_avg = np.mean(block1)
        block_avg_list.append(block_avg)
        block_std = np.std(block1)
        block_std_list.append(block_std)
        block1 = []
    
    deviation = np.sum(np.mean(np.asarray(block_avg_list)**2) - np.asarray(np.mean(block_avg_list))**2)

    sig_b = np.sqrt((1/(num_blocks)) * deviation)
    
    #s_tau = np.mean(block_std_list)
    s_tau = sig_b
    
    if verbose == True:
        print('Standard Deviation: ', s_tau)
        # print('Sig_sq_A_b: ', sig_ab_sq)
        # print('Sig_sq_A: ', sig_a_sq)
        # print('Sim_Avg: ', sim_avg)
    
    return s_tau

def get_RMS(data, t_total):
    data_sq = data ** 2
    RMS = np.sqrt((1/t_total)*np.sum(data_sq))
    return RMS

def block_average(data, tau, t_total, verbose: bool):
    tau_list = []
    for k in range(iter_tau):
        s_tau = statistical_inefficiency(tau, data, verbose)
        s.append(s_tau)
        tau_list.append(tau)
        tau = tau + delta_tau
    return tau, tau_list

def objective(x, a, b):
    return a * np.arctan(b*x)

def fit_line(tau_list, s, RMS, t_total, verbose: bool):
    popt, _ = curve_fit(objective, tau_list, s)
    a, b = popt
    x_line = np.arange(min(tau_list), max(tau_list), 1)
    y_line = objective(x_line, a, b)

    s_estimate = max(y_line)
    t_estimate = x_line[np.where(y_line == s_estimate)]
    if verbose == True:
        print("S Estimate: ", s_estimate)
        print("T Estimate: ", t_estimate)
    s_est_list = []
    for i in range(len(s)):
        s_est_list.append(s_estimate)
    
    if verbose == True:
        #Plot averages
        fig, ax = plt.subplots()
        ax.plot(tau_list, s, color='slategray', marker="s", ls="")
        ax.plot(x_line, y_line, ls='--', color = 'mediumseagreen')
        ax.plot(tau_list, s_est_list, ls='--', color = 'seagreen')
        # ax.annotate(round(s_estimate, 2), xy=(0, s_estimate + 0.5))
        plt.autoscale()
        ax.xaxis.set_major_locator(ticker.AutoLocator())
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.yaxis.set_major_locator(ticker.AutoLocator())
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        plt.xticks(fontsize=small, **hfont)
        plt.yticks(fontsize=small, **hfont)
        plt.title("Block Averaging", **hfont, fontsize=big)
        plt.xlabel(r"$\tau_{b}$ [Block Length]", **hfont, fontsize=big)
        plt.ylabel("$\sigma_{b}$", **hfont, fontsize=big)
        # ax.legend(bbox_to_anchor=(1, 1))
        plt.tight_layout()
        plt.show()
    
    return s_estimate



