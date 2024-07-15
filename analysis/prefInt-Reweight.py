#%% Import Stuff!!
import sys
sys.path.append('C:/Research/.Projects/Universal-Files')
import blockaveragingstdv
import dG
import pandas as pd
import numpy as np
import copy
import time
from tqdm import tqdm
from numpy import savetxt

begin = time.monotonic()

#User Inputs
nW = 12
nR = 3

#%% 
suffix="glu"
titles=["47{}".format(suffix), "93{}".format(suffix), "185{}".format(suffix)]
windows=np.arange(0,nW,1)
runs=np.arange(1,nR+1,1)

ppath="C:/Research/.Projects/Virus-Stabilization/Hydrophobic-Polymer/REUS/arginine/avg-pmfs"

binaryExc = False # Two excipients? I.e., arg + lys...
binaryIon = False # Two ions? I.e., arg + glu...

#%% Define some functions!!

def openFile(title, window, run):
    data = pd.read_csv("output/prefInt-{}-{}-{}.dat".format(title, window, run), sep=" ")
    rg = pd.read_table("../../gyrate/gyrate-{}-reus-{}-{}.xvg".format(title, window, run), sep=" ", skiprows=27, skipinitialspace=True, names=['time', 'rg', 'rgx', 'rgy', 'rgz'])
    return data, rg

def getPDF(title):
    data = pd.read_table("{}/average_pmf_{}.dat".format(ppath, title), skiprows=1, sep=" ", names=['rg', 'ffrg', 'stdv'])
    pmf = copy.deepcopy(data)
    pdf = np.exp(-1 * data.ffrg) / np.sum(np.exp(-1 * data.ffrg))
    pdf = pd.concat((data.rg, pdf), axis=1)
    pdf = pdf.set_index(np.arange(0, len(pdf.rg), 1))
    pdf.rename(columns={"ffrg" : "p"}, inplace=True)
    return pmf, pdf

def binData(cvData, pdf):
    outputCV, outputVal, outputErr = [], [], []
    lowerWalls = pdf.rg - (0.5*(pdf.rg[1] - pdf.rg[0]))
    upperWalls = pdf.rg + (0.5*(pdf.rg[1] - pdf.rg[0]))
    i=0
    for i in range(len(pdf)):
        currBin = cvData[(cvData.rg >= lowerWalls[i]) & (cvData.rg < upperWalls[i])]
        outputCV.append(pdf.rg[i])
        outputVal.append(np.mean(currBin.n))
        outputErr.append(np.std(currBin.n))
    outputData = pd.DataFrame(np.vstack((outputCV, outputVal, outputErr)).T, columns=["rg", "val", "stdv"])
    return outputData

def unbiasData(cvData, pdf):
    pdf = pdf[(pdf.rg >= 0.35) & (pdf.rg < 0.9)]
    pdf.p = pdf.p / np.sum(pdf.p)
    pdf = pdf.set_index(np.arange(0, len(pdf.rg), 1))
    
    binnedData = binData(cvData, pdf)
    
    folded = binnedData[binnedData['rg'] < rcutPMF]
    unfolded = binnedData[binnedData['rg'] >= rcutPMF]
    
    foldedProbability = pdf[pdf['rg'] < rcutPMF].p
    unfoldedProbability = pdf[pdf['rg'] >= rcutPMF].p
    
    foldedReweight = np.sum(folded.val*foldedProbability) / np.sum(foldedProbability)
    unfoldedReweight = np.sum(unfolded.val*unfoldedProbability) / np.sum(unfoldedProbability)
    
    # fErr = np.sum(folded.stdv*foldedProbability) / np.sum(foldedProbability) # error in bins of folded state
    # uErr = np.sum(unfolded.stdv*unfoldedProbability) / np.sum(unfoldedProbability) # error in bins of folded state
    
    fErr = np.sqrt(np.sum((folded.stdv**2)*(foldedProbability**2))) / np.sum(foldedProbability) # error in bins of folded state
    uErr = np.sqrt(np.sum((unfolded.stdv**2)*(unfoldedProbability**2))) / np.sum(unfoldedProbability) # error in bins of folded state

    
    delta = unfoldedReweight - foldedReweight
    err = np.sqrt(fErr**2 + uErr**2) # total dG error
    
    return delta, err, foldedReweight, fErr, unfoldedReweight, uErr


#%% Initialize data, split into folded and unfolded, and work through unbiasing protocl

for title in titles:
    runData=[]
    for run in runs:
        catData=[]
        catRg=[]
        rcutData=[]
        for window in windows:
            data, rg = openFile(title, window, run)
            catData.append(data)
            catRg.append(rg.rg)
        catData = pd.concat(catData, axis=0)
        catRg = pd.concat(catRg, axis=0)
        catRg.reset_index(drop=True, inplace=True)
        catRg = catRg.to_frame()
        
        for x in np.unique(data.rcut):
            tempData = catData[catData.rcut == x]
            tempData.reset_index(drop=True, inplace=True)
            rcutData.append(pd.concat((catRg, tempData), axis=1))
        runData.append(pd.concat(rcutData, axis=0))
    
    pmf, pdf = getPDF(title)
    dG_val, dG_err, rcutPMF = dG.get_dG(pmf, 1, 0.35, 0.90)
    exc1Output=pd.DataFrame(None, columns=["rcut", "delta", "err", "fRW", "fErr", "uRW", "uErr"])
    exc2Output=pd.DataFrame(None, columns=["rcut", "delta", "err", "fRW", "fErr", "uRW", "uErr"])

    ion1Output=pd.DataFrame(None, columns=["rcut", "delta", "err", "fRW", "fErr", "uRW", "uErr"])
    ion2Output=pd.DataFrame(None, columns=["rcut", "delta", "err", "fRW", "fErr", "uRW", "uErr"])

    for x in tqdm(np.unique(data.rcut)):
        for i in range(len(runs)):
            currData = runData[i][runData[i].rcut == x]
            
            cvData = pd.concat((currData.rg, currData.gammaExc1), axis=1)
            cvData.columns.values[1] = "n"
            unbiasedData = np.vstack((x, np.vstack((unbiasData(cvData, pdf))))).T
            output = pd.DataFrame(unbiasedData, columns=["rcut", "delta", "err", "fRW", "fErr", "uRW", "uErr"])
            exc1Output = pd.concat((exc1Output, output), axis=0)
            
            if binaryExc == True:
                cvData = pd.concat((currData.rg, currData.gammaExc2), axis=1)
                cvData.columns.values[1] = "n"
                unbiasedData = np.vstack((x, np.vstack((unbiasData(cvData, pdf))))).T
                output = pd.DataFrame(unbiasedData, columns=["rcut", "delta", "err", "fRW", "fErr", "uRW", "uErr"])
                exc2Output = pd.concat((exc2Output, output), axis=0)
            
            currData = runData[i][runData[i].rcut == x]
            cvData = pd.concat((currData.rg, currData.gammaIon1), axis=1)
            cvData.columns.values[1] = "n"
            unbiasedData = np.vstack((x, np.vstack((unbiasData(cvData, pdf))))).T
            output = pd.DataFrame(unbiasedData, columns=["rcut", "delta", "err", "fRW", "fErr", "uRW", "uErr"])
            ion1Output = pd.concat((ion1Output, output), axis=0)
            
            if binaryIon == True:
                currData = runData[i][runData[i].rcut == x]
                cvData = pd.concat((currData.rg, currData.gammaIon2), axis=1)
                cvData.columns.values[1] = "n"
                unbiasedData = np.vstack((x, np.vstack((unbiasData(cvData, pdf))))).T
                output = pd.DataFrame(unbiasedData, columns=["rcut", "delta", "err", "fRW", "fErr", "uRW", "uErr"])
                ion2Output = pd.concat((ion2Output, output), axis=0)

    avgExc1Output=pd.DataFrame(None, columns=["rcut", "delta", "err", "fRW", "fErr", "uRW", "uErr"])
    errExc1Output=pd.DataFrame(None, columns=["rcut", "delta", "err", "fRW", "fErr", "uRW", "uErr"])

    if binaryExc == True:
        avgExc2Output=pd.DataFrame(None, columns=["rcut", "delta", "err", "fRW", "fErr", "uRW", "uErr"])
        errExc2Output=pd.DataFrame(None, columns=["rcut", "delta", "err", "fRW", "fErr", "uRW", "uErr"])

    avgIon1Output=pd.DataFrame(None, columns=["rcut", "delta", "err", "fRW", "fErr", "uRW", "uErr"])
    errIon1Output=pd.DataFrame(None, columns=["rcut", "delta", "err", "fRW", "fErr", "uRW", "uErr"])

    if binaryIon == True:
        avgIon2Output=pd.DataFrame(None, columns=["rcut", "delta", "err", "fRW", "fErr", "uRW", "uErr"])
        errIon2Output=pd.DataFrame(None, columns=["rcut", "delta", "err", "fRW", "fErr", "uRW", "uErr"])


    for x in tqdm(np.unique(data.rcut)):
        currData = exc1Output[exc1Output.rcut == x]
        avgData = pd.DataFrame(currData.mean(axis=0)).transpose()
        avgExc1Output = pd.concat((avgExc1Output, avgData), axis=0)
        errData = pd.DataFrame(currData.std(axis=0)).transpose()
        errExc1Output = pd.concat((errExc1Output, errData), axis=0)
        
        if binaryExc == True:
            currData = exc2Output[exc2Output.rcut == x]
            avgData = pd.DataFrame(currData.mean(axis=0)).transpose()
            avgExc2Output = pd.concat((avgExc2Output, avgData), axis=0)
            errData = pd.DataFrame(currData.std(axis=0)).transpose()
            errExc2Output = pd.concat((errExc2Output, errData), axis=0)
        
        currData = ion1Output[ion1Output.rcut == x]
        avgData = pd.DataFrame(currData.mean(axis=0)).transpose()
        avgIon1Output = pd.concat((avgIon1Output, avgData), axis=0)
        errData = pd.DataFrame(currData.std(axis=0)).transpose()
        errIon1Output = pd.concat((errIon1Output, errData), axis=0)
        
        if binaryIon == True:
            currData = ion2Output[ion2Output.rcut == x]
            avgData = pd.DataFrame(currData.mean(axis=0)).transpose()
            avgIon2Output = pd.concat((avgIon2Output, avgData), axis=0)
            errData = pd.DataFrame(currData.std(axis=0)).transpose()
            errIon2Output = pd.concat((errIon2Output, errData), axis=0)

    avgExc1Output.err = errExc1Output.delta
    avgExc1Output.fErr = errExc1Output.fRW
    avgExc1Output.uErr = errExc1Output.uRW
    
    if binaryExc == True:
        avgExc2Output.err = errExc2Output.delta
        avgExc2Output.fErr = errExc2Output.fRW
        avgExc2Output.uErr = errExc2Output.uRW
    
    avgIon1Output.err = errIon1Output.delta
    avgIon1Output.fErr = errIon1Output.fRW
    avgIon1Output.uErr = errIon1Output.uRW
    
    if binaryIon == True:
        avgIon2Output.err = errIon2Output.delta
        avgIon2Output.fErr = errIon2Output.fRW
        avgIon2Output.uErr = errIon2Output.uRW
    
    header = "rcut delta err fRW fErr uRW uErr"
    
    savetxt('reweighted/prefInt-reweighted-{}-exc1.dat'.format(title), avgExc1Output, header=header, delimiter=' ', fmt='%.3f', comments='') 
    if binaryExc == True:
        savetxt('reweighted/prefInt-reweighted-{}-exc2.dat'.format(title), avgExc2Output, header=header, delimiter=' ', fmt='%.3f', comments='') 

    savetxt('reweighted/prefInt-reweighted-{}-ion1.dat'.format(title), avgIon1Output, header=header, delimiter=' ', fmt='%.3f', comments='') 
    if binaryIon == True:
        savetxt('reweighted/prefInt-reweighted-{}-ion2.dat'.format(title), avgIon2Output, header=header, delimiter=' ', fmt='%.3f', comments='') 

print(f"Time elapsed: {time.monotonic() - begin:.2f} second(s)")
