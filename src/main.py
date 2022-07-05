" Data Analysis Pipeline for NR reflectometry "
" Author: @S.Winnall "

import sys
import csv
import pandas as pd
import numpy as np
import config
from fitRefnx import refnxOptimisation
from plotFits import plotNR

def getFile(path,nSkip,delim):
    return pd.read_csv(path, skiprows=nSkip, sep=delim, comment='#', na_values =' ', skip_blank_lines=True, encoding = "utf-8") # on_bad_lines='skip',


def getSampleInfo():
    
    instructionsDir = '../input/fitNR-Instructions.txt'
    
    with open(instructionsDir, newline='') as f:
        title = list(csv.reader(f))[0][0].split('=')[1]
    
    instructionsFile = getFile(path=instructionsDir, nSkip=1, delim=',')

    nFiles = len(instructionsFile)
    
    file_paths   = {}
    inputLabels  = {}
    lipids       = {}
    ratios       = {}
    contrastList = []
    
    for i in range(nFiles):
        
        contrast = instructionsFile["contrast"][i]
        contrastList.append(contrast)
        
        file_paths[contrast] = "../input/" + instructionsFile["fname"][i] + ".txt"
        inputLabels[contrast] = instructionsFile["fname"][i] 
        
        lipids[contrast] = instructionsFile["lipid"][i] 
        ratios[contrast] = instructionsFile["ratio"][i] 
    
    return title, file_paths, inputLabels, contrastList, lipids, ratios



def main():
    
    # get sample instructions 
    title, file_paths, inputLabels, contrastList, lipids, ratios = getSampleInfo()

    # define number of models to fit 
    nContrasts = len(contrastList)
    
    # initialise variables 
    Q         = {}
    expNR     = {}
    expNR_err = {}
    labels    = {}
    
    
    # for each contrast, store the experimental information 
    for i in range(nContrasts):
        contrast = contrastList[i] 
        
        if file_paths.get(contrast) is not None:
            experimentFile = getFile(path=file_paths.get(contrast), nSkip=1, delim='\t')
             
            Q[i] = experimentFile[experimentFile.columns.values[0]]
            expNR[i] = experimentFile[experimentFile.columns.values[1]]
            expNR_err[i] = experimentFile[experimentFile.columns.values[2]]
            labels[i] = inputLabels.get(contrast)       

 
    # get model datasets   
    qmin = min(Q.get(0))
    qmax = max(Q.get(0)) # update to max value in any dict 
    modelQ, modelNR, global_objective, redChiSq = refnxOptimisation(title, file_paths, qmin, qmax, lipids, ratios)

    # plot fits via own code for greater flexibility 
    if config.plotFinalObj == True:
        plotNR(Q,expNR,expNR_err,modelQ,modelNR,labels,title)
        
    
    # write global objective output to txt file 
    if config.writeGlobalObj == True:
        global_objective_file = open('../output/global-objective-'+title+'.txt', 'w')
        global_objective_file.write('\nReduced ChiSq = %f\n' %redChiSq)
        global_objective_file.write('\nCaclulated from: global_objective.chisqr()/(nPoints-nPars)\n')
        global_objective_file.write(str(global_objective))
        global_objective_file.close()
        
    
    return




if __name__ == '__main__':
    print('--------------------------------------------------------------------')
    print('Program NRfits - Fitting Program for Neutron reflectivity')
    print('Version 0.0.3, April 2022')
    print('Developed by Samuel Winnall. @ UoM')
    print('--------------------------------------------------------------------\n\n')
    main()
