" Data Analysis Pipeline for the NR reflectometry "
" Author: @S.Winnall "

import pandas as pd
import numpy as np
import config
from fitting import geneticAlgo, leastSquares
from plotting import geneticAnalysis, plotCost, plotThickness


def getFile(path,nSkip,delim):
    return pd.read_csv(path, skiprows=nSkip, sep=delim, comment='#', na_values =' ', skip_blank_lines=True, encoding = "utf-8") # on_bad_lines='skip',


def main():
    
    if config.verbose == False:
        print("Verbose = False, terminal silenced.")
    
    # read file into memory data
    inputFiles  = ['../input/S11_excel.txt']#,'../input/S14_excel.txt']
    inputLabels = ['d-MC3 ACMW']#,'h-MC3 D2O']
   
    nModels = len(inputFiles) 
       
    Q     = {init: [] for init in range(nModels)}
    expNR = {init: [] for init in range(nModels)}
    
    # for each model store list of genetic outputs
    macroData = {init: [] for init in range(nModels)}
    

    # define which head group thicknesses are to be studied 
    headParList = []
    for num in np.arange(5.0, 15.5, 1.0):
        headParList.append(num)
    nThickness  = len(headParList)

    # iterate over every file to get experiment data
    for modelNum, fileDIR in enumerate(inputFiles):
    
        # get sample data as pandas df
        experimentFile = getFile(path=fileDIR, nSkip=0, delim='\t')
    
        # parse variables
        Q[modelNum]     = experimentFile[experimentFile.columns.values[0]]
        expNR[modelNum] = experimentFile[experimentFile.columns.values[1]]
        
        # for each file/model, fit d1 for given d2, N times 
        for thicknessIdx, headPar in enumerate(headParList):
            
            # numnber of repeats 
            N = 1
            for n in range(N):  

                # run fitting module
                maxIter = 1
                res = geneticAlgo(maxIter, Q.get(modelNum), expNR.get(modelNum), modelNum, headPar)
                
                # initialise geneticOutput dictionary 
                geneticOutput = {}
                
                # model number 
                geneticOutput["modelNum"] = modelNum 
                
                # head group thickness 
                geneticOutput["d2"] = headPar
                
                # tail thickness 
                geneticOutput["d1"] = res.x[7]
                
                # parameter solution 
                geneticOutput["parSolution"] = res.x
                
                #print('\n\ParSolution:\n%s\n\n' %res.x)
                
                # associated cost
                expData = (Q.get(modelNum), expNR.get(modelNum))
                geneticOutput["Cost"] = leastSquares(res.x, *expData)
                
                # append dict to a list in a macro dict of each model 
                #print('\ngeneticOutput:\n%s\n' %geneticOutput)
                macroData[modelNum].append(geneticOutput)
                #print('macro:\n%s\n' %macroData.get(modelNum))
    
    #print('\n\nMacroData:\n%s\n\n' %macroData)


    # print output parameters, statistics and a figure
    #geneticAnalysis(geneticOutputData, Q, expNR, inputLabels)
    
    plotCost(N, macroData, Q, expNR, inputLabels)
    
    # for fits where  d2 = [6, 6.1, 6.2 ... 15] and 
    plotThickness(N, macroData, Q, expNR, headParList, inputLabels)

    return




if __name__ == '__main__':
    print('--------------------------------------------------------------------')
    print('Program NRfits - Fit Program for Neutron reflection datasets')
    print('Version 0.0.2, April 2022')
    print('Developed by Samuel Winnall. @ UoM')
    print('--------------------------------------------------------------------\n\n')
    main()
