" Data Analysis Pipeline for NR reflectometry "
" Author: @S.Winnall "

import pandas as pd
import numpy as np
import config
from fitting import geneticAlgo, leastSquares
from genModelData import genModelNR
from plotting import geneticAnalysis, plotCost, plotThickness, plotColorMap, plotFits


def getFile(path,nSkip,delim):
    return pd.read_csv(path, skiprows=nSkip, sep=delim, comment='#', na_values =' ', skip_blank_lines=True, encoding = "utf-8") # on_bad_lines='skip',


def getMinModelData(plotAllModels,nModels,modelNumber,costList,macroData,Q):

    # minimise macroData with respect to cost
    minModelQ  = {init: [] for init in range(nModels)}
    minModelNR = {init: [] for init in range(nModels)}

    if plotAllModels == True:
        for modelNum in range(nModels):

            justCosts = []
            justIDs   = []
            for cost in costList.get(modelNum):
                justCosts.append(cost[0])
                justIDs.append(cost[1])

            minCost    = min(justCosts)
            minCostID  = justCosts.index(minCost)
            minMacroID = justIDs[minCostID]

            justCosts.clear()
            justIDs.clear()

            # iterate along the list of macros for given model and get list index of corresponding min macro
            for id, macro in enumerate(macroData.get(modelNum)):
                    if macro.get("macroID") == minMacroID:
                        minMacroID = id

            minMacro = macroData.get(modelNum)[minMacroID]
            minPar   = minMacro.get("parSolution")

            minModelQ[modelNum], minModelNR[modelNum] = genModelNR(minPar,Q.get(modelNum))

    else:
        justCosts = []
        justIDs   = []
        modelNum  = modelNumber
        for cost in costList.get(modelNum):
            justCosts.append(cost[0])
            justIDs.append(cost[1])

        minCost    = min(justCosts)
        minCostID  = justCosts.index(minCost)
        minMacroID = justIDs[minCostID]

        justCosts.clear()
        justIDs.clear()

        # iterate along the list of macros for given model and get list index of corresponding min macro
        for id, macro in enumerate(macroData.get(modelNum)):
                if macro.get("macroID") == minMacroID:
                    minMacroID = id

        minMacro = macroData.get(modelNum)[minMacroID]
        minPar   = minMacro.get("parSolution")

        print("\nMinimum solution: %s\n" %minMacro)

        minModelQ[0], minModelNR[0] = genModelNR(minPar,Q.get(modelNum))

    return minModelQ, minModelNR


def main():

    if config.verbose == False:
        print("Verbose = False, terminal silenced.")

    # prelim
    inputFiles  = ['../input/S14_excel.txt','../input/S13.txt','../input/S11_excel.txt',]
    inputLabels = ['h-MC3 D2O','d-MC3 D2O','d-MC3 ACMW',]

    d2_0, d2_1, d2_step = 5, 10.5, 1
    d1_0, d1_1, d1_step = 10, 15.5, 1

    # define which head group thicknesses are to be studied
    d2List = []
    for num in np.arange(d2_0, d2_1, d2_step):
        d2List.append(num)

    # define tail thicknesses for case fixedD1=True
    d1List = []
    for num in np.arange(d1_0, d1_1, d1_step):
        d1List.append(num)

    # define number of models
    nModels = len(inputFiles)

    # initialise experimental data dictionaries
    Q     = {init: [] for init in range(nModels)}
    expNR = {init: [] for init in range(nModels)}

    # for each model store list of genetic outputs
    costList  = {init: [] for init in range(nModels)}
    macroData = {init: [] for init in range(nModels)}

    # iterate over every file to get experiment data
    macroID = 0
    for modelNum, fileDIR in enumerate(inputFiles):
        print("\nmodel = %d/%d" %(modelNum+1,nModels))

        # get sample data as pandas df
        experimentFile = getFile(path=fileDIR, nSkip=0, delim='\t')

        # parse variables
        Q[modelNum]     = experimentFile[experimentFile.columns.values[0]]
        expNR[modelNum] = experimentFile[experimentFile.columns.values[1]]

        # for each file/model, fit d1 for given d2, N times
        for d2ID, d2 in enumerate(d2List):
            print("d2 = %d/%d" %(d2ID+1,len(d2List)))

            # matrix structure where for each d1 and d2, backingSLD is fitted
            if config.fixedD1 == True:
                for d1ID, d1 in enumerate(d1List):
                    print(" d1 = %d/%d" %(d1ID+1,len(d1List)))
                    # numnber of repeats
                    N = 1
                    for n in range(N):

                        # run fitting module
                        maxIter = 1
                        res = geneticAlgo(maxIter, Q.get(modelNum), expNR.get(modelNum), modelNum, (d2,d1))

                        # initialise geneticOutput dictionary
                        geneticOutput = {}

                        # macro ID number
                        geneticOutput["macroID"] = macroID
                        macroID += 1

                        # model number
                        geneticOutput["modelNum"] = modelNum

                        # head group thickness
                        geneticOutput["d2"] = d2

                        # tail thickness
                        geneticOutput["d1"] = d1

                        # parameter solution
                        geneticOutput["parSolution"] = res.x

                        # associated cost
                        expData = (Q.get(modelNum), expNR.get(modelNum))
                        geneticOutput["Cost"] = leastSquares(res.x, *expData)
                        costList[modelNum].append( (geneticOutput.get("Cost"), geneticOutput.get("macroID")))

                        # append dict to a list in a macro dict of each model
                        macroData[modelNum].append(geneticOutput)



            # else case is where d1 is left as a free parameter
            else:
                # numnber of repeats
                N = 1
                for n in range(N):

                    maxIter = 1
                    res = geneticAlgo(maxIter, Q.get(modelNum), expNR.get(modelNum), modelNum, d2)

                    geneticOutput = {}
                    geneticOutput["modelNum"] = modelNum
                    geneticOutput["d2"] = d2
                    geneticOutput["d1"] = res.x[7]
                    geneticOutput["parSolution"] = res.x

                    # associated cost
                    expData = (Q.get(modelNum), expNR.get(modelNum))
                    geneticOutput["Cost"] = leastSquares(res.x, *expData)
                    macroData[modelNum].append(geneticOutput)

        # plot colour map for every model
        plotColorMap(inputLabels[modelNum], macroData.get(modelNum), Q, d2List, d1List)

        nModelsToPlot, plotAllModels = 1, False
        minModelQ, minModelNR = getMinModelData(plotAllModels,nModelsToPlot,modelNum,costList,macroData,Q)
        plotFits(plotAllModels,modelNum,Q,expNR,minModelQ,minModelNR,inputLabels)

    # get model data for the lowest cost solutions and plot
    plotAllModels = True
    minModelQ, minModelNR = getMinModelData(plotAllModels,nModels,modelNum,costList,macroData,Q)
    plotFits(plotAllModels,modelNum,Q,expNR,minModelQ,minModelNR,inputLabels)


    #print('\n\nMacroData:\n%s\n\n' %macroData)
    # print output parameters, statistics and a figure
    #geneticAnalysis(geneticOutputData, Q, expNR, inputLabels)

    #plotCost(N, macroData, Q, expNR, inputLabels)

    # for fits where  d2 = [6, 6.1, 6.2 ... 15] and
    #plotThickness(N, macroData, Q, expNR, d2List, inputLabels)

    return




if __name__ == '__main__':
    print('--------------------------------------------------------------------')
    print('Program NRfits - Fitting Program for Neutron reflectivity')
    print('Version 0.0.2, April 2022')
    print('Developed by Samuel Winnall. @ UoM')
    print('--------------------------------------------------------------------\n\n')
    main()
