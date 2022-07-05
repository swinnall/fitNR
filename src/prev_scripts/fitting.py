" Fit NR data "

import scipy.optimize as opt
from genModelData import genModelNR
from sldAnalysis import mainCalcSolvFrac
import numpy as np
import sys
import config


def getModel(i,thickPar):

    roughness = 3.5
    Im_SLD    = 0

    Re_SLD_Air = 0
    Thick_Air  = 0
    Rough_Air  = 0
    Solv_Air   = 0

    Re_SLD_Tails = [-0.0730e-6,6.1932e-6,6.1932e-6] # d-tails / h-tails
    Thick_Tails = thickPar[1]
    Solv_Tails   = 0

    Re_SLD_Heads = 0.7262e-6
    Thick_Heads  = thickPar[0]
    Solv_Heads   = 0 # init, calc@next step

    Re_SLD_Buffer = [6.1e-6,6.1e-6,0] # D2O~6.1, ACMW~0
    Thick_Buffer  = 0
    Solv_Buffer   = 0

    modelPar = [ Re_SLD_Air,          Im_SLD, Thick_Air,      Rough_Air, Solv_Air,
                    Re_SLD_Tails[i],  Im_SLD, Thick_Tails,    roughness, Solv_Tails,
                    Re_SLD_Heads,     Im_SLD, Thick_Heads,    roughness, Solv_Heads,
                    Re_SLD_Buffer[i], Im_SLD, Thick_Buffer,   roughness, Solv_Buffer
                    ]

    return modelPar


def getBounds(i,headBounds,solvHeadBounds,tailBounds):

    roughness = (3.5,3.5)
    Im_SLD    = (0,0)

    Re_SLD_Air = (0,0)
    Thick_Air  = (0,0)
    Rough_Air  = (0,0)
    Solv_Air   = (0,0)

    Re_SLD_Tails = [(-0.0730e-6,-0.0730e-6),(6.1932e-6,6.1932e-6),(6.1932e-6,6.1932e-6)] # d-tails / h-tails
    Thick_Tails  = tailBounds
    Solv_Tails   = (0,0)

    Re_SLD_Heads = (0.7262e-6,0.7262e-6)
    Thick_Heads  = headBounds
    Solv_Heads   = solvHeadBounds

    Re_SLD_Buffer = [(6.1e-6,6.1e-6),(6.1e-6,6.1e-6),(0,0)] # ACMW / D2O (0.0,0.1)
    Thick_Buffer  = (0,0)
    Solv_Buffer   = (0,0)


    bounds = [ Re_SLD_Air,       Im_SLD, Thick_Air,      Rough_Air, Solv_Air,
               Re_SLD_Tails[i],  Im_SLD, Thick_Tails,    roughness, Solv_Tails,
               Re_SLD_Heads,     Im_SLD, Thick_Heads,    roughness, Solv_Heads,
               Re_SLD_Buffer[i], Im_SLD, Thick_Buffer,   roughness, Solv_Buffer
              ]

    return bounds



# function to be minimised (residuals via least squared)
def leastSquares(modelPar, *expData):

    # unpack data
    Q, expNR = expData

    # calculate % solvent in monolayer and update modelPar
    modelPar[14] = mainCalcSolvFrac(modelIdx, modelPar[7],modelPar[12])

    # generate associated model data
    modelQ, modelNR = genModelNR(modelPar,Q)

    # sum over diffSq for each model
    diffSq = []
    for i, j in zip(expNR,modelNR):
        diffSq.append( np.power(i-j,2) )

    return np.sum(diffSq)



def geneticAlgo(maxIter, Q, expNR, modelNum, thickPar):

    # create global modelIdx variable to bypass limited ways of telling SLDAnalysis which model to study
    # modelNum can't be global as it's already in the local namespace
    global modelIdx
    modelIdx = modelNum

    # get model parameters
    modelPar = getModel(modelNum,thickPar)

    # update modelPar with correct solvent fraction
    modelPar[14] = mainCalcSolvFrac(modelIdx, modelPar[7], modelPar[12])

    # get parameter bounds
    headBounds     = (modelPar[12],modelPar[12])
    solvHeadBounds = (modelPar[14],modelPar[14])

    if config.fixedD1 == True:
        tailBounds = (modelPar[7],modelPar[7])
    else:
        print("Fatal Error: no tail bounds input.")
        sys.exit()

    bounds = getBounds(modelNum,headBounds,solvHeadBounds, tailBounds)

    # define constraints; x[7] (d1) > x[12] (d2) and < d1_ub
    #lc = opt.LinearConstraint(np.array(modelPar[7]), modelPar[12], bounds[7][1])
    #print(lc)

    # putting experimental data into args
    args = (Q,expNR)

    # genetic algorithm; might need args=*par or make par global constraints=(lc),
    geneticOutput = opt.differential_evolution(leastSquares, bounds, args, maxiter=maxIter)

    return geneticOutput
