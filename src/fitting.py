" Fit NR data "

import scipy.optimize as opt
from genModelData import genModelNR
import numpy as np
import config


# least square condition
def leastsquare(expNR, modelNR):

    ## need to update this function such that it iterates over all models
    ## and sums chiSq accordingly

    diffSq = []
    for [modelIdx] in len(modelNR):

        for i, j in zip(expNR[modelIdx],modelNR[modelIdx]):
            diffSq.append( np.power(i-j,2) )

    return np.sum(diffSq)


# residual function for genetic algorithm
def residuals(Q, expNR, modelList):

    # initialise dict everytime model data is calculated
    # data is stored for each model
    modelQ  = {init: [] for i in range(nModel)}
    modelNR = {init: [] for i in range(nModel)}
    for i in range(nModel):
        modelQ[i], modelNR[i] = genModelNR(modelList[i],Q[i])

    return leastsquare(expNR, modelNR)


def geneticAlgo(Q, expNR, modelList):

    # associated parameter bounds; could fix pars by defining in model function
    d1_lb  = 10
    d1_ub  = 20
    d2_lb  = 5
    d2_ub  = 10
    bounds = [(d1_lb,d1_ub),(d2_lb,d2_ub)]

    # define constraints; x[1] (d2) < d2_lb and > x[0] (d1)
    #print(par[1])
    #print(d2_lb)
    #print(par[0])
    #lc = opt.LinearConstraint(np.ones(1)*par[1], d2_lb, par[0])
    #print(lc)

    # genetic algorithm; might need args=*par or make par global
    geneticOutput = opt.differential_evolution(residuals, bounds, args=(Q, expNR, modelList), maxiter=1000)

    return geneticOutput
