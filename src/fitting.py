" Fit NR data "

import scipy.optimize as opt
from genModelData import genModelNR
import numpy as np
import config


# least square condition
def leastsquare(expNR, modelNR):

    diffSq = []
    for i, j in zip(expNR,modelNR):
        diffSq.append( np.power(i-j,2) )

    return np.sum(diffSq)


# residual function for genetic algorithm
def residuals(par, Q, expNR):

    modelQ, modelNR = genModelNR(par,Q)
    #print(len(Q))
    #print(len(expNR))
    #print(len(modelQ))
    #print(len(modelNR))

    return leastsquare(expNR, modelNR)


def geneticAlgo(Q, expNR):

    # define input parameters; d1, d2
    d1  = 12.4130
    d2  = 6.0
    par = [d1, d2]

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
    geneticOutput = opt.differential_evolution(residuals, bounds, args=(Q, expNR), maxiter=1000)

    return geneticOutput
