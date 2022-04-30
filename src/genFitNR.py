" Generate NR model data "

from anaklasis import ref
import matplotlib.pyplot as plt
import scipy.optimize as opt
import numpy as np
import config
from genFunc import getFile


def genModelNR(par,Q):

    # unpack parameters
    d1 = par[0]
    d2 = par[1]

    input_file = '../input/S14.mft' # input curve
    units = ['A'] # Q units in Angstrom

    project='MC3 PBS ref test'

    # We have a single uniform layer with full coverage
    patches=[1.0]

    # Create single model(patch) list
    model=[
     # Re_sld Im_sld thk rough solv description
    	[ 0.000e-6,   0.00e-6,    0, 0.0, 0.00, 'Air'],
    	[ -0.0730e-6, 0.00e-6,   d1, 3.5, 0.00, 'tails'],
    	[ 0.7262e-6,  0.00e-6,   d2, 3.5, 0.52, 'inner_heads'],
    # 	[ 9.41e-6, 0.00e-6, 20, 3.5, 1.0, 'nucleic_acid'],
    	[ 6.10e-6,    0.00e-6,    0, 3.5, 0.00, 'D2O'],
        ]

    system = [model]
    global_param = []
    resolution = [0.05]
    background = [5.0e-7]
    scale = [1.0]
    qmax = [Q[len(Q)-1]]

    # generate model data
    res = ref.calculate(project, resolution, patches, system, global_param,
            background, scale, qmax, plot=False)

    modelNR = res[("reflectivity")][:,1]


    return modelNR


# least square condition
def leastsquare(expNR, modelNR):

    diffSq = []
    for i, j in zip(expNR,modelNR):
        diffSq.append( np.power(i-j,2) )

    return np.sum(diffSq)


# residual function for genetic algorithm
def residuals(par, Q, expNR):

    # where sinModel is the sinusoidal function
    modelNR = genModelNR(par,Q)

    return leastsquare(expNR, modelNR)


def geneticAlgo():

    # read file into memory data
    fileDIR = '../input/S14_excel.txt'

    # get sample data as pandas df
    data = getFile(path=fileDIR, nSkip=0, delim='\t')

    # parse variables
    Q     = data[data.columns.values[0]]
    expNR = data[data.columns.values[1]]

    # define input parameters; d1, d2
    par = [12.0, 6.0]

    # associated parameter bounds; could fix pars by defining in model function
    bounds = [(10,20),(5,15)]

    # genetic algorithm; might need args=*par or make par global
    geneticOutput = opt.differential_evolution(residuals, bounds, args=(Q, expNR), maxiter=10)

    return geneticOutput, Q, expNR


def printGeneticOutput(geneticOutput, Q, expNR):

    # parameter solution
    solution = geneticOutput.x
    print("\nParameter solution (d1, d2):\n %s" %solution)

    # associated cost
    lstsq = residuals(solution, Q, expNR)
    print("\nCost of chosen solution: %f" %lstsq)

    # number of iterations
    print("\nNumber of iterations: %s" %geneticOutput.nit)

    # bool of success
    print("\nOptimisation status: %s" %geneticOutput.success)

    # termination message
    print("\nTermination message: %s" %geneticOutput.message)

    return


def main():

    # fit the model to the experimental data
    geneticOutput, Q, expNR = geneticAlgo()

    # print output parameters and associated statistics
    printGeneticOutput(geneticOutput, Q, expNR)

    return



if __name__ == '__main__':
    print("~Running genFitNR.py~")
    main()
