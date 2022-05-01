" Fit NR data "

from anaklasis import ref
from contextlib import contextmanager
import matplotlib.pyplot as plt
import scipy.optimize as opt
import numpy as np
import sys, os
import config
from genFunc import getFile


@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout


def genModelNR(par,Q):

    # project name ('none' = no save)
    project='none'

    # We have a single uniform layer with full coverage
    patches=[1.0]

    # unpack parameters
    d1 = par[0]
    d2 = par[1]

    # model parameters; col = air, tails, heads, ACMW
    Re_SLD = [0.000e-6, 6.1932e-6, 0.7262e-6, 0.00e-6]
    Im_SLD = [0, 0, 0, 0,]
    Thick  = [0, d1, d2, 0]
    Rough  = [0, 3.5, 3.5, 3.5]
    Solv   = [0, 0, 0.522, 0]

    # Create single model(patch) list
    model=[
    	[ Re_SLD[0], Im_SLD[0], Thick[0], Rough[0], Solv[0], 'air'    ],
    	[ Re_SLD[1], Im_SLD[1], Thick[1], Rough[1], Solv[1], 'd-tails'],
    	[ Re_SLD[2], Im_SLD[2], Thick[2], Rough[2], Solv[2], 'heads'  ],
    	[ Re_SLD[3], Im_SLD[3], Thick[3], Rough[3], Solv[3], 'ACMW'   ],
        ]
    # 	[ 0.7262e-6,  0.00e-6,   d2, 3.5, 0.52, 'h-tails'],
    # 	[ 9.41e-6, 0.00e-6, 20, 3.5, 1.0, 'nucleic acid'],
    # 	[ 6.10e-6,    0.00e-6,    0, 3.5, 0.00, 'D2O'],

    system       = [model]
    global_param = []
    resolution   = [0.05]
    background   = [4.5e-7]
    scale        = [1.0]
    qmaxIDX      = len(Q)-1
    qmax         = [Q[qmaxIDX]]

    # generate model data
    if config.verbose == True:
        res = ref.calculate(project, resolution, patches, system, global_param,
            background, scale, qmax, plot=False)
    else:
        with suppress_stdout():
            res = ref.calculate(project, resolution, patches, system, global_param,
                background, scale, qmax, plot=False)

    # close all figures to prevent runtime warning
    plt.close('all')

    modelQ  = res[("reflectivity")][:,0]
    modelNR = res[("reflectivity")][:,1]

    # find the closest value in the model array
    def closest_value(input_list, input_value):
        arr = np.asarray(input_list)
        idx = (np.abs(arr - input_value)).argmin()
        return arr[idx], idx

    reducedModelQ  = []
    reducedModelNR = []
    for qExpVal in Q:
        closeModelqVal, closeModelqValIDX = closest_value(modelQ,qExpVal)
        reducedModelQ.append(closeModelqVal)
        reducedModelNR.append(modelNR[closeModelqValIDX])

    return reducedModelQ, reducedModelNR


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


def geneticAlgo():

    # read file into memory data
    fileDIR = '../input/S11_excel.txt'

    # get sample data as pandas df
    data = getFile(path=fileDIR, nSkip=0, delim='\t')

    # parse variables
    Q     = data[data.columns.values[0]]
    expNR = data[data.columns.values[1]]

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

    return geneticOutput, Q, expNR



def geneticAnalysis(geneticOutput, Q, expNR):

    # parameter solution
    solution = geneticOutput.x
    print("\nParameter solution (d1, d2):\n %s" %solution)

    # associated cost
    lstsq = residuals(solution, Q, expNR)
    print("\nCost of chosen solution: %.8e" %lstsq)

    # number of iterations
    print("\nNumber of iterations: %s" %geneticOutput.nit)

    # bool of success
    print("\nOptimisation status: %s" %geneticOutput.success)

    # termination message
    print("\nTermination message: %s" %geneticOutput.message)

    # generate figure
    fig, ax = plt.subplots()

    with suppress_stdout():
        modelQ, modelNR = genModelNR(solution,Q)

    # fontsize
    fs = 14

    # Rmodel vs Q
    plt.plot(Q, expNR, 'o', label='Experiment')
    plt.plot(modelQ, modelNR, '-', label='Model')

    plt.yscale('log')

    # set axis labels
    plt.xlabel("Q ($\AA^{-1}$)", fontsize=fs, fontweight='bold')
    plt.ylabel("R", fontsize=fs, fontweight='bold')

    plt.tick_params(axis='x', labelsize=fs, which='major', size=5, width=1, direction='in', top='on')
    plt.tick_params(axis='y', labelsize=fs, which='major', size=5, width=1, direction='in', right='on')
    plt.tick_params(axis='y', labelsize=fs, which='minor', size=5, width=1, direction='in', right='on')

    # legend
    plt.legend(prop={'size': fs, 'weight':'bold'}, frameon = False, loc='upper right')

    # grid
    plt.grid(False)

    # chiSq annotation
    chiSqText = 'ChiSq = ' + "{:.4e}".format(lstsq) + ''
    props     = dict(boxstyle='none', facecolor='none', alpha=0.5)
    plt.text(0.80, 3.0, chiSqText, transform=ax.transAxes, fontsize=fs, fontweight='bold')

    # tight layout function
    plt.tight_layout()
    fig.subplots_adjust(wspace=0.05, hspace=0.05)

    # save the plot as a file
    plt.savefig('../output/NRfit.png',
            format='png',
            dpi=300,
            bbox_inches='tight')

    # show plot
    #plt.show()

    return



def main():

    # fit the model to the experimental data
    geneticOutput, Q, expNR = geneticAlgo()

    # print output parameters, statistics and a figure
    geneticAnalysis(geneticOutput, Q, expNR)

    return



if __name__ == '__main__':
    print("~Running genFitNR.py~")
    main()
