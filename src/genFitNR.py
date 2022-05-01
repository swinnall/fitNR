" Fit NR data "

from anaklasis import ref
import matplotlib.pyplot as plt
import scipy.optimize as opt
import numpy as np
import config
from genFunc import getFile

from contextlib import contextmanager
import sys, os

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

    # unpack parameters
    d1 = par[0]
    d2 = par[1]

    # project name ('none' = no save)
    project='none'

    # We have a single uniform layer with full coverage
    patches=[1.0]

    # Create single model(patch) list; Re_sld Im_sld thk rough solv description
    model=[
    	[ 0.000e-6,   0.00e-6,    0, 0.0, 0.00, 'Air'],
    	[ 6.19302e-6, 0.00e-6,   d1, 3.5, 0.00, 'tails'],
    	[ 0.7262e-6,  0.00e-6,   d2, 3.5, 0.52, 'inner_heads'],
    # 	[ 9.41e-6, 0.00e-6, 20, 3.5, 1.0, 'nucleic_acid'],
    	[ 0.00e-6,    0.00e-6,    0, 3.5, 0.00, 'ACMW'],
        ]

    # model=[ D2O
    # 	[ 0.000e-6,   0.00e-6,    0, 0.0, 0.00, 'Air'],
    # 	[ -0.0730e-6, 0.00e-6,   d1, 3.5, 0.00, 'tails'],
    # 	[ 0.7262e-6,  0.00e-6,   d2, 3.5, 0.52, 'inner_heads'],
    # # 	[ 9.41e-6, 0.00e-6, 20, 3.5, 1.0, 'nucleic_acid'],
    # 	[ 6.10e-6,    0.00e-6,    0, 3.5, 0.00, 'D2O'],
    #     ]

    system       = [model]
    global_param = []
    resolution   = [0.05]
    background   = [5.0e-7]
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

    # extract model data
    #qminIDX = np.where(res[("reflectivity")][:,0] == Q[0])
    modelQ  = res[("reflectivity")][:,0]
    modelNR = res[("reflectivity")][:,1]

    # close all figures to prevent runtime warning
    plt.close('all')

    return modelQ, modelNR


# least square condition
def leastsquare(expNR, modelNR):

    diffSq = []
    for i, j in zip(expNR,modelNR):
        diffSq.append( np.power(i-j,2) )

    return np.sum(diffSq)


# residual function for genetic algorithm
def residuals(par, Q, expNR):

    modelQ, modelNR = genModelNR(par,Q)
    print(len(Q))
    print(len(expNR))
    print(len(modelQ))
    print(len(modelNR))

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
    d1  = 12
    d2  = 6
    par = [d1, d2]

    # associated parameter bounds; could fix pars by defining in model function
    d1_lb  = 15
    d1_ub  = 25
    d2_lb  = 5
    d2_ub  = 15
    bounds = [(d1_lb,d1_ub),(d2_lb,d2_ub)]

    # x[1] (d2) < d2_lb and > x[0] (d1)
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
    print("\nCost of chosen solution: %f" %lstsq)

    # number of iterations
    print("\nNumber of iterations: %s" %geneticOutput.nit)

    # bool of success
    print("\nOptimisation status: %s" %geneticOutput.success)

    # termination message
    print("\nTermination message: %s" %geneticOutput.message)

    # generate figure
    fig, ax = plt.subplots()

    #for spine in ['top', 'right', 'bottom', 'left']:
    #    ax.spines[spine].set_linewidth(2)

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
    chiSqText = 'ChiSq = ' + "{:.3f}".format(lstsq) + ''
    props     = dict(boxstyle='none', facecolor='none', alpha=0.5)
    plt.text(1.0, 3.0, chiSqText, transform=ax.transAxes, fontsize=fs, fontweight='bold')

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
