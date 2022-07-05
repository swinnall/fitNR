from fitting import leastSquares, getModel
from genModelData import suppress_stdout, genModelNR
import matplotlib.pyplot as plt
import config
import numpy as np


def plotColorMap(title, macroData, Q, d2List, d1List):

    nModels, nRows, nCols = len(Q), len(d1List), len(d2List)
    costMat = [[0 for x in range(nCols)] for y in range(nRows)]

    for macro in macroData:
        d2IDX = d2List.index(macro.get("d2"))
        d1IDX = d1List.index(macro.get("d1"))
        costMat[d1IDX][d2IDX] = macro.get("Cost")

    #print("\ncostMat = %s" %costMat)

    # generate figure
    fig, ax = plt.subplots()

    # extentfloats = (left, right, bottom, top)
    plt.imshow(costMat, origin='upper', extent=[d2List[0],d2List[-1],d1List[-1],d1List[0]]) #, vmin=0, vmax=1
    cb = plt.colorbar()

    # fontsize
    fs = 14

    # set axis labels
    plt.xlabel("d2 - heads ($\AA$)", fontsize=fs, fontweight='bold')
    plt.ylabel("d1 - tails ($\AA$)", fontsize=fs, fontweight='bold')
    cb.set_label(label='Cost', fontsize=fs, fontweight='bold', labelpad=10)

    # title
    plt.title(title, fontsize=fs, fontweight='bold')

    # save the plot as a file
    plt.savefig('../output/NRfitColorMap - '+ title +'.png',
            format='png',
            dpi=300,
            bbox_inches='tight')

    return


def plotFits(plotAllModels,modelNum,Q,expNR,modelQ,modelNR,inputLabels):

    # generate figure
    fig, ax = plt.subplots()

    # fontsize
    fs = 14

    # Rmodel vs Q
    if plotAllModels == True:
        nModels = len(Q)
        for i in range(nModels):
            plt.plot(Q.get(i), expNR.get(i), 'o', label=inputLabels[i])
            plt.plot(modelQ.get(i), modelNR.get(i), '-')

    else:
        plt.plot(Q.get(modelNum), expNR.get(modelNum), 'o', label=inputLabels[modelNum])
        plt.plot(modelQ.get(0), modelNR.get(0), '-')

    # set y scale
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

    # chiSq annotation - error: costList not defined
    #chiSqText = '\u03A3\u03C7$^{2}$ = ' + "{:.4e}".format(sum(costList)) + ''
    #plt.text(0.1, 0.90, chiSqText, transform=ax.transAxes, fontsize=fs, fontweight='bold')

    # tight layout function
    plt.tight_layout()
    fig.subplots_adjust(wspace=0.05, hspace=0.05)

    # save the plot as a file
    if plotAllModels == True:
        title = 'all models'
    else:
        title = inputLabels[modelNum]

    plt.savefig('../output/NRfit - ' + title + '.png',
            format='png',
            dpi=300,
            bbox_inches='tight')

    # show plot
    #plt.show()

    return



def geneticAnalysis(geneticOutputData, Q, expNR, inputLabels):

    nModels = len(Q)

    modelQ  = {init: [] for init in range(nModels)}
    modelNR = {init: [] for init in range(nModels)}

    costList = []
    for idx, result in enumerate(geneticOutputData):

        ### this no longer works as geneticOutputData is a dict of lists not a list

        # parameter solution
        solution = result.x
        print("\nParameter solution (d1, d2):\n %s" %solution)

        # associated cost
        lstsq = leastSquares(solution, Q.get(idx), expNR.get(idx))
        costList.append(lstsq)
        print("\nCost of chosen solution: %.8e" %lstsq)

        # number of iterations
        print("\nNumber of iterations: %s" %result.nit)

        # bool of success
        print("\nOptimisation status: %s" %result.success)

        # termination message
        print("\nTermination message: %s" %result.message)

        # generate associated model data and store for plotting
        # never interested in re-printing this output
        with suppress_stdout():
            modelQ[idx], modelNR[idx] = genModelNR(solution,Q.get(idx))


    print("I\nndividual Model Costs: %s" %costList)

    return


def plotCost(N, macroData, Q, expNR, inputLabels):

    nModels = len(Q)

    # generate figure
    fig, ax = plt.subplots()

    # fontsize
    fs = 14

    #
    d1   = {init: [] for init in range(nModels)}
    d2   = {init: [] for init in range(nModels)}
    cost = {init: [] for init in range(nModels)}


    for modelNum in range(nModels):

        # calculate cost of each stored parameter solution
        for macro in macroData.get(modelNum):

            #d1[modelNum].append(macro.get("d1"))
            d2[modelNum].append(macro.get("d2"))
            cost[modelNum].append(macro.get("Cost"))


        plt.plot(d2.get(modelNum), cost.get(modelNum), linestyle='--', marker='o', label=inputLabels[modelNum])


    # set axis labels
    plt.xlabel("d2 - heads (A)", fontsize=fs, fontweight='bold')
    plt.ylabel("Cost", fontsize=fs, fontweight='bold')

    plt.tick_params(axis='x', labelsize=fs, which='major', size=5, width=1, direction='in', top='on')
    plt.tick_params(axis='y', labelsize=fs, which='major', size=5, width=1, direction='in', right='on')
    plt.tick_params(axis='y', labelsize=fs, which='minor', size=5, width=1, direction='in', right='on')

    # legend
    plt.legend(prop={'size': fs, 'weight':'bold'}, frameon = False, loc='upper right')

    # grid
    plt.grid(False)

    # tight layout function
    plt.tight_layout()
    fig.subplots_adjust(wspace=0.05, hspace=0.05)

    # save the plot as a file
    plt.savefig('../output/NRfitCost-.png',
            format='png',
            dpi=300,
            bbox_inches='tight')

    return



def plotThickness(N, macroData, Q, expNR, headParList, inputLabels):

    nModels = len(Q)

    # generate figure
    fig, ax = plt.subplots()

    # fontsize
    fs = 14

    #
    #print(macroData)

    d1 = {init: [] for init in range(nModels)}
    for modelNum in range(nModels):

        # calculate cost of each stored parameter solution
        for macro in macroData.get(modelNum):

            #print(macro.get("d1"))
            d1[modelNum].append(macro.get("d1"))


        plt.plot(headParList, d1.get(modelNum), linestyle='--', marker='o', label=inputLabels[modelNum])


    # set axis labels
    plt.xlabel("d2 - heads (A)", fontsize=fs, fontweight='bold')
    plt.ylabel("d1 - tails (A)", fontsize=fs, fontweight='bold')

    plt.tick_params(axis='x', labelsize=fs, which='major', size=5, width=1, direction='in', top='on')
    plt.tick_params(axis='y', labelsize=fs, which='major', size=5, width=1, direction='in', right='on')
    plt.tick_params(axis='y', labelsize=fs, which='minor', size=5, width=1, direction='in', right='on')

    # legend
    plt.legend(prop={'size': fs, 'weight':'bold'}, frameon = False, loc='upper right')

    # grid
    plt.grid(False)

    # tight layout function
    plt.tight_layout()
    fig.subplots_adjust(wspace=0.05, hspace=0.05)

    # save the plot as a file
    plt.savefig('../output/NRfitThick-.png',
            format='png',
            dpi=300,
            bbox_inches='tight')

    return
