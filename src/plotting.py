from fitting import residuals
from genModelData import suppress_stdout, genModelNR
import matplotlib.pyplot as plt
import config


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

    # never interested in re-printing this output 
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
