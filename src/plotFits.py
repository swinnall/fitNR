" Plotting code for refnx solutions "

import matplotlib.pyplot as plt
import config 

def plotNR(Q,expNR,expNR_err,modelQ,modelNR,inputLabels,title):

    # fontsize
    fs = 15
    
    # border width 
    bw = 3    

    # generate figure
    fig, ax = plt.subplots()

    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(bw)

    # Rmodel vs Q
    nModels = len(Q)

    for i in range(nModels):
        plt.errorbar(x=Q.get(i), y=expNR.get(i), yerr=expNR_err.get(i), fmt='', marker='o', color=config.c[0][i], ms=8, mec=config.c[0][i], mfc=config.c[0][i], ecolor=config.c[0][i], elinewidth=4, lw=0, label=inputLabels.get(i), zorder=1)
        plt.plot(modelQ.get(i), modelNR.get(i), '-', color=config.c1[0][i], lw=4, zorder=2)


    # set y scale
    plt.yscale('log')

    # set axis labels
    plt.xlabel("Q ($\AA^{-1}$)", fontsize=fs, fontweight='bold')
    plt.ylabel("R", fontsize=fs, fontweight='bold')

    plt.tick_params(axis='x', labelsize=fs, which='major', size=6, width=bw, direction='in', top='on')
    plt.tick_params(axis='y', labelsize=fs, which='major', size=6, width=bw, direction='in', right='on')
    plt.tick_params(axis='y', labelsize=fs, which='minor', size=6, width=bw, direction='in', right='on')

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

    plt.savefig('../output/global-objective-'+title+'.png',
            format='png',
            dpi=300,
            bbox_inches='tight')

    # show plot
    plt.show()

    return