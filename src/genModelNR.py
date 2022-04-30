" Generate NR model data "

from anaklasis import ref
import matplotlib.pyplot as plt
import config


# plot test figure
def genSimplePlot(res):

    # initialise figure
    fig = plt.figure()
    ax  = plt.axes()

    # fontsize
    fs = 10

    # Rmodel vs Q
    plt.plot(res1[("reflectivity")][:,0],res1[("reflectivity")][:,1])
    plt.yscale('log')

    # set axis labels
    ax.set_xlabel("$Q$ ($\AA^{-1}$)",fontsize=fs)
    ax.set_ylabel("$R$", color="black", fontsize=fs)

    # set y axis range
    #plt.ylim([0, 1.1])

    # legend
    ax.legend()

    # grid
    plt.grid(False)

    # tight layout function
    plt.tight_layout()

    # show plot
    plt.show()

    # save the plot as a file
    fig.savefig('output/R_modelTest.png',
            format='png',
            dpi=200,
            bbox_inches='tight')

    return





def main():

    input_file = '../input/S14.mft' # input curve
    units = ['A'] # Q units in Angstrom

    project='MC3 PBS ref test'
    #project='2layers'

    # We have a single uniform layer with full coverage
    patches=[1.0]

    # Create single model(patch) list
    model=[
     # Re_sld Im_sld thk rough solv description
    	[ 0.000e-6,   0.00e-6,    0, 0.0, 0.00, 'Air'],
    	[ -0.0730e-6, 0.00e-6, 12.4, 3.5, 0.00, 'tails'],
    	[ 0.7262e-6,  0.00e-6,  6.0, 3.5, 0.52, 'inner_heads'],
    # 	[ 9.41e-6, 0.00e-6, 20, 3.5, 1.0, 'nucleic_acid'],
    	[ 6.10e-6,    0.00e-6,    0, 3.5, 0.00, 'D2O'],
        ]


    system=[model]
    global_param = []
    resolution=[0.05]
    background = [5.0e-7]
    scale = [1.0]
    qmax = [0.3]

    # generate model data  
    res = ref.calculate(project, resolution, patches, system, global_param,
            background, scale, qmax, plot=True)

    # for comparing datasets
    #res = ref.compare(project, input_file, units, resolution, patches, system,
            #global_param,background, scale, qmax, experror=True, plot=True)

     #genSimplePlot(res)

    return res



if __name__ == '__main__':
    print("~Running genModelNR.py~")
    main()
