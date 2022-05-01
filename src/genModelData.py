from anaklasis import ref
from contextlib import contextmanager
import matplotlib.pyplot as plt
import numpy as np
import sys, os, config

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
