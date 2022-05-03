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


def genModelNR(model,Q):

    # project name ('none' = no save)
    project='none'

    # We have a single uniform layer with full coverage
    patches=[1.0]

    # calculation parameters
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
