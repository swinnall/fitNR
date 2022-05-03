" Data Analysis Pipeline for the NR reflectometry "
" Author: @S.Winnall "

import glob, os, sys, csv, shutil
from shutil import copyfile
import pandas as pd
import config
from genModelData import genModelNR
from fitting import geneticAlgo
from plotting import geneticAnalysis


def getFile(path,nSkip,delim):
    return pd.read_csv(path, skiprows=nSkip, sep=delim, comment='#', na_values =' ', skip_blank_lines=True, encoding = "utf-8") # on_bad_lines='skip',


def getModel(Re_SLD,Im_SLD,Thick,Rough,Solv):

    # Create single model(patch) list
    model0=[
    	[ Re_SLD[0], Im_SLD[0], Thick[0], Rough[0], Solv[0], 'air'    ],
    	[ Re_SLD[1], Im_SLD[1], Thick[1], Rough[1], Solv[1], 'd-tails'],
    	[ Re_SLD[2], Im_SLD[2], Thick[2], Rough[2], Solv[2], 'heads'  ],
    	[ Re_SLD[3], Im_SLD[3], Thick[3], Rough[3], Solv[3], 'ACMW'   ],
        ]
    # 	[ 0.7262e-6,  0.00e-6,   d2, 3.5, 0.52, 'h-tails'],
    # 	[ 9.41e-6, 0.00e-6, 20, 3.5, 1.0, 'nucleic acid'],
    # 	[ 6.10e-6,    0.00e-6,    0, 3.5, 0.00, 'D2O'],

    return model

def main():

    # define input parameters; d1, d2
    d1  = 12.4130
    d2  = 6.0

    # model parameters; col = air, tails, heads, solvent
    # each row in the list of lists is a separate model
    Re_SLD = [[0.000e-6, 6.1932e-6, 0.7262e-6, 0.00e-6],[0.000e-6, -0.0730, 0.7262e-6, 0.00e-6]]
    Im_SLD = [[0, 0, 0, 0,],[0, 0, 0, 0,]]
    Thick  = [[0, d1, d2, 0],[0, d1, d2, 0]]
    Rough  = [[0, 3.5, 3.5, 3.5],[0, 3.5, 3.5, 3.5]]
    Solv   = [[0, 0, 0.522, 0],[0, 0, 0.522, 0]]

    # read file into memory data
    inputFiles = ['../input/S11_excel.txt','../input/S14_excel.txt']
    nFiles = len(inputFiles)

    # initialise dicts
    Q     = {init: [] for i in range(nFiles)}
    expNR = {init: [] for i in range(nFiles)}

    # iterate over every file to get experiment data
    for idx, fileDIR in enumerate(inputFiles):

        # get sample data as pandas df
        experimentFile = getFile(path=fileDIR, nSkip=0, delim='\t')

        # parse variables
        Q[idx]     = experimentFile[experimentFile.columns.values[0]]
        expNR[idx] = experimentFile[experimentFile.columns.values[1]]

        # get model
        modelList[idx] = getModel(Re_SLD[idx],Im_SLD[idx],Thick[idx],Rough[idx],Solv[idx])

    # running fitting module
    geneticOutput = geneticAlgo(Q, expNR, modelList)

    # print output parameters, statistics and a figure
    # geneticAnalysis(geneticOutput, Q, expNR)

    ### Next Steps ###
    # Extend to fitting two datasets at once
    # Read in two files, add another function in residuals before generating model data
    # that gets the relevant (SLD) parameters (eventually from txt input file)
    # update genModelData such that it takes the full model matrix AND which pars to be fitted
    # store each model in instructions as dict and calc leastSquares (need to update)
    # remember pars (d1,d2) must be able to be varied according to the bounds
    # e.g. pars come from geneticAlgo and values come from user

    # alternative is to fit one at a time, sum chiSq but I don't think this works long term
    # would need extra steps to validate the parameters of one fit against the others etc

    # then implement constraints
    # then start batch fit protocol, refine workflow, then incorporate SLDanalysis

    return




if __name__ == '__main__':
    print('--------------------------------------------------------------------')
    print('Program NRfits - Fit Module for Neutron reflection datasets')
    print('Version 0.0.1, April 2022')
    print('Developed by Samuel Winnall. @ UoM')
    print('--------------------------------------------------------------------\n\n')
    main()
