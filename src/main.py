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


def main():

    print("hi")

    # read file into memory data
    fileDIR = '../input/S11_excel.txt'

    # get sample data as pandas df
    experimentFile = getFile(path=fileDIR, nSkip=0, delim='\t')

    # parse variables
    Q     = experimentFile[experimentFile.columns.values[0]]
    expNR = experimentFile[experimentFile.columns.values[1]]

    # running fitting module
    geneticOutput = geneticAlgo(Q, expNR)

    # print output parameters, statistics and a figure
    geneticAnalysis(geneticOutput, Q, expNR)

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
