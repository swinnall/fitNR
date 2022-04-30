" Data Analysis Pipeline for the NR reflectometry "
" Author: @S.Winnall "

import glob, os, sys
import pandas as pd
import csv
import shutil
from shutil import copyfile
import config
import isoAnalysis
import ellipsAnalysis
import chemFormulations
import sldAnalysis
import surfaceExcess
from genFunc import modSelection, getFile


def programSelection(analysisOptions):

    # root input directory
    inputDir  = config.inputDir

    # root output directory
    outputDir = config.outputDir

    # ask user to pick one of the analysisOptions
    analysisType, analysisRunning = modSelection(analysisOptions)

    # calls from config database
    instructionsName = config.pathNames.get(analysisType)[0]

    # input/output folder names to be accessed within root in/output directory
    inputOutputStem  = config.pathNames.get(analysisType)[1]

    # instructions file path
    instructionsPath = "" + inputDir + instructionsName + ".txt"

    # input data files path
    inputDataPath    = '' + inputDir + '/' + inputOutputStem + ''

    # output data files path
    outputDataPath   = '' + outputDir + '/' + inputOutputStem + ''

    return analysisType, instructionsName, instructionsPath, inputDataPath, outputDataPath



def organisePaths(analysisType, instructionsName, instructionsPath, outputDataPath):

    # get analysis title
    with open(instructionsPath, newline = '') as f:
        title = list(csv.reader(f))[0][0].split('=')[1]


    # read instructions file with pandas
    instructionsFile = getFile(path=instructionsPath,nSkip=1,delim=',')


    # update output path folder dir to include title of chosen analysis
    # e.g. output/chemFormulations/analysisTitle/
    outputDataPath = '' + outputDataPath + '/' + title + ''


    try:
        # delete folder if exists and create it again
        shutil.rmtree(outputDataPath)
        os.mkdir(outputDataPath)
    except FileNotFoundError:
        os.mkdir(outputDataPath)


    try:
        # copy input instructions instructionsFilermation to outputDataPath directory
        shutil.copyfile(instructionsPath, outputDataPath + '/' + instructionsName + '.txt')

        # rename copied instructions file to include analysis title
        # this is the file which will be appended in the relevant analysis
        old_name              = outputDataPath + '/' + instructionsName + '.txt'
        outputInstructionFile = outputDataPath + '/' + title + ' - ' + instructionsName + '.txt'
        os.rename(old_name, outputInstructionFile)


    # if source and destination are same
    except shutil.SameFileError:
        print("Instructions Copying Error:\n  Source and destination represents the same file.")

    # if destination is a directory
    except IsADirectoryError:
        print("Instructions Copying Error:\n  Destination is a directory.")

    # if there is any permission issue
    except PermissionError:
        print("Instructions Copying Error:\n  Permission denied.")

    # other errors
    except:
        print("Instructions Copying Error:\n  Error occurred while copying file.")


    # clean output file by getting file (auto skips commented lines)
    File = getFile(path=outputInstructionFile,nSkip=0,delim=',')

    # write new file without the comments (isolates just the files that were run)
    File.to_csv(outputInstructionFile)

    return instructionsFile, title, outputDataPath, outputInstructionFile




def main():

    analysisOptions = ['fitNR']

    mainRunning = True
    while mainRunning:

        # selects module and creates relevant file paths
        analysisType, instructionsName, instructionsPath, inputDataPath, outputDataPath = programSelection(analysisOptions)

        # reads instructions instructionsFile, gets title and paths
        instructionsFile, title, outputDataPath, outputInstructionFile = organisePaths(analysisType, instructionsName, instructionsPath, outputDataPath)

        # calls analysis module
        if analysisType == analysisOptions[0]:
            NRanalysis.main(instructionsFile, title, inputDataPath, outputDataPath)

    return


if __name__ == '__main__':
	print('--------------------------------------------------------------------')
	print('Program NRfits - Fit Module for Neutron reflection datasets')
	print('version 0.0.1, April 2022')
	print('developed by Samuel Winnall. @ UoM')
	print('--------------------------------------------------------------------\n\n')

    main()
