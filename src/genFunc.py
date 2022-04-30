" Module for general functions within NRfits "
# Next: add in geneneral smooth data functions

import csv
import sys
import pandas as pd


def modSelection(analysisOptions):

    # ask user to pick one of the analysisOptions
    print("\n~~~\nAnalysis Options:")
    for i,option in enumerate(analysisOptions):
        print("%d: %s" %(i+1,option))
    print("~~~\n")

    analysisChoice = input("Which analysis would you like to do? Pick the associated number (1-%d):\n  " %len(analysisOptions) )

    if analysisChoice.upper() == 'Q':
        print("Session closed.")
        sys.exit()

    elif analysisChoice in [str(i) for i in range(len(analysisOptions)+1)]:
        analysisType = analysisOptions[int(analysisChoice)-1]
        print("You picked %s.py\n" %analysisType)
        analysisRunning = True

    elif analysisChoice not in [str(i) for i in range(len(analysisOptions)+1)] and analysisOptions[0]!="isoAnalysis":
        print("Not a valid response. Returning to Pluto landing page.\n\n")
        analysisType    = 'n/a'
        analysisRunning = False

    else:
        print("Not a valid response. Session closed.")
        sys.exit()

    return analysisType, analysisRunning



# general function for Pluto, isotherm and surface excess modules
def getFile(path,nSkip,delim):
    return pd.read_csv(path, skiprows=nSkip, sep=delim, comment='#', na_values =' ', skip_blank_lines=True, encoding = "utf-8") # on_bad_lines='skip',
