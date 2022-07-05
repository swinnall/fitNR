" Calculate parameters needed for refnx "

import csv
import glob, os
import re
import sys
import config

class Membrane:

    def __init__(self, lipidNames, molRatios, thickness, monolayerPar):

        # read system parameters
        self.lipidNames  = lipidNames
        self.molRatios   = molRatios
        self.thickness   = thickness
        self.nLipids     = len(self.lipidNames)

        # unpack monolayer parameters (0s if 1st run; each tuple ele is dict[struct])
        self.monolayerMolVol = monolayerPar[0]
        self.monolayerSL     = monolayerPar[1]
        self.monolayerSLD    = monolayerPar[2]

        # import filepath
        #self.outputFilePath = outputFilePath

        # import databases
        self.lipidStruct = config.lipidStruct
        self.atomSL      = config.atomSL
        self.lipidMolVol = config.lipidMolVol

        # initialise new variables
        self.headVolFrac        = 0
        self.twoSolv            = 0
        self.d3                 = 0
        self.threeSolv          = 0
        self.twoSLD_H2O         = 0
        self.twoSLD_D2O         = 0
        self.drugOverlapSLD_H2O = 0
        self.drugOverlapSLD_D2O = 0
        self.normMolRatios      = []
        self.volFrac            = {}
        self.lipidSL            = {}
        self.sumAvSL            = {}
        self.sumAvSLD           = {}
        self.totalLipidVol      = {'head': 0, 'tails': 0}
        self.avLipidVol         = {'head': 0, 'tails': 0}
        self.avSL               = {'head': 0, 'tails': 0}
        self.avSLD              = {'head': 0, 'tails': 0}


    # converts molar ratio txt input to a normalised version e.g. 4:5 ==> 4/9:5/9
    def normaliseMolarRatios(self):

        totMol = 0
        for i, value in enumerate(self.molRatios):
            totMol += float(value)

        for i, value in enumerate(self.molRatios):
            self.normMolRatios.append( float(value) / totMol  )

        if config.verbose == True and "Monolayer" not in self.lipidNames:
            print("\nLipid names:\n%s\n\nInput molar ratios:\n%s" %(self.lipidNames, self.molRatios))

        if config.very_verbose == True:
            print("\nNormalised molar ratios:\n%s" %self.normMolRatios)


    # calculates the total lipid volume
    def calcTotalLipidVol(self):

        for i, lipid in enumerate(self.lipidNames):

            # check lipid exists in database
            if self.lipidStruct.get(lipid) == None and lipid == "Monolayer": pass
            elif self.lipidStruct.get(lipid) == None:
                print("\nError: Lipid type not found in Lipid Molecular Formula Database.")
                print("Lipid: %s" %lipid)
                sys.exit()

            # calculate total lipid volumes
            for j, struct in enumerate(['head','tails']):

                if lipid == "Monolayer":
                    self.totalLipidVol[struct] += self.normMolRatios[i] * self.monolayerMolVol[struct]
                else:
                    self.totalLipidVol[struct] += self.normMolRatios[i] * self.lipidMolVol.get(lipid)[j]


        # save Monolayer monolayer struct volumes on the first iteration
        if "Monolayer" not in self.lipidNames:
            self.monolayerMolVol['head']  = self.totalLipidVol.get('head')
            self.monolayerMolVol['tails'] = self.totalLipidVol.get('tails')

        if config.very_verbose == True:
            print('\nTotal Volume:\n%s' %self.totalLipidVol)



    # calculates volume fraction
    def calcVolFrac(self):

        warningChoice = input("\nWarning: You have set volFrac = True. Do you want to continue? (y/n)\n ")

        if warningChoice.upper() == 'N':
            print("Session closed, you must change config parameter.")
            sys.exit()


        # component volume calculation
        for i, lipid in enumerate(self.lipidNames):
            self.volFrac[lipid] = {'head': 0, 'tails': 0}

            for j, struct in enumerate(['head','tails']):
                if lipid == "Monolayer":
                    self.volFrac[lipid][struct] = self.normMolRatios[i] * self.monolayerMolVol.get(struct) / self.totalLipidVol[struct]
                else:
                    self.volFrac[lipid][struct] = self.normMolRatios[i] * self.lipidMolVol.get(lipid)[j] / self.totalLipidVol[struct]

        if config.verbose == True:
            print('\nComponent Volume Fraction:\n%s' %self.volFrac)



    # calculates the scattering length of each lipid component
    def calcSL(self):

        for i, lipid in enumerate(self.lipidNames):
            self.lipidSL[lipid]       = {'head': 0, 'tails': 0}

            for j, struct in enumerate(['head','tails']):

                if lipid == "Monolayer":
                    self.lipidSL[lipid][struct] = self.monolayerSL.get(struct)

                else:
                    # splits head/tail into list of constituent atoms
                    splitStruct = re.split('-', self.lipidStruct.get(lipid)[j])

                    # this loop iterates across elements in a given head/tail and sums the atomic scattering lengths
                    for ele in splitStruct:

                        # multiply scattering length of atom by number of atoms
                        if hasNumbers(ele) == True:

                            # x is split into ['atom','number of given atom']
                            x = list(filter(None, re.split(r'(\d+)', ele)))

                            self.lipidSL[lipid][struct] += self.atomSL.get(x[0]) * int(x[1])


                        # add the scattering length of the single atom identified
                        elif hasNumbers(ele) == False:
                            self.lipidSL[lipid][struct] += self.atomSL.get(ele[0])


                # multiply total lipid scattering length of a given lipid's head/tail by corresponding vol frac
                if config.useVolFrac == True:
                    self.avSL[struct] += self.volFrac[lipid][struct] * self.lipidSL[lipid][struct]
                else:
                    self.avSL[struct] += self.normMolRatios[i] * self.lipidSL[lipid][struct]


        # calculate the average SL of the whole averaged lipid
        self.sumAvSL = self.avSL.get("head") + self.avSL.get("tails")

        # save Monolayer monolayer struct volumes on the first iteration
        if "Monolayer" not in self.lipidNames:
            self.monolayerSL['head']  = self.avSL.get('head')
            self.monolayerSL['tails'] = self.avSL.get('tails')

        if config.very_verbose == True:
            print("\nLipid scattering lengths:\n%s" %self.lipidSL)

        if config.verbose == True:
            print("\nAverage SL:\n%s" %self.avSL)

        if config.very_verbose == True:
            print("\nSummed average SL:\n%s" %self.sumAvSL)


    def calcAvLipidVol(self):

        for i, lipid in enumerate(self.lipidNames):

            # calculate total lipid volumes
            for j, struct in enumerate(['head','tails']):

                if config.useVolFrac == True:
                    if lipid == "Monolayer":
                        self.avLipidVol[struct] += self.volFrac[lipid][struct] * self.monolayerMolVol[struct]
                    else:
                        self.avLipidVol[struct] += self.volFrac[lipid][struct] * self.lipidMolVol.get(lipid)[j]

                else:
                    if lipid == "Monolayer":
                        self.avLipidVol[struct] += self.normMolRatios[i] * self.monolayerMolVol[struct]
                    else:
                        self.avLipidVol[struct] += self.normMolRatios[i] * self.lipidMolVol.get(lipid)[j]


        # if accounting for chain compaction
        if config.compactChains == True:
            chainCompactFactor = config.chainCompactFactor
            self.avLipidVol['tails'] = chainCompactFactor * self.avLipidVol['tails']


        if config.verbose == True:
            print("\nAverage Lipid Head/Tail Volume:\n%s" %self.avLipidVol)



    # calculates the scattering length density of each lipid component
    def calcSLD(self):

        # take avSL and divide by avStructVol (tails, heads)
        for i, struct in enumerate(['head','tails']):
            self.avSLD[struct] = 10 * self.avSL[struct] / self.avLipidVol[struct]

        # calculate the average SLD of the whole averaged lipid
        self.sumAvSLD = self.avSLD.get("head") + self.avSLD.get("tails")

        # save Monolayer monolayer struct volumes on the first iteration
        if "Monolayer" not in self.lipidNames:
            self.monolayerSLD['head']  = self.avSLD.get('head')
            self.monolayerSLD['tails'] = self.avSLD.get('tails')

        if config.verbose == True:
            print("\nAverage SLD:\n%s" %self.avSLD)

        if config.very_verbose == True:
            print("\nSummed average SLD:\n%s" %self.sumAvSLD)


    # calculates volume fraction of the head group based on SLD
    def calcHeadVolumeFraction(self):

        # call SL
        SL1 = self.avSL.get("tails")
        SL2 = self.avSL.get("head")

        # call SLD
        rho1 = self.avSLD.get("tails")
        rho2 = self.avSLD.get("head")

        # call membrane thickness
        if "Monolayer" not in self.lipidNames:
            d1 = self.thickness.get("tails")
            d2 = self.thickness.get("head")

        # set new membrane thickness if intended
        elif "Monolayer" in self.lipidNames:
            if config.updateMonolayerThickness == True:
                d1 = config.new_d1
                d2 = config.new_d2
                print("\nYou have selected to update the monolayer thickness such that d1 = %.3f and d2 = %.3f." %(d1,d2))

            else:
                d1 = self.thickness.get("tails")
                d2 = self.thickness.get("head")
                print("\nMonolayer thicknesses are as in instructions file; d1 = %.3f and d2 = %.3f." %(d1,d2))

        # calculate volume fraction
        self.headVolFrac = (rho1 * d1 * SL2 ) / ( rho2 * d2 * SL1)
        self.twoSolv     = (1-self.headVolFrac)*100

        # solvent volume in head group
        if config.verbose == True:
            print("\nHead volume fraction: %f" %self.headVolFrac)
            print("\n2-solv = %f" %self.twoSolv)

        return (1-self.headVolFrac)



    # method that returns monolayer information for adding extra lipid components
    def getMonolayerPar(self):
        return (self.monolayerMolVol, self.monolayerSL, self.monolayerSLD)


# function to check whether a string contains a number
def hasNumbers(inputString):
    return bool(re.search(r'\d', inputString))



def mainCalcPar(calcType,membrane,lipidRatio,t_thick,h_thick):

    # get list of lipids within membrane with ratio
    lipids = re.split(':',membrane)
    ratios = re.split(':',str(lipidRatio))

    # structure thickness information into dict
    thickness = {
        'head':  float(h_thick),
        'tails': float(t_thick)
    }

    # initialise monolayer information
    monolayerMolVol = {'head': 0, 'tails': 0}
    monolayerSL     = {'head': 0, 'tails': 0}
    monolayerSLD    = {'head': 0, 'tails': 0}
    monolayerPar    = (monolayerMolVol, monolayerSL, monolayerSLD)


    # create class instance with input variables
    m = Membrane(lipids, ratios, thickness, monolayerPar)

    # converts 3:5 to 3/8:5/8
    m.normaliseMolarRatios()

    # calculates the total lipid volume of the monolayer (needed for addLipidToMonolayer)
    m.calcTotalLipidVol()

    # convert molar fraction to component volume fraction
    if config.useVolFrac == True: m.calcVolFrac()

    # calculate coherent scattering lengths
    m.calcSL()

    # calculate average lipid structure volumes
    m.calcAvLipidVol()

    # divide scattering length by the molecular volume
    m.calcSLD()

    # calculate volume fraction of the headgroups
    solvFrac = m.calcHeadVolumeFraction()
    
    # get monolayer parameters: (self.monolayerMolVol, self.monolayerSL, self.monolayerSLD)
    monolayerPar = m.getMonolayerPar()
    
    if calcType.upper() == 'SL':
        return monolayerPar[1]

    if calcType.upper() == 'SLD':
        return monolayerPar[2]
    
    if calcType.upper() == 'solvFrac':
        return solvFrac


