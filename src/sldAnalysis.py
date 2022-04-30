" Calculate SLD profiles "

import csv
import glob, os
import re
import sys
import config

class Membrane:

    def __init__(self, lipidNames, molRatios, thickness, monolayerPar, outputFilePath):

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
        self.outputFilePath = outputFilePath

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



    # calculates 2-SLD for Igor Motofit; used when merging L3 drug with L2 headgroup
    def calcSLD_drugToMonolayer(self):

        # check only one drug component being added to the monolayer
        if len(config.injectedDrugNames) > 1:
            warningChoice = input("\n\n\nConfig Warning: Attempted insertion of > 1 drug component into headgroup.\
                \nThis feature is not yet available!\
                \nPluto will only take first component.\
                \nDo you want to continue? (y/n)\n")

            if warningChoice.upper() == 'N':
                print("Session closed.")
                sys.exit()

        # update thickness as part of the drug enters the second layer
        self.d3 = config.injectedDrugSizes[0] - self.thickness.get("head")

        # get solvent amount in third layer from neutron fit, set in config
        self.threeSolv = config.threeSolv

        # calculate the volume fraction of drug in the third layer
        vf_drug = (100 - self.threeSolv) / 100

        # update 2-Solv where it is assumed the vol.frac. of drug is the same in both layers
        self.twoSolv = self.twoSolv - (vf_drug*100)

        # get SLD values of drug in H2O and D2O
        drugSLD_H2O  = config.injectedDrugSLD_H2O.get(config.injectedDrugNames[0])
        drugSLD_D2O  = config.injectedDrugSLD_D2O.get(config.injectedDrugNames[0])

        # averages the existing SLDs of layer 2 (headgroup layer) with the added drug
        self.twoSLD_H2O = (self.headVolFrac*self.avSLD.get("head") + vf_drug*drugSLD_H2O ) / (self.headVolFrac + vf_drug)
        self.twoSLD_D2O = (self.headVolFrac*self.avSLD.get("head") + vf_drug*drugSLD_D2O ) / (self.headVolFrac + vf_drug)

        # print values to terminal
        print("\n\nYou have chosen to add the drug (%s) to the SECOND layer:" %config.injectedDrugNames[0])
        print("\n3-thick = %f" %self.d3)
        print("\n2-solv = %f\n3-solv = %f" %(self.twoSolv, self.threeSolv))
        print("\n2-SLD_H2O = %f\n2-SLD_D2O = %f" %(self.twoSLD_H2O, self.twoSLD_D2O))


    def calcSLD_drugToThirdLayer(self):

        # checks there are more than 1 injected drug components
        if len(config.injectedDrugNames) < 2:
            print("\nFatal Config Error: Insufficient drug components to calculate.")
            sys.exit()

        # get drug name and size
        drugName = config.injectedDrugNames

        tot = 0
        for i, value in enumerate(config.injectedDrugRatios):
            tot += float(value)

        normDrugMolRatios = []
        for i, value in enumerate(config.injectedDrugRatios):
            normDrugMolRatios.append( float(value) / tot  )

        # calculates SLD for the overlapping region of the drug components
        for contrast in ["H2O", "D2O"]:
            for i, drugName in enumerate(config.injectedDrugNames):

                if contrast == "H2O":
                    self.drugOverlapSLD_H2O += normDrugMolRatios[i] * config.injectedDrugSLD_H2O.get(drugName)

                if contrast == "D2O":
                    self.drugOverlapSLD_D2O += normDrugMolRatios[i] * config.injectedDrugSLD_D2O.get(drugName)

        print("\nSLD of overlapping drug in H2O: %.4f" %self.drugOverlapSLD_H2O)
        print("SLD of overlapping drug in D2O: %.4f" %self.drugOverlapSLD_D2O)




    # method that returns monolayer information for adding extra lipid components
    def getMonolayerPar(self):
        return (self.monolayerMolVol, self.monolayerSL, self.monolayerSLD)


    # write output to file
    def appendFile_Lipid(self):

        with open(self.outputFilePath, 'a', newline = '') as f:

            if "Monolayer" not in self.lipidNames:
                f.write('\n~~~~\nCalculated Membrane: %s\n' %self.lipidNames)
            else:
                f.write('\nCalculated Membrane: %s\n' %self.lipidNames)

            f.write('Corresponding ratio: %s\n' %self.molRatios)
            f.write("Average SL;  Head = %.4f; Tail = %.4f\n" %(self.avSL.get("head"),self.avSL.get("tails")))
            f.write("Average SLD; Head = %.4f; Tail = %.4f\n" %(self.avSLD.get("head"),self.avSLD.get("tails")))
            f.write("Thickness;   Head = %.4f; Tail = %.4f\n" %(self.thickness.get("head"),self.thickness.get("tails")))
            f.write("Head vol frac = %.4f\n" %self.headVolFrac)
            f.write("2-solv = %.4f\n" %self.twoSolv)


    def appendFile_drugBinding(self):
        with open(self.outputFilePath, 'a', newline = '') as f:

            if config.addDrugToMonolayer == True:
                f.write("\nDrug (%s) added to layer 2 of monolayer:\n" %config.injectedDrugNames[0])
                f.write("d3 = %f\n" %self.d3)
                f.write("2-solv = %f\n" %self.twoSolv)
                f.write("3-solv = %.4f\n" %self.threeSolv)
                f.write("2-SLD_H2O = %f\n" %self.twoSLD_H2O)
                f.write("2-SLD_D2O = %f\n" %self.twoSLD_D2O)

            if config.addDrugToThirdLayer == True:
                f.write("\nDrug %s added to layer 3 of monolayer:\n" %config.injectedDrugNames)
                f.write('Corresponding ratio: %s\n' %config.injectedDrugRatios)
                f.write("SLD of overlapping drug in H2O: %.4f\n" %self.drugOverlapSLD_H2O)
                f.write("SLD of overlapping drug in D2O: %.4f" %self.drugOverlapSLD_D2O)




def importSampleData(instructionsFile, sampleNum):

    membrane   = instructionsFile["membranes"][sampleNum]
    lipidRatio = instructionsFile["lipidRatios"][sampleNum]
    t_thick    = instructionsFile["d1"][sampleNum]
    h_thick    = instructionsFile["d2"][sampleNum]
    label      = instructionsFile["label"][sampleNum]

    return membrane, lipidRatio, t_thick, h_thick, label


# function to check whether a string contains a number
def hasNumbers(inputString):
    return bool(re.search(r'\d', inputString))



def main(instructionsFile, outputFilePath):

    # number of membranes to calculate
    nMemb = len(instructionsFile)

    # calculate component volumes
    for sampleNum in range(nMemb):

        # import sample data
        membrane, lipidRatio, t_thick, h_thick, label = importSampleData(instructionsFile, sampleNum)

        # print membrane label to terminal
        print("\n\n\n~~~\nMembrane %d - %s" %(sampleNum+1,label))

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
        m = Membrane(lipids, ratios, thickness, monolayerPar, outputFilePath)

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
        m.calcHeadVolumeFraction()

        # append file with initial monolayer calculation information
        m.appendFile_Lipid()


        # repeat for incorporating an injected lipid component into the membrane
        if config.addLipidToMonolayer == True:

            # get monlayer parameter information for next iteration
            monolayerPar = m.getMonolayerPar()

            if config.very_verbose == True: print("\nMonolayerPar:\n%s" %(monolayerPar,))

            # get added lipids and associated ratios from config file
            lipids = config.injectedLipidNames
            ratios = config.injectedLipidRatios

            print("\n\n\nYou have chosen to add components from injected sample to the monolayer at the following ratios:")
            sumRatiosTest = 0
            for ele, lipid in enumerate(lipids):
                sumRatiosTest += ratios[ele]
                if lipid == "Monolayer": print("Averaged %s: %d%%." %(lipid, ratios[ele]))
                else: print("Added Lipid: %s: %d%%." %(lipid, ratios[ele]))

            if sumRatiosTest != 100:
                print("Fatal Config Error: Injected lipid molar ratios defined in config do not equal 100.")
                sys.exit()


            # creates a new class instance to pass new config params to
            m = Membrane(lipids, ratios, thickness, monolayerPar, outputFilePath)

            m.normaliseMolarRatios()

            m.calcTotalLipidVol()

            if config.useVolFrac == True: m.calcVolFrac()

            m.calcSL()

            m.calcAvLipidVol()

            m.calcSLD()

            m.calcHeadVolumeFraction()

            m.appendFile_Lipid()


        #if config.addDrugToThirdLayer == True and config.addDrugToMonolayer == True:
        #    print("\nFatal Config Error: Both addDrug config commands are true.")
        #    sys.exit()

        # calculates SLD for mixing drug in layer 3 to layer 2
        if config.addDrugToMonolayer == True:
            m.calcSLD_drugToMonolayer()

        # calculatese SLD for adding drug to layer 3 when multiple drug components
        if  config.addDrugToThirdLayer == True:

            print("\n\n\nYou have chosen to add drug components to the THIRD layer at the following ratios:")
            sumRatiosTest = 0
            for ele, drug in enumerate(config.injectedDrugNames):
                sumRatiosTest += config.injectedDrugRatios[ele]
                print("Added drug: %s: %d%%." %(drug, config.injectedDrugRatios[ele]))

            if sumRatiosTest != 100:
                print("\nFatal Config Error: Injected drug molar ratios defined in config do not equal 100.")
                sys.exit()

            m.calcSLD_drugToThirdLayer()

        # write drug binding information to file, this indentation prevents writing twice
        m.appendFile_drugBinding()


    sys.exit()
    return



if __name__ == '__main__':
    print("~Running sldAnalysis.py~")
    main()
