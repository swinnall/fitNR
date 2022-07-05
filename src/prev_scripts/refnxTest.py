import os.path
import numpy as np
import matplotlib.pyplot as plt
import scipy

import refnx
from refnx.dataset import ReflectDataset, Data1D
from refnx.analysis import Transform, CurveFitter, Objective, GlobalObjective, Model, Parameter, process_chain
from refnx.reflect import SLD, Slab, ReflectModel, Structure, LipidLeaflet

from sldAnalysis import mainCalcSolvFrac


def getFilePaths():

    file_paths = ['../input/S14_excel.txt','../input/S11_excel.txt']

    expProperties = [{'solvType': 'd2o','monolayerType': 'h-MC3'},\
                     {'solvType': 'acmw','monolayerType': 'd-MC3'}]

    return file_paths, expProperties



def buildModel(expProperties):

    # # parameter for the area per molecule each molecule occupies at the surface
    # apm = Parameter(60, 'area per molecule', vary=True, bounds=(50, 100))


    # # the sum of scattering lengths for the lipid head and tail in Angstrom
    # if expProperties.get("monolayerType") == 'h-MC3':
    #     b_heads = Parameter(18.8810e-4, 'b_heads')
    #     b_tails = Parameter(-7.5180e-4, 'b_tails')

    # elif expProperties.get("monolayerType") == 'd-MC3':
    #     b_heads = Parameter(18.8810e-4, 'b_heads')
    #     b_tails = Parameter(637.9020e-4, 'b_tails')

    # # the volume occupied by the head and tail groups in cubic Angstrom
    # v_heads = Parameter(260, 'v_heads')
    # v_tails = Parameter(1030, 'v_tails')

    # # check SLD
    # #print(b_heads/v_heads)
    # #print(b_tails/v_tails)

    # # the head and tail group thicknesses.
    # head_thickness = Parameter(6.5, 'head_thickness', vary=False, bounds=(5, 10))
    # tail_thickness = Parameter(10.9, 'tail_thickness', vary=False, bounds=(10, 17))

    # # define the head solvent fraction
    # head_solvent = Parameter(mainCalcSolvFrac(0,tail_thickness,head_thickness), 'head_solvent', vary=False, bounds=(0,1))
    #   # in version 0.1.30 you can set constraints
    #   # head_solvent.set_constraint(mainCalcSolvFrac(0,tail_thickness,head_thickness))
    # tail_solvent = Parameter(0, 'tail_solvent')

    # # roughness is calculated for all layers
    # roughness = 3.5

    # # tail groups face upwards
    # leaflet = LipidLeaflet(apm,
    #                              b_heads, v_heads, head_thickness,
    #                              b_tails, v_tails, tail_thickness,
    #                              roughness, roughness,
    #                              #head_solvent=(head_solvent), tail_solvent=(tail_solvent),
    #                              # tail_solvent=(tail_solvent),
    #                              reverse_monolayer=False)


    ###### TEST ######
    # roughness is calculated for all layers
    roughness = 3.5

    # solvent contrasts; vary `real` parameter with uniform bounds; `setp` parameter method changes many aspects at once
    if expProperties.get("solvType") == 'd2o':
        solvent = SLD(6.36 + 0j)
        solvent.real.setp(vary=True, bounds=(6.1, 6.4))
        solvent.real.name='d2o SLD'

    elif expProperties.get("solvType") == 'acmw':
        solvent = SLD(0 + 0j)
        solvent.real.setp(vary=True, bounds=(0, 0.02))
        solvent.real.name='acmw SLD'
    solvent.name = 'buffer'

    # define tail as slab object
    if expProperties.get("monolayerType") == 'h-MC3':
        tail = SLD(-0.0730 + 0j)
        tail.real.setp(vary=False)

    elif expProperties.get("monolayerType") == 'd-MC3':
        tail = SLD(6.1932 + 0j)
        tail.real.setp(vary=False)

    tail.name = 'tails'
    tailSlab = tail(10.9,roughness)
    tailSlab.thick.setp(10.9,vary=False,bounds=(10,17))
    tailSlab.thick.name = 'tails thickness'
    tailSlab.rough.setp(roughness,vary=False)
    tailSlab.rough.name = 'tails roughness'
    tailSlab.vfsolv.setp(0,vary=False)
    tailSlab.vfsolv.name = 'tails solvation fraction'


    head = SLD(0.7262 + 0j)
    head.name = 'heads'
    headSlab = head(6.5,roughness)
    headSlab.thick.setp(6.5,vary=False,bounds=(5,10))
    headSlab.thick.name = 'heads thickness'
    headSlab.rough.setp(roughness,vary=False)
    headSlab.rough.name = 'heads roughness'
    headSlab.vfsolv.setp(mainCalcSolvFrac(0,tailSlab.thick,headSlab.thick),vary=False)
    headSlab.vfsolv.name = 'heads solvation fraction'

    # make air slab
    air = SLD(0 + 0j)
    air.name = 'air'

    # make structure
    #s_ = air(0,roughness) | leaflet | solvent(0,roughness) # (10.9,roughness)
    s_ = air(0,roughness) | tailSlab | headSlab | solvent(0,roughness)

    ###
    'Problem: only by creating models with common slabs within structures \
        can you constrain parameters to one another. E.g. only solvent differing \
        This cant be done as SL always differs so always need a new slab/leaflet\
        Would need to constrain underpinning parameters... E.g. thickness \
        Could couple hMC3 in D2O and ACMW but not to d-MC3 in the same contrasts'

    # generate model data
    model_ = ReflectModel(s_)

    model_.scale.setp(vary=False, bounds=(0.9, 1.1))
    model_.bkg.setp(vary=True, bounds=(4.5e-7, 5.5e-7))

    return s_, model_


## monolayer at air-liquid interface with few contrasts
## lipid leaflet is one half of bilayer
## can give refnx a .analysis.parameter
## e.g. from refnx.analysis import Parameter; my_particular_parameter(55.0, name="APM")
## leaflet1 = LipidLeaflet(my_particular_parameter,..)
## leaflet2 = LipidLeaflet(my_particular_parameter,..) two leaflets that use the same APM
## acmw = SLD(0.0)
## d2O = SLD(6.36)
## air = SLD(0.0)
## apm = Parameter(55.0, name='APM')
## b_heads = Parameter(200.0)
## vm_heads = Parameter(300.0)
## b_tails_d = Parameter(300)
## b_tails_h = Parameter(-100.0)
## vm_tails = Parameter(700.0)

# thickness_tails = Parameter(25.0)
# thickness_heads = Parameter(15.0)
# rough_everything = Parameter(3.5)
# head_solvent

# fronting medium is air therefore reverse_monolayer = True

#### leaflet_acwm_dlip = LipidLeaflet(apm, b_heads, vm_heads, thickness_heads, \
####    b_tails_d, vm_tails, thickness_tails, rough_everything, rough_everything, acmw, air, reverse_monolayer=True, name='acwm_dlip')

#### leaflet_d2o_h = LipidLeaflet(apm, b_heads, vm_heads, thickness_heads, \
####    b_tails_h, vm_tails, thickness_tails, rough_everything, rough_everything, d2o, air, reverse_monolayer=True, name='acwm_dlip')

#### leaflet_acwm_hlip = LipidLeaflet(apm, b_heads, vm_heads, thickness_heads, \
####    b_tails_h, vm_tails, thickness_tails, rough_everything, rough_everything, acmw, air, reverse_monolayer=True, name='acwm_dlip')

#### leaflet_d2o_d = LipidLeaflet(apm, b_heads, vm_heads, thickness_heads, \
####    b_tails_d, vm_tails, thickness_tails, rough_everything, rough_everything, d2o, air, reverse_monolayer=True, name='acwm_dlip')


# s_acmw_dlip = air | leaflet_acwm_dlip | acmw(0, rough_everything)
# s_d2O_hlip = air | leaflet_d2o_h | d2o(0, rough_everything)
# s_acmw_hlip = air | leaflet_acwm_hlip | acmw(0, rough_everything)
# s_d2O_dlip = air | leaflet_d2o_d | d2o(0, rough_everything)

## sio2 = SLD(3.47)
#### sio2_slab = slab(10,sio2,3)     can either give it a 10 or do sio2_thickness

### leaflet_d2o_d.apm.constraint = 2 * leaflet_acwm_dlip example constraint
### def my_constraint(x): ...
### then leaflet_d2O_hlip.amp.constraint

### leaflet_acmw_dlip.logp()

## least squares fit impleemnt inequality fit. Get the four contrast
## check least squares if reasonable, sampling estimates the uncertainty in the system and should be final step
## sampling - covariance matrix is calculated from least square sfit which gives uncertainty in the parameters and generate series of solutions with same statistics and use that to initalise the walks in the markov chain
## recommended jitter option which improves bad walker position at the beginning if volFrac is in a bad situation
##

## leaflet works for single component, preferable for users.
# when you go to multicomponent systems, average values of scattering lengths etc can you convince yourself apm is valid

# if not go back to slabs

# would have to calculate head volume if have insertion into

## leaflet for average components possible
## fit with slabs for. Drug structure is in between d2o and leaflet

## vase data, single angle single wavelength, can still fit !!
## use same model for refnx to fit ellipsometry
## refnx tools/app/lipids_sld_calculator.xlsx repository of lipid properties. DPPC/DMPC etc. Lipid browser inside refnx gui


def plotFit(global_objective, s_List):

    plotObjective = True
    if plotObjective == True:
        global_objective.plot()
        plt.yscale('log')
        plt.yscale('log')
        plt.xlabel('Q / $\AA^{-1}$')
        plt.ylabel('Reflectivity')
        plt.legend()

    # plt.savefig('global.png',
    #     format='png',
    #     dpi=400,
    #     bbox_inches='tight')

    # plot SLD
    plotSLD = True
    if plotSLD == True:
        fig, ax = plt.subplots()
        for i in range(len(s_List)):
            ax.plot(*s_List[i].sld_profile(), label='')
        ax.set_ylabel("$\\rho$ / $10^{-6} \AA^{-2}$")
        ax.set_xlabel("z / $\AA$")
        ax.legend()

    return



def main():

    # initialise variables
    objList = []; s_List = []

    # get path names of the NR experimental datasets
    file_paths, expProperties = getFilePaths()

    for idx, path in enumerate(file_paths):

        # convert to reflect dataset object
        data_ = ReflectDataset(path)

        # generate model
        s_, model_ = buildModel(expProperties[idx])

        # generate objective for given model (model, experiment)
        objective_ = Objective(model_, data_)

        # append objective to list
        objList.append(objective_)
        s_List.append(s_)

    # combine all individual objectives into a GlobalObjective
    global_objective = GlobalObjective(objList)

    # curve fitter object performs leaste squares fitting or MCMC sampling
    fitter = CurveFitter(global_objective, nwalkers=1e4)

    # do fit
    fitter.fit('differential_evolution')#, constraints=(constraint,));
    print(global_objective)

    # plot the output fits of all the samples
    plotFit(global_objective, s_List)

    return



if __name__ == '__main__':
    print(f"refnx: {refnx.version.version}\n"
      f"scipy: {scipy.version.version}\n"
      f"numpy: {np.version.version}\n")
    main()
