import os.path
import numpy as np
import matplotlib.pyplot as plt
import scipy

import refnx
from refnx.dataset import ReflectDataset, Data1D
from refnx.analysis import Transform, CurveFitter, Objective, GlobalObjective, Model, Parameter, process_chain
from refnx.reflect import SLD, Slab, ReflectModel, Structure, LipidLeaflet

from sldAnalysis import mainCalcSolvFrac

print(f"refnx: {refnx.version.version}\n"
  f"scipy: {scipy.version.version}\n"
  f"numpy: {np.version.version}\n")

file_paths = ['../input/S14_excel.txt','../input/S11_excel.txt']

expProperties = [{'solvType': 'd2o','monolayerType': 'h-MC3'},\
                 {'solvType': 'acmw','monolayerType': 'd-MC3'}]

# create shared sld slabs 
acmw = SLD(0.0, name='acmw')
d2o  = SLD(6.36, name='d2o')
air  = SLD(0.0, name='air')

# create parameters to be shared between leaflets 
# APM = (SL1/SLD1)d1
apm       = Parameter(55.0, name='APM', bounds=(0,1000), vary=True)

b_heads   = Parameter(18.881, name='SL2')
vm_heads  = Parameter(260.0, name='VM2')

b_tails_h = Parameter(-7.5180, name='SL1_h')
b_tails_d = Parameter(637.9020, name='SL1_d')
vm_tails  = Parameter(1030.0, name='VM1')

thickness_heads = Parameter(10.0, bounds=(5,15), vary=True, name='d2')
thickness_tails = Parameter(20.0, bounds=(15,25), vary=True, name='d1')

roughness = Parameter(3.5, name='Roughness')

tail_solvent = Parameter(0, name='Solv1')

# define lipid membrane leaflets 
# fronting medium is air therefore reverse_monolayer = True
leaflet_d2o_hlip = LipidLeaflet(apm, b_heads, vm_heads, thickness_heads, \
                             b_tails_h, vm_tails, thickness_tails, \
                             roughness, roughness, tail_solvent, \
                             reverse_monolayer=True, name='d2o_hlip')


leaflet_acwm_dlip = LipidLeaflet(apm, b_heads, vm_heads, thickness_heads, \
                                 b_tails_d, vm_tails, thickness_tails, \
                                 roughness, roughness, tail_solvent, \
                                 reverse_monolayer=True, name='acwm_dlip')


#### leaflet_acwm_hlip = LipidLeaflet(apm, b_heads, vm_heads, thickness_heads, \
####    b_tails_h, vm_tails, thickness_tails, rough_everything, rough_everything, acmw, air, reverse_monolayer=True, name='acwm_dlip')

#### leaflet_d2o_d = LipidLeaflet(apm, b_heads, vm_heads, thickness_heads, \
####    b_tails_d, vm_tails, thickness_tails, rough_everything, rough_everything, d2o, air, reverse_monolayer=True, name='acwm_dlip')


# set constraints on apm  (SL1/SLD1)/d1
#apm.constraint = vm_heads/thickness_tails
#print(apm)

# define structure objects 
s_d2O_hlip  = air | leaflet_d2o_hlip | d2o(0, roughness)
s_acmw_dlip = air | leaflet_acwm_dlip | acmw(0, roughness)

# s_acmw_hlip = air | leaflet_acwm_hlip | acmw(0, roughness)
# s_d2O_dlip = air | leaflet_d2o_d | d2o(0, roughness)

# generate model data and set background 
model_d2O_hlip  = ReflectModel(s_d2O_hlip)
model_d2O_hlip.bkg.setp(vary=True, bounds=(4.5e-7, 5.5e-7))

model_acmw_dlip = ReflectModel(s_acmw_dlip)
model_acmw_dlip.bkg.setp(vary=True, bounds=(4.5e-7, 5.5e-7))


# convert to reflect dataset object
data_d2O_hlip  = ReflectDataset(file_paths[0])
data_acmw_dlip = ReflectDataset(file_paths[1])

# generate objective for given model (model, experiment)
obj_d2O_hlip  = Objective(model_d2O_hlip, data_d2O_hlip)
obj_acmw_dlip = Objective(model_acmw_dlip, data_acmw_dlip)

# append objective to list
obj_list = [obj_d2O_hlip,obj_acmw_dlip]

# combine all individual objectives into a GlobalObjective
global_objective = GlobalObjective(obj_list)

# curve fitter object performs leaste squares fitting or MCMC sampling
fitter = CurveFitter(global_objective, nwalkers=1e4)

# do fit
fitter.fit('differential_evolution')#, constraints=(constraint,));
print(global_objective)





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
# plotSLD = True
# if plotSLD == True:
#     fig, ax = plt.subplots()
#     ax.plot(*s_List[i].sld_profile(), label='')
#     ax.set_ylabel("$\\rho$ / $10^{-6} \AA^{-2}$")
#     ax.set_xlabel("z / $\AA$")
#     ax.legend()




## monolayer at air-liquid interface with few contrasts
## lipid leaflet is one half of bilayer
## can give refnx a .analysis.parameter
## e.g. from refnx.analysis import Parameter; my_particular_parameter(55.0, name="APM")
## leaflet1 = LipidLeaflet(my_particular_parameter,..)
## leaflet2 = LipidLeaflet(my_particular_parameter,..) two leaflets that use the same APM


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



