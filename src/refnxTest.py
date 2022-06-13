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

# get test sample 
file_path = '../input/S14_excel.txt'
data_d2o = ReflectDataset(file_path)

# solvent contrasts
d2o   = SLD(6.36 + 0j)
acmw  = SLD(0.00 + 0j)

# vary `real` parameter with uniform bounds
# `setp` parameter method changes many aspects at once
d2o.real.setp(vary=True, bounds=(6.1, 6.4))
d2o.real.name='d2o SLD'

# Parameter for the area per molecule each DMPC molecule occupies at the surface. We
# use the same area per molecule for the inner and outer leaflets.
apm = Parameter(100, 'area per molecule', vary=False, bounds=(90, 110))

# the sum of scattering lengths for the lipid head and tail in Angstrom
b_heads = Parameter(18.8810e-4, 'b_heads')
b_tails = Parameter(-7.5180e-4, 'b_tails')

# the volume occupied by the head and tail groups in cubic Angstrom
v_heads = Parameter(260, 'v_heads')
v_tails = Parameter(1030, 'v_tails')

# the head and tail group thicknesses.
head_thickness = Parameter(6.5, 'head_thickness', vary=False, bounds=(5, 10))
tail_thickness = Parameter(10.9, 'tail_thickness', vary=False, bounds=(10, 17))

# define the head solvent fraction 
head_solvent = Parameter(mainCalcSolvFrac(0,tail_thickness,head_thickness), 'head_solvent', vary=False, bounds=(0,1))
# in version 0.1.30 you can set constraints 
#head_solvent.set_constraint(mainCalcSolvFrac(0,tail_thickness,head_thickness))
tail_solvent = Parameter(0, 'tail_solvent')

# check SLD
#print(b_heads/v_heads)
#print(b_tails/v_tails)

# roughness is calculated for all layers 
roughness = 3.5

# tail groups face upwards
leaflet = LipidLeaflet(apm,
                             b_heads, v_heads, head_thickness,
                             b_tails, v_tails, tail_thickness,
                             roughness, roughness, 
                             head_solvent=(head_solvent), tail_solvent=(tail_solvent),
 #                            tail_solvent=(tail_solvent),
                             reverse_monolayer=False)

# make air slab 
air = SLD(0 + 0j)

# make structure 
s_d2o = air(0,roughness) | leaflet | d2o(0,roughness)

# generate model data 
model_d2o = ReflectModel(s_d2o)

model_d2o.scale.setp(vary=False, bounds=(0.9, 1.1))
model_d2o.bkg.setp(vary=True, bounds=(4.5e-7, 5.5e-7))

objective_d2o = Objective(model_d2o, data_d2o)

# combine all individual objectives into a GlobalObjective 
global_objective = GlobalObjective([objective_d2o])

# curve fitter object performs leaste squares fitting or MCMC sampling 
fitter = CurveFitter(global_objective, nwalkers=1e3)



# from scipy.optimize import NonlinearConstraint

# class DEC(object):
#     def __init__(self, pars, objective):
#         # we'll store the parameters and objective in this object
#         # this will be necessary for pickling in the future
#         self.pars = pars
#         self.objective = objective

#     def __call__(self, x):
#         # we need to update the varying parameters in the objective first
#         self.objective.setp(x)
#         #print(self.pars[0],self.pars[1])
#         # print(mainCalcSolvFrac(0,self.pars[0],self.pars[1]))
        
#         #print(pars)
#         headSolvFrac = mainCalcSolvFrac(0,self.pars[0],self.pars[1])
        
#         if self.pars[2] == headSolvFrac:
#             return self.pars[2]
#         else: 
#             return 2
        
        
    
#     ###return mainCalcSolvFrac(0,self.pars[0],self.pars[1])

# pars = (s_d2o[1].thickness_tails, s_d2o[1].thickness_heads, s_d2o[1].head_solvent)
# dec = DEC(pars, global_objective)

# constraint = NonlinearConstraint(dec, 0, 1)

# do fit and plot output 
fitter.fit('differential_evolution')#, constraints=(constraint,));
print(global_objective)




# plot fit 
global_objective.plot()
plt.yscale('log')
plt.yscale('log')
plt.xlabel('Q / $\AA^{-1}$')
plt.ylabel('Reflectivity')
plt.legend()

# plot SLD 
fig, ax = plt.subplots()
ax.plot(*s_d2o.sld_profile(), label='d2o')
ax.set_ylabel("$\\rho$ / $10^{-6} \AA^{-2}$")
ax.set_xlabel("z / $\AA$")
ax.legend();


