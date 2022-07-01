import os.path
import numpy as np
import matplotlib.pyplot as plt
import scipy

import corner 
import refnx
from refnx.dataset import ReflectDataset, Data1D
from refnx.analysis import Transform, CurveFitter, Objective, GlobalObjective, Model, Parameter, process_chain
from refnx.reflect import SLD, Slab, ReflectModel, Structure, LipidLeaflet
from scipy.optimize import NonlinearConstraint

from sldAnalysis import mainCalcSolvFrac

print(f"refnx: {refnx.version.version}\n"
  f"scipy: {scipy.version.version}\n"
  f"numpy: {np.version.version}\n")

file_paths = ['../input/hMC3-D2O.txt','../input/dMC3-ACMW.txt',\
              '../input/dMC3-D2O.txt']

expProperties = [{'solvType': 'd2o','monolayerType': 'h-MC3'},\
                 {'solvType': 'acmw','monolayerType': 'd-MC3'},\
                 {'solvType': 'd2o','monolayerType': 'd-MC3'}]

# create sld slabs 
acmw = SLD(0.0, name='acmw')
d2o  = SLD(6.36,name='d2o')
air  = SLD(0.0, name='air')

# SLD values set as parameters for headSolv constraint 
rho_tails_h = Parameter(-0.0730, name='SLD1_h')
rho_tails_d = Parameter(6.19320, name='SLD1_d')
rho_heads   = Parameter(0.72620, name='SLD2')

# SLD values set as parameters for headSolv constraint 
b_tails_h = Parameter(-7.5180, name='SL1_h')
b_tails_d = Parameter(637.9020, name='SL1_d')
b_heads   = Parameter(18.881, name='SL2')

tails_h = SLD(rho_tails_h, name='tails_h')
tails_d = SLD(rho_tails_d, name='tails_d')
heads   = SLD(rho_heads, name='heads')

# create parameters to be shared between slabs 
roughness       = Parameter(3.5, name='Roughness')

thickness_tails = Parameter(11.388, bounds=(10,35), vary=True, name='Tails Thickness')
tail_solvent    = Parameter(0, name='Solv1')

thickness_heads = Parameter(6.0, bounds=(5,15), vary=True, name='Heads Thickness')
head_solvent    = Parameter(0.0, name='Solv2')


# constrain head_solvent = 1-headVolFrac = 1-(rho1*d1*SL2)/(rho2*d2*SL1)
head_solvent.constraint = 1-((rho_tails_d*thickness_tails*b_heads)/(rho_heads*thickness_heads*b_tails_d))
#print(head_solvent)

#test_headSolv = mainCalcSolvFrac(modelNum=0,t_thick=thickness_tails,h_thick=thickness_heads)
#print(test_headSolv)
#print('\n\n')

# create the tails slab
tails_layer_h = tails_h(thickness_tails,roughness,tail_solvent) 
tails_layer_d = tails_d(thickness_tails,roughness,tail_solvent)  
heads_layer   = heads(thickness_heads,roughness,head_solvent)  


# define structure objects 
s_d2O_hlip  = air | tails_layer_h | heads_layer | d2o(0, roughness)
s_acmw_dlip = air | tails_layer_d | heads_layer | acmw(0, roughness)
s_d2o_dlip  = air | tails_layer_d | heads_layer | d2o(0, roughness)


# generate model data and set background 
model_d2O_hlip  = ReflectModel(s_d2O_hlip)
model_d2O_hlip.bkg.setp(5e-7,vary=False, bounds=(4.5e-7, 5.5e-7))

model_acmw_dlip = ReflectModel(s_acmw_dlip)
model_acmw_dlip.bkg.setp(5e-7,vary=False, bounds=(4.5e-7, 5.5e-7))

model_d2o_dlip = ReflectModel(s_d2o_dlip)
model_d2o_dlip.bkg.setp(5e-7,vary=False, bounds=(4.5e-7, 5.5e-7))

# convert to reflect dataset object
data_d2O_hlip  = ReflectDataset(file_paths[0])
data_acmw_dlip = ReflectDataset(file_paths[1])
data_d2o_dlip  = ReflectDataset(file_paths[2])

# generate objective for given model (model, experiment)
obj_d2O_hlip  = Objective(model_d2O_hlip, data_d2O_hlip)
obj_acmw_dlip = Objective(model_acmw_dlip, data_acmw_dlip)
obj_d2o_dlip  = Objective(model_d2o_dlip, data_d2o_dlip)


# append objective to list
obj_list = [obj_d2O_hlip,obj_acmw_dlip,obj_d2o_dlip]

# combine all individual objectives into a GlobalObjective
global_objective = GlobalObjective(obj_list)

# curve fitter object performs leaste squares fitting
fitter = CurveFitter(global_objective)

#constrain t_thick > h_thick 
class DEC(object):
    def __init__(self, pars, objective):
        # we'll store the parameters and objective in this object
        # this will be necessary for pickling in the future
        self.pars = pars
        self.objective = objective

    def __call__(self, x):
        # we need to update the varying parameters in the
        # objective first
        self.objective.setp(x)
        return float(self.pars[0] - self.pars[1])
    
pars = (thickness_tails, thickness_heads)
dec = DEC(pars, global_objective)

thickness_constraint = NonlinearConstraint(dec, 0, np.inf)



# do fit
fitter.fit('differential_evolution', constraints=(thickness_constraint));
print(global_objective)
print('Reduced chiSq = %f' %(global_objective.chisqr()/(3*len(data_d2O_hlip)-4)))

# plot fits 
plotObjective = True
if plotObjective == True:
    fig, ax = global_objective.plot()
    
    plt.yscale('log')
    plt.xlabel('Q / $\AA^{-1}$', fontsize=18, fontweight='bold')
    plt.ylabel('R', fontsize=18, fontweight='bold')
    plt.legend(prop={'size': 16, 'weight':'bold'}, frameon = False, loc='best')
      
    plt.tick_params(axis='both',which='major', labelsize=13, size=8, width=2, direction='in', top='on')
    plt.tick_params(axis='both',which='minor', labelsize=13, size=4, width=2, direction='in', right='on')
   
    plt.savefig('../output/globalObj.png',
        format='png',
        dpi=400,
        bbox_inches='tight') 



class LogpExtra(object):
    def __init__(self, pars):
        # we'll store the parameters and objective in this object
        # this will be necessary for pickling in the future
        self.pars = pars

    def __call__(self, model, data):
        if float(self.pars[0] - self.pars[1]) > 0:
            print('d1-d2>0')
            return 0
        else: print('not satisfied')
        return -np.inf
    
lpe = LogpExtra(pars)

# set the log_extra attribute of the Objective with our extra log-probability term.
global_objective.logp_extra = lpe

    


doMCMC = True  
if doMCMC == True:
    
    # MCMC - burn first 400 as initial chain might not be representative of equilibrated system  
    fitter.sample(steps=400, random_state=1, pool=-1)
    fitter.sampler.reset()
    
    # MCMC - production run - save 1 in 100 samples to remove autocorrelation, save 15 stes giving 15*200 samples (200 walkers by default)
    res = fitter.sample(15, nthin=100, random_state=1, pool=-1)
    
    print('\n\n\n')
    print(global_objective)
    print('Reduced chiSq = %f' %(global_objective.chisqr()/(3*len(data_d2O_hlip)-4)))
    #print('\nMCMC Parameters:\nd1=%s\nd2=%s' %(global_objective.thickness_tails, global_objective.thickness_heads))
    #print(global_objective.logpost())
    
    # check headSolvFrac is correct 
    calc_headSolv = mainCalcSolvFrac(modelNum=0,t_thick=thickness_tails.value,h_thick=thickness_heads.value)
    print('refnx headSolv = %f\ntrue headSolv = %f' %(head_solvent.value, calc_headSolv))
    print('Values match = %s' %bool(head_solvent.value == calc_headSolv))
    
    # plot spread of the data 
    global_objective.plot(samples=300)
    plt.yscale('log')
    plt.yscale('log')
    plt.xlabel('Q / $\AA^{-1}$')
    plt.ylabel('Reflectivity')
    plt.legend()
    
    # plot corner plot to show covariance between parameters 
    global_objective.corner()

    plt.savefig('../output/corner.png',
        format='png',
        dpi=400,
        bbox_inches='tight')    
    
    
    # plot SLD
    plotSLD = False
    if plotSLD == True:
        s_List = [s_d2O_hlip,s_acmw_dlip,s_d2o_dlip]
        fig, ax = plt.subplots()
        ax.plot(*s_List, label='')
        ax.set_ylabel("$\\rho$ / $10^{-6} \AA^{-2}$")
        ax.set_xlabel("z / $\AA$")
        ax.legend()    
        # s_d2O_hlip.plot(samples=300)
        # s_acmw_dlip.plot(samples=300)
        # s_d2o_dlip.plot(samples=300)

