import os.path
import numpy as np
import matplotlib.pyplot as plt
import scipy
import pandas as pd
import sys
import functools 


import corner 
import refnx
from refnx.dataset import ReflectDataset, Data1D
from refnx.analysis import Transform, CurveFitter, Objective, GlobalObjective, Model, Parameter, process_chain
from refnx.reflect import SLD, Slab, ReflectModel, Structure, LipidLeaflet
from scipy.optimize import NonlinearConstraint

from sldAnalysis import mainCalcSolvFrac    
from calcLipidPars import mainCalcPar
import config


def calcReducedChiSq(exp_list, global_objective):
    
    nPoints  = functools.reduce(lambda count, l: count + len(l), exp_list, 0)
    nPars    = 4 
    redChiSq = global_objective.chisqr()/(nPoints-nPars)
    print('Reduced chiSq = %f' %redChiSq)
    
    return redChiSq


def refnxOptimisation(title, file_paths, qmin, qmax, lipids, ratios):
    
    print(f"refnx: {refnx.version.version}\n"
      f"scipy: {scipy.version.version}\n"
      f"numpy: {np.version.version}\n")


#####################
# Define Parameters #
#####################    
    
    if file_paths.get("d2O_dlip") is not None: 
        d_lipid = lipids.get('d2O_dlip')
        d_ratio = ratios.get('d2O_dlip')
    elif file_paths.get("acmw_dlip") is not None:
        d_lipid = lipids.get('acmw_dlip')
        d_ratio = ratios.get('acmw_dlip')
    else: 
        d_lipid = 'NA'
        d_ratio = 'NA'

    if file_paths.get("d2O_hlip") is not None: 
        h_lipid = lipids.get('d2O_hlip')
        h_ratio = ratios.get('d2O_hlip')
    elif file_paths.get("acmw_hlip") is not None:
        h_lipid = lipids.get('acmw_hlip')
        h_ratio = ratios.get('acmw_hlip')
    else: 
        h_lipid = 'NA'
        h_ratio = 'NA'
      

    # define universal roughness 
    roughness = Parameter(3.5, name='Roughness')

    # create sld slabs 
    rho_acmw = Parameter(0.00, bounds=(0.0,0.1), vary=False, name='SLD_acmw')
    rho_d2o  = Parameter(6.36, bounds=(6.3,6.4), vary=False, name='SLD_d2o')
    
    acmw = SLD(rho_acmw, name='acmw')
    d2o  = SLD(rho_d2o,name='d2o')
    air  = SLD(0.0, name='air')
    
    # SLD values    
    calculatedSLD_hlip = mainCalcPar("SLD",h_lipid,h_ratio,1,1)
    calculatedSLD_dlip = mainCalcPar("SLD",d_lipid,d_ratio,1,1)

    rho_tails_h = Parameter(calculatedSLD_hlip.get('tails'), name='SLD1_h')
    rho_tails_d = Parameter(calculatedSLD_dlip.get('tails'), name='SLD1_d')
    rho_heads   = Parameter(calculatedSLD_hlip.get('head'), name='SLD2')
    
    # SL values 
    calculatedSL_hlip = mainCalcPar("SL",h_lipid,h_ratio,1,1)
    calculatedSL_dlip = mainCalcPar("SL",d_lipid,d_ratio,1,1)
    
    b_tails_h = Parameter(calculatedSL_hlip.get('tails'), name='SL1_h')
    b_tails_d = Parameter(calculatedSL_dlip.get('tails'), name='SL1_d')
    b_heads   = Parameter(calculatedSL_hlip.get("head"), name='SL2')
    
    # SLD objects for tails and heads 
    tails_h = SLD(rho_tails_h, name='tails_h')
    tails_d = SLD(rho_tails_d, name='tails_d')
    heads   = SLD(rho_heads, name='heads')
    
    # define tail and head thicknesses 
    thickness_tails = Parameter(11.388, bounds=(10,35), vary=True, name='Tails Thickness')
    thickness_heads = Parameter(6.0, bounds=(5,15), vary=config.varyHead, name='Heads Thickness')
    
    # constrain head_solvent = 1-headVolFrac = 1-(rho1*d1*SL2)/(rho2*d2*SL1)
    head_solvent = Parameter(0.0, name='Solv2')
    head_solvent.constraint = 1-((rho_tails_d*thickness_tails*b_heads)/(rho_heads*thickness_heads*b_tails_d))
    tail_solvent = Parameter(0, name='Solv1')
  
    # create the tails and head slabs
    tails_layer_h = tails_h(thickness_tails,roughness,tail_solvent) 
    tails_layer_d = tails_d(thickness_tails,roughness,tail_solvent)  
    heads_layer   = heads(thickness_heads,roughness,head_solvent)  
    
    # create drug layer, SLD depends on buffer (acmw/d2o)
    rho_drug_h     = Parameter(3.67, name='SLD3_h')
    rho_drug_d     = Parameter(4.46, name='SLD3_d') 
    thickness_drug = Parameter(20, bounds=(15,25), vary=True, name='Drug Thickness')
    drug_solvent   = Parameter(0.8675, bounds=(0,1), vary=True, name='Solv3')
    
    drug_h       = SLD(rho_drug_h, name='drug_h')
    drug_d       = SLD(rho_drug_d, name='drug_d')

    drug_layer_h = drug_h(thickness_drug,roughness,drug_solvent)
    drug_layer_d = drug_d(thickness_drug,roughness,drug_solvent)
   

####################
# Build Structures #
####################

    # get expData, create structures, define model and set background, 
    # calc and store objective
    # order: '0: hMC3-D2O-polyA-PBS','1: dMC3-ACMW-polyA-PBS','2: dMC3-D2O-polyA-PBS','3: hMC3-ACMW-polyA-PBS'
    exp_list = []  
    obj_list = []
    if file_paths.get("d2O_dlip") is not None:
 
        # get experimental data 
        data_d2o_dlip  = ReflectDataset(file_paths.get("d2O_dlip")  )
        exp_list.append(data_d2o_dlip)
        
        # build structure with/out 3rd layer 
        if config.withDrugLayer == True: 
            s_d2o_dlip  = air | tails_layer_d | heads_layer | drug_layer_d | d2o(0, roughness)

        else: 
            s_d2o_dlip  = air | tails_layer_d | heads_layer | d2o(0, roughness)
        
        # build model from structure and set background 
        model_d2o_dlip = ReflectModel(s_d2o_dlip)
        model_d2o_dlip.bkg.setp(5e-7,vary=False, bounds=(4.5e-7, 5.5e-7))
        
        # set & store objective from model and experimental data 
        obj_d2o_dlip  = Objective(model_d2o_dlip, data_d2o_dlip)
        obj_list.append(obj_d2o_dlip)
        
        
    if file_paths.get("d2O_hlip") is not None:
          
        # get experimental data 
        data_d2O_hlip  = ReflectDataset(file_paths.get("d2O_hlip"))
        exp_list.append(data_d2O_hlip)
        
        # build structure with/out 3rd layer 
        if config.withDrugLayer == True: 
            s_d2O_hlip  = air | tails_layer_h | heads_layer | drug_layer_d | d2o(0, roughness)
        else: 
            s_d2O_hlip  = air | tails_layer_h | heads_layer | d2o(0, roughness)
      
        # build model from structure and set background 
        model_d2O_hlip  = ReflectModel(s_d2O_hlip)
        model_d2O_hlip.bkg.setp(5e-7,vary=False, bounds=(4.5e-7, 5.5e-7))
        
        # set & store objective from model and experimental data 
        obj_d2O_hlip  = Objective(model_d2O_hlip, data_d2O_hlip)
        obj_list.append(obj_d2O_hlip)
    
    
    if file_paths.get("acmw_dlip") is not None:
        
        # get experimental data 
        data_acmw_dlip = ReflectDataset(file_paths.get("acmw_dlip"))
        exp_list.append(data_acmw_dlip)
        
        # build structure with/out 3rd layer 
        if config.withDrugLayer == True: 
            s_acmw_dlip = air | tails_layer_d | heads_layer | drug_layer_h | acmw(0, roughness)
        else: 
            s_acmw_dlip  = air | tails_layer_d | heads_layer | acmw(0, roughness)

        # build model from structure and set background 
        model_acmw_dlip = ReflectModel(s_acmw_dlip)
        model_acmw_dlip.bkg.setp(5e-7,vary=False, bounds=(4.5e-7, 5.5e-7))
        
        # set & store objective from model and experimental data 
        obj_acmw_dlip = Objective(model_acmw_dlip, data_acmw_dlip)
        obj_list.append(obj_acmw_dlip)
        
    if file_paths.get("acmw_hlip") is not None:
        
        # get experimental data 
        data_acmw_hlip = ReflectDataset(file_paths.get("acmw_hlip"))
        exp_list.append(data_acmw_hlip)
        
        # build structure with/out 3rd layer 
        if config.withDrugLayer == True:
            s_acmw_hlip = air | tails_layer_h | heads_layer | drug_layer_h | acmw(0, roughness)
        else: 
            s_acmw_hlip = air | tails_layer_h | heads_layer | acmw(0, roughness)

        # build model from structure and set background 
        model_acmw_hlip = ReflectModel(s_acmw_hlip)
        model_acmw_hlip.bkg.setp(5e-7,vary=False, bounds=(4.5e-7, 5.5e-7))
        
        # set & store objective from model and experimental data 
        obj_acmw_hlip = Objective(model_acmw_hlip, data_acmw_hlip)
        obj_list.append(obj_acmw_hlip)


##########
# Do Fit #
##########
 
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
    
    # only print global objective if MCMC is not going to replace it 
    if config.doMCMC == False:
        print(global_objective)
        
    # calculate reduced chi square value 
    redChiSq = calcReducedChiSq(exp_list, global_objective)
    
    # plot fits 
    if config.plotObjective == True:
        fig, ax = global_objective.plot()
        
        plt.yscale('log')
        plt.xlabel('Q / $\AA^{-1}$', fontsize=18, fontweight='bold')
        plt.ylabel('R', fontsize=18, fontweight='bold')
        plt.legend(prop={'size': 16, 'weight':'bold'}, frameon = False, loc='best')
          
        plt.tick_params(axis='both',which='major', labelsize=13, size=8, width=2, direction='in', top='on')
        plt.tick_params(axis='both',which='minor', labelsize=13, size=4, width=2, direction='in', right='on')
       
        plt.savefig('../output/globalObj-' +title+'.png',
            format='png',
            dpi=400,
            bbox_inches='tight') 
        plt.show()
    
    
    
    ### Error: can't pickle local object 'refnxOptimisation.<locals>.LogpExtra'
    
    # class LogpExtra(object):
    #     def __init__(self, pars):
    #         # we'll store the parameters and objective in this object
    #         # this will be necessary for pickling in the future
    #         self.pars = pars
    
    #     def __call__(self, model, data):
    #         if float(self.pars[0] - self.pars[1]) > 0:
    #             print('d1-d2>0')
    #             return 0
    #         else: print('not satisfied')
    #         return -np.inf
        
        
    # pars = (thickness_tails, thickness_heads) 
    # lpe = LogpExtra(pars)
    
    # set the log_extra attribute of the Objective with our extra log-probability term.
    #global_objective.logp_extra = lpe
    
 
#################
# MCMC Sampling #
#################       
    
    if config.doMCMC == True:
        
        # MCMC - burn first 400 as initial chain might not be representative of equilibrated system  
        fitter.sample(steps=400, random_state=1, pool=-1)
        fitter.sampler.reset()
        
        # MCMC - production run - save 1 in 100 samples to remove autocorrelation, save 15 stes giving 15*200 samples (200 walkers by default)
        res = fitter.sample(15, nthin=100, random_state=1, pool=-1)
        
        print('\n\n\n')
        print(global_objective)
        redChiSq = calcReducedChiSq(exp_list, global_objective)

        #print(global_objective.logpost())
        
        # check headSolvFrac is correct 
        if config.compareSolv2 == True:
            calc_headSolv = mainCalcSolvFrac(modelNum=0,t_thick=thickness_tails.value,h_thick=thickness_heads.value)
            print('refnx headSolv = %f\ntrue headSolv = %f' %(head_solvent.value, calc_headSolv))
        
        
        # plot spread of the data 
        if config.plotGlobalObjSpread == True:
            global_objective.plot(samples=300)
            plt.yscale('log')
            plt.yscale('log')
            plt.xlabel('Q / $\AA^{-1}$')
            plt.ylabel('Reflectivity')
            plt.legend()
            plt.savefig('../output/GlobalObjSpread-' +title+'.png',
                format='png',
                dpi=400,
                bbox_inches='tight')
        
        # plot corner plot to show covariance between parameters 
        if config.plotCorner == True:
            global_objective.corner()
            plt.savefig('../output/corner-' +title+'.png',
                format='png',
                dpi=400,
                bbox_inches='tight')    
        
        # plot SLD
        if config.plotSLD == True:
            s_List = [s_d2O_hlip,s_acmw_dlip,s_d2o_dlip]
            fig, ax = plt.subplots()
            ax.plot(*s_List, label='')
            ax.set_ylabel("$\\rho$ / $10^{-6} \AA^{-2}$")
            ax.set_xlabel("z / $\AA$")
            ax.legend()    
            # s_d2O_hlip.plot(samples=300)
            # s_acmw_dlip.plot(samples=300)
            # s_d2o_dlip.plot(samples=300)
    

#######################
# Generate Model Data #
#######################
    
    # generate model datasets for plotting 
    q       = list(np.linspace(qmin, qmax, 1001))
    modelQ  = {}
    modelNR = {}
    if file_paths.get("d2O_hlip") is not None:
        modelQ[len(modelQ)]   = q
        modelNR[len(modelNR)] = list(model_d2O_hlip(q))
    
    if file_paths.get("acmw_dlip") is not None:
        modelQ[len(modelQ)]   = q
        modelNR[len(modelNR)] = list(model_acmw_dlip(q))
        
    if file_paths.get("d2O_dlip") is not None:
        modelQ[len(modelQ)]   = q
        modelNR[len(modelNR)] = list(model_d2o_dlip(q))
        
    if file_paths.get("acmw_hlip") is not None:
        qmax = 0.05
        q = list(np.linspace(qmin, qmax, 1001))
        modelQ[len(modelQ)]   = q
        modelNR[len(modelNR)] = list(model_acmw_hlip(q))
        
        
    return modelQ, modelNR, global_objective, redChiSq
    

