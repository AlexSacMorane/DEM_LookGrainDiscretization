# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This is the file where the user can change the different parameters for the simulation.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import math
import numpy as np

#-------------------------------------------------------------------------------
#User
#-------------------------------------------------------------------------------

def All_parameters():
    """
    This function is called in main() to have all the parameters needed in the simulation

        Input :
            Nothing
        Output :
            an algorithm dictionnary (a dict)
            a geometry dictionnary (a dict)
            an initial condition dictionnary (a dict)
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
            a sollicitations dictionnary (a dict)
    """
    #---------------------------------------------------------------------------
    #Geometric parameters

    N_grain = 300 #total number of grains

    #Disk
    R_mean = 350 #µm radius to compute the grain distribution. Then recomputed
    L_R = [1.2*R_mean,1.1*R_mean,0.9*R_mean,0.8*R_mean] #from larger to smaller
    L_percentage_R = [1/6,1/3,1/3,1/6] #distribution of the different radius
    #Recompute the mean radius
    R_mean = 0
    for i in range(len(L_R)):
        R_mean = R_mean + L_R[i]*L_percentage_R[i]

    #Compute number of grain (square or disk)
    N_grain_disk = N_grain

    #write dict
    dict_geometry = {
    'N_grain' : N_grain,
    'N_grain_disk' : N_grain_disk,
    'R_mean' : R_mean,
    'L_R' : L_R,
    'L_percentage_R' : L_percentage_R
    }

    #---------------------------------------------------------------------------
    #Material parameters

    Y = 70*(10**9)*(10**6)*(10**(-12)) #Young Modulus µN/µm2
    nu = 0.3 #Poisson's ratio
    rho = 2500*10**(-6*3) #density kg/µm3
    rho_surf_disk = 4/3*rho*R_mean #kg/µm2
    mu_friction_gg = 0.5 #grain-grain
    mu_friction_gw = 0 #grain-wall
    coeff_restitution = 0.2 #1 is perfect elastic
    # PF parameters
    M_pf = 1 # mobility
    kc_pf = 3 #graident coefficient

    #write dict
    dict_material = {
    'Y' : Y,
    'nu' : nu,
    'rho' : rho,
    'rho_surf' : rho_surf_disk,
    'mu_friction_gg' : mu_friction_gg,
    'mu_friction_gw' : mu_friction_gw,
    'coeff_restitution' : coeff_restitution,
    'M_pf' : M_pf,
    'kc_pf' : kc_pf
    }

    #---------------------------------------------------------------------------
    #Sample definition

    Lenght_mean = R_mean*N_grain_disk /N_grain #mean characteristic lenght

    #Box définition
    x_box_min = 0 #µm
    x_box_max = 2*Lenght_mean*math.sqrt(N_grain/0.6) #µm 0.6 from Santamarina, 2014 to avoid boundaries effect
    y_box_min = 0 #µm

    #write dict
    dict_sample = {
    'x_box_min' : x_box_min,
    'x_box_max' : x_box_max,
    'y_box_min' : y_box_min
    }

    #---------------------------------------------------------------------------
    #Algorithm parameters

    #Phase field
    dt_PF = 0.01 #s time step during MOOSE simulation
    n_t_PF = 10 #number of iterations PF-DEM
    factor_distribution_etai = 1.5 #margin to distribute etai
    MovePF_selector = 'DeconstructRebuild' #Move PF
    n_local = 40 #number of node inside local PF simulation
    dx_local = 2*min(dict_geometry['L_R'])/(n_local-1)
    dy_local = 2*min(dict_geometry['L_R'])/(n_local-1)
    #add into material dict from this data
    w = 4*math.sqrt(dx_local**2+dy_local**2)
    double_well_height = 10*dict_material['kc_pf']/w/w
    dict_material['w'] = w
    dict_material['double_well_height'] = double_well_height

    #DEM parameters
    dt_DEM_crit = math.pi*min(L_R)/(0.16*nu+0.88)*math.sqrt(rho*(2+2*nu)/Y) #s critical time step from O'Sullivan 2011

    dt_DEM = dt_DEM_crit/8 #s time step during DEM simulation
    factor_neighborhood = 1.5 #margin to detect a grain into a neighborhood
    i_update_neighborhoods = 100 #the frequency of the update of the neighborhood of the grains and the walls
    Spring_type = 'Ponctual' #Kind of contact
    #Stop criteria of the DEM
    i_DEM_stop = 3000 #maximum iteration for one DEM simulation
    Ecin_ratio = 0.0002
    n_window_stop = 50
    dk0_stop = 0.05
    dy_box_max_stop = 0.5

    #PF-DEM
    n_t_PFDEM = 10 #number of cycle PF-DEM

    #Number of processor
    np_proc = 4

    #Debugging
    Debug = True #plot configuration before and after DEM simulation
    Debug_DEM = False #plot configuration inside DEM
    i_print_plot = 200 #frenquency of the print and plot (if Debug_DEM) in DEM step
    clean_memory = True #delete Data, Input, Output at the end of the simulation
    SaveData = True #save simulation
    main_folder_name = 'Data_Compa' #where data are saved
    template_simulation_name = 'Run_' #template of the simulation name

    #write dict
    dict_algorithm = {
    'dt_PF' : dt_PF,
    'n_t_PF' : n_t_PF,
    'dt_DEM_crit' : dt_DEM_crit,
    'n_local' : n_local,
    'dx_local' : dx_local,
    'dy_local' : dy_local,
    'dt_DEM' : dt_DEM,
    'i_update_neighborhoods': i_update_neighborhoods,
    'i_DEM_stop' : i_DEM_stop,
    'Ecin_ratio' : Ecin_ratio,
    'n_window_stop' : n_window_stop,
    'dk0_stop' : dk0_stop,
    'dy_box_max_stop' : dy_box_max_stop,
    'n_t_PFDEM' : n_t_PFDEM,
    'MovePF_selector' : MovePF_selector,
    'Spring_type' : Spring_type,
    'np_proc' : np_proc,
    'Debug' : Debug,
    'Debug_DEM' : Debug_DEM,
    'SaveData' : SaveData,
    'main_folder_name' : main_folder_name,
    'template_simulation_name' : template_simulation_name,
    'i_print_plot' : i_print_plot,
    'factor_neighborhood' : factor_neighborhood,
    'factor_distribution_etai' : factor_distribution_etai,
    'clean_memory' : clean_memory
    }

    #---------------------------------------------------------------------------
    #Initial condition parameters

    n_generation = 2 #number of grains generation
    factor_ymax_box = 1.8 #margin to generate grains
    N_test_max = 5000 # maximum number of tries to generate a grain without overlap
    i_DEM_stop_IC = 4000 #stop criteria for DEM during IC
    Debug_DEM_IC = False #plot configuration inside DEM during IC
    i_print_plot_IC = 400 #frenquency of the print and plot (if Debug_DEM_IC) for IC
    dt_DEM_IC = dt_DEM_crit/6 #s time step during IC
    Ecin_ratio_IC = 0.0005
    factor_neighborhood_IC = 1.5 #margin to detect a grain into a neighborhood
    i_update_neighborhoods_gen = 100 #the frequency of the update of the neighborhood of the grains and the walls during IC generations
    i_update_neighborhoods_com = 100 #the frequency of the update of the neighborhood of the grains and the walls during IC combination

    #write dict
    dict_ic = {
    'n_generation' : n_generation,
    'i_update_neighborhoods_gen': i_update_neighborhoods_gen,
    'i_update_neighborhoods_com': i_update_neighborhoods_com,
    'factor_ymax_box' : factor_ymax_box,
    'i_DEM_stop_IC' : i_DEM_stop_IC,
    'Debug_DEM' : Debug_DEM_IC,
    'dt_DEM_IC' : dt_DEM_IC,
    'Ecin_ratio_IC' : Ecin_ratio_IC,
    'i_print_plot_IC' : i_print_plot_IC,
    'factor_neighborhood_IC' : factor_neighborhood_IC,
    'N_test_max' : N_test_max
    }

    #---------------------------------------------------------------------------
    #External sollicitations

    Vertical_Confinement_Linear_Force = Y*2*R_mean/1000 #µN/µm used to compute the Vertical_Confinement_Force
    Vertical_Confinement_Force = Vertical_Confinement_Linear_Force*(x_box_max-x_box_min) #µN
    gravity = 0 #µm/s2

    #Add energy to dissolved grain
    Dissolution_Energy = 0.2

    #write dict
    dict_sollicitations = {
    'Dissolution_Energy' : Dissolution_Energy,
    'Vertical_Confinement_Force' : Vertical_Confinement_Force,
    'gravity' : gravity
    }

    #---------------------------------------------------------------------------

    return dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations

#-------------------------------------------------------------------------------

def Criteria_StopSimulation(dict_algorithm):
    """
    Criteria to stop simulation (PF and DEM).

        Input :
            an algorithm dictionnary (a dict)
        Output :
            a Booean (True if the simulation must be stopped)
    """
    Criteria_Verified = False
    if dict_algorithm['i_PF'] >= dict_algorithm['n_t_PFDEM']:
        Criteria_Verified = True
    return Criteria_Verified
