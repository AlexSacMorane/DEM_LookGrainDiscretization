# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

The goal of this file is to add a discrete study from an already initial condition.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import os
import shutil
from datetime import datetime
from pathlib import Path
import pickle
import glob
import matplotlib.pyplot as plt

#Own function and class
import Create_IC
import Create_IC_Polygonal
import Report
import User

#-------------------------------------------------------------------------------
#User
#-------------------------------------------------------------------------------

grain_discretization = 60 #new discretization to compute
name_run = 'Run_1' #name of the compute to rework

#-------------------------------------------------------------------------------
#Plan simulation
#-------------------------------------------------------------------------------

if not Path('ICs/'+str(grain_discretization)).exists():
    os.mkdir('ICs/'+str(grain_discretization))
simulation_report = Report.Report('Report',datetime.now())

#-------------------------------------------------------------------------------
#Main
#-------------------------------------------------------------------------------

#load perfect sphere data
infile = open('ICs/'+name_run+'_dict_ic','rb')
dict_ic = pickle.load(infile,encoding ='byte')
infile.close()

#convert into discrete grains
dict_ic_discrete = Create_IC_Polygonal.Discretize_Grains(dict_ic, grain_discretization)

#load discrete grains
Create_IC_Polygonal.DEM_loading(dict_algorithm, dict_ic_discrete, dict_material, dict_sample, dict_sollicitations, simulation_report)

#save
outfile = open('ICs/'+str(grain_discretization)+'/'+name_run+'_dict_ic','wb')
pickle.dump(dict_ic_discrete,outfile)
outfile.close()

#prepare plot
plt.figure(1,figsize=(16,9))

#work on discretization
for name_to_load in glob.glob('ICs/*/'+name_run+'_dict_ic') :
    dict_geometry['grain_discretization'] = grain_discretization

    #load data
    infile = open(name_to_load,'rb')
    dict_ic_discrete = pickle.load(infile,encoding ='byte')
    infile.close()

    #plot
    plt.subplot(221)
    plt.plot(dict_ic_discrete['Ecin_tracker'], label = 'n border = '+name_to_load[4:6])
    plt.subplot(222)
    plt.plot(dict_ic_discrete['k0_tracker'], label = 'n border = '+name_to_load[4:6])
    plt.subplot(223)
    plt.plot(dict_ic_discrete['Ymax_tracker'], label = 'n border = '+name_to_load[4:6])
    plt.subplot(224)
    plt.plot(dict_ic_discrete['Fv_tracker'], label = 'n border = '+name_to_load[4:6])

#close plot
plt.subplot(221)
plt.title('Ecin')
plt.legend()
plt.subplot(222)
plt.title(r'k0 = $\sigma_2$ / $\sigma_1$')
plt.legend()
plt.subplot(223)
plt.title('Upper wall position')
plt.legend()
plt.subplot(224)
plt.title('Force on upper wall')
plt.legend()

plt.savefig('ICs/'+name_run+'.png')
plt.close(1)
