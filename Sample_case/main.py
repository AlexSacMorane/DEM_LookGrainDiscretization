# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This is the main file.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import os
import shutil
from datetime import datetime
from pathlib import Path
import pickle
import matplotlib.pyplot as plt

#Own function and class
import Create_IC
import Create_IC_Polygonal
import Report
import User

#-------------------------------------------------------------------------------
#User
#-------------------------------------------------------------------------------

grain_discretization_L = [20,30,40]

#-------------------------------------------------------------------------------
#Plan simulation
#-------------------------------------------------------------------------------

for grain_discretization in grain_discretization_L :
    if not Path('ICs/'+str(grain_discretization)).exists():
        os.mkdir('ICs/'+str(grain_discretization))

simulation_report = Report.Report('Report',datetime.now())

#get data
dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations = User.All_parameters()

#find name
i_run = 1
folderpath = Path('ICs/'+dict_algorithm['template_simulation_name']+str(i_run)+'_dict_ic')
while folderpath.exists():
    i_run = i_run + 1
    folderpath = Path('ICs/'+dict_algorithm['template_simulation_name']+str(i_run)+'_dict_ic')
dict_algorithm['name_folder'] = dict_algorithm['template_simulation_name']+str(i_run)

#initial ic
Create_IC.LG_tempo(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report)

#save
outfile = open('ICs/'+dict_algorithm['name_folder']+'_dict_ic','wb')
pickle.dump(dict_ic,outfile)
outfile.close()

#prepare plot
plt.figure(1,figsize=(16,9))

#work on discretization
for grain_discretization in grain_discretization_L :
    dict_geometry['grain_discretization'] = grain_discretization

    #load perfect sphere data
    infile = open('ICs/'+dict_algorithm['name_folder']+'_dict_ic','rb')
    dict_ic = pickle.load(infile,encoding ='byte')
    infile.close()

    #convert into discrete grains
    dict_ic_discrete = Create_IC_Polygonal.Discretize_Grains(dict_ic, grain_discretization)

    #load discrete grains
    Create_IC_Polygonal.DEM_loading(dict_algorithm, dict_ic_discrete, dict_material, dict_sample, dict_sollicitations, simulation_report)

    #plot
    plt.subplot(221)
    plt.plot(dict_ic_discrete['Ecin_tracker'], label = 'n border = '+str(grain_discretization))
    plt.subplot(222)
    plt.plot(dict_ic_discrete['k0_tracker'], label = 'n border = '+str(grain_discretization))
    plt.subplot(223)
    plt.plot(dict_ic_discrete['Ymax_tracker'], label = 'n border = '+str(grain_discretization))
    plt.subplot(224)
    plt.plot(dict_ic_discrete['Fv_tracker'], label = 'n border = '+str(grain_discretization))

    #save
    outfile = open('ICs/'+str(grain_discretization)+'/'+dict_algorithm['name_folder']+'_dict_ic','wb')
    pickle.dump(dict_ic_discrete,outfile)
    outfile.close()

    #Interface user
    print()
    print('Loading discretization '+str(grain_discretization)+' done')
    print()

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

plt.savefig('ICs/'+dict_algorithm['name_folder']+'.png')
plt.close(1)
