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
import Report
import User


#-------------------------------------------------------------------------------
#Plan simulation
#-------------------------------------------------------------------------------

simulation_report = Report.Report('Report',datetime.now())

#get data
dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations = User.All_parameters()

if not Path('ICs/'+str(dict_geometry['grain_discretisation'])).exists():
    os.mkdir('ICs/'+str(dict_geometry['grain_discretisation']))

#find name
i_run = 1
folderpath = Path('ICs/'+str(dict_geometry['grain_discretisation'])+'/'+dict_algorithm['template_simulation_name']+str(i_run)+'_dict_ic')
while folderpath.exists():
    i_run = i_run + 1
    folderpath = Path('ICs/'+str(dict_geometry['grain_discretisation'])+'/'+dict_algorithm['template_simulation_name']+str(i_run)+'_dict_ic')
dict_algorithm['name_folder'] = dict_algorithm['template_simulation_name']+str(i_run)

#-------------------------------------------------------------------------------
#Initial conditions
#-------------------------------------------------------------------------------

Create_IC.LG_tempo(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report)

#save
outfile = open('ICs/'+str(dict_geometry['grain_discretisation'])+'/'+dict_algorithm['name_folder']+'_dict_ic','wb')
pickle.dump(dict_ic,outfile)
outfile.close()

#plot
plt.close(1)
plt.figure(1,figsize=(16,9))

plt.subplot(221)
plt.plot(dict_ic['Ecin_tracker'])
plt.title('Kinetic energy')

plt.subplot(222)
plt.plot(dict_ic['Ymax_tracker'])
plt.title('Upper wall position')

plt.subplot(223)
plt.plot(dict_ic['Fv_tracker'])
plt.title("Force on upper wall")

plt.subplot(224)
plt.plot(dict_ic['k0_tracker'])
plt.ylim(0,1)
plt.title("k0")


plt.savefig('Trackers.png')
