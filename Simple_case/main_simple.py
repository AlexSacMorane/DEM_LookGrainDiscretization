# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 09:58:04 2023

@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

The goal of this file is to see the influence of the grain discretisation on the DEM result.
"""
#------------------------------------------------------------------------------
# Librairies
#------------------------------------------------------------------------------

import math 
import random
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

#------------------------------------------------------------------------------
# Class
#------------------------------------------------------------------------------

class Grain :
    
    def __init__(self, Center, Radius, N_border):
        """
        Define the grain class.
        
            Input :
                a center (a 2 x 1 numpy array)
                a radius (a float)
                a discretisation (a int)
            Output :
                itself (a grain)
        """
        self.center = Center
        L_border = []
        L_border_x = []
        L_border_y = []
        L_theta = []
        L_r = [] 
        for i in range(N_border):
            theta = 2*math.pi*i/N_border
            L_border.append(Center + Radius*np.array([math.cos(theta), math.sin(theta)]))
            L_border_x.append(Center[0] + Radius*math.cos(theta))
            L_border_y.append(Center[1] + Radius*math.sin(theta))
            L_theta.append(theta)
            L_r.append(Radius)
        L_border.append(L_border[0])
        L_border_x.append(L_border_x[0])
        L_border_y.append(L_border_y[0])
        self.l_border = L_border
        self.l_border_x = L_border_x
        self.l_border_y = L_border_y
        self.l_theta = L_theta
        self.l_r = L_r
        self.r_max = max(L_r)
    
    #------------------------------------------------------------------------------
    
    def random_rotation(self):
        """
        Rotate a grain with a random angle.
            
            Input : 
                itself (a grain)
            Output :
                Nothing, but the grain is rotate
        """
        angle = random.uniform(0,2*math.pi)
        for i_p in range(len(self.l_border)):
            M_rot = np.array([[math.cos(angle), -math.sin(angle)],
                              [math.sin(angle),  math.cos(angle)]])
            u = self.l_border[i_p] - self.center
            u = np.dot(M_rot,u)
            self.l_border[i_p] = u + self.center
            self.l_border_x[i_p] = u[0] + self.center[0]
            self.l_border_y[i_p] = u[1] + self.center[1]
            if i_p != len(self.l_border) - 1 :
                theta = self.l_theta[i_p]
                theta = theta + angle
                if theta >= 2*math.pi :
                    theta = theta - 2*math.pi
                self.l_theta[i_p] = theta
        
#------------------------------------------------------------------------------
# Functions
#------------------------------------------------------------------------------

def Compute_Overlap(g1, g2):
    """
    Compute the overlap between two grains.
    
        Input: 
            two grains (grains)
        Output :
            an overlap (a float)
    """
    #compute angle between grains
    g1_to_g2 = g2.center - g1.center
    if g1_to_g2[1] >= 0 :
        angle_g1_to_g2 = math.acos(g1_to_g2[0]/np.linalg.norm(g1_to_g2))
        angle_g2_to_g1 = angle_g1_to_g2 + math.pi
    else :
        angle_g1_to_g2 = math.pi + math.acos(-g1_to_g2[0]/np.linalg.norm(g1_to_g2))
        angle_g2_to_g1 = angle_g1_to_g2 - math.pi

    #extract
    L_i_vertices_1 = extract_vertices(g1, angle_g1_to_g2)
    L_i_vertices_2 = extract_vertices(g2, angle_g2_to_g1)

    #looking for the nearest nodes
    d_virtual = max(g1.r_max,g2.r_max)
    ij_min = [0,0]
    d_ij_min = 100*d_virtual #Large
    for i in L_i_vertices_1:
        for j in L_i_vertices_2:
            d_ij = np.linalg.norm(g2.l_border[:-1][j]-g1.l_border[:-1][i]+d_virtual*(g2.center-g1.center)/np.linalg.norm(g2.center-g1.center))
            if d_ij < d_ij_min :
                d_ij_min = d_ij
                ij_min = [i,j]

    #-----------------------------------------------------------------------------
    #Computing CP
    #-----------------------------------------------------------------------------

    M = (g1.l_border[:-1][ij_min[0]]+g2.l_border[:-1][ij_min[1]])/2
    # 5 candidates for CP
    N = np.array([g1.l_border[:-1][ij_min[0]][0] - g2.l_border[:-1][ij_min[1]][0],
                  g1.l_border[:-1][ij_min[0]][1] - g2.l_border[:-1][ij_min[1]][1]])
    N = N/np.linalg.norm(N)
    PB = np.array([-N[1] ,N[0]])
    PB = PB/np.linalg.norm(PB)

    #candidats from grain 1
    if ij_min[0] <len(g1.l_border[:-1]) - 1:
        M1 = g1.l_border[:-1][ij_min[0]+1]-g1.l_border[:-1][ij_min[0]]
    else :
        M1 = g1.l_border[:-1][0]-g1.l_border[:-1][ij_min[0]]
    M1 = M1/np.linalg.norm(M1)
    M3 = g1.l_border[:-1][ij_min[0]-1]-g1.l_border[:-1][ij_min[0]]
    M3 = M3/np.linalg.norm(M3)
    #reorganize the candidats
    if np.dot(M1,PB) < 0:
        Mtempo = M1.copy()
        M1 = M3.copy()
        M3 = Mtempo.copy()

    #candidats from grain 2
    if ij_min[1] <len(g2.l_border[:-1]) - 1:
        M2 = g2.l_border[:-1][ij_min[1]+1]-g2.l_border[:-1][ij_min[1]]
    else :
        M2 = g2.l_border[:-1][0]-g2.l_border[:-1][ij_min[1]]
    M2 = M2/np.linalg.norm(M2)
    M4 = g2.l_border[:-1][ij_min[1]-1]-g2.l_border[:-1][ij_min[1]]
    M4 = M4/np.linalg.norm(M4)
    #reorganize the candidats
    if np.dot(M2,PB) < 0:
      Mtempo = M2.copy()
      M2 = M4.copy()
      M4 = Mtempo.copy()

    #compute the different angles
    theta_PB = math.pi/2
    theta_M1 =  math.acos(np.dot(M1,N))
    theta_M2 =  math.acos(np.dot(M2,N))
    theta_M3 = -math.acos(np.dot(M3,N))
    theta_M4 = -math.acos(np.dot(M4,N))

    #find the PC
    if theta_M2 < theta_PB and theta_PB < theta_M1\
       and theta_M3 < -theta_PB and -theta_PB < theta_M4:
       PC = PB
    else:
      L_Mi = [M1,M2,M3,M4]
      L_theta_Mi_PB=[theta_M1-theta_PB, theta_PB-theta_M2, -theta_M3-theta_PB, theta_PB+theta_M4]
      PC = L_Mi[L_theta_Mi_PB.index(min(L_theta_Mi_PB))]

    #-----------------------------------------------------------------------------
    # Compute the normal and tangential planes
    #-----------------------------------------------------------------------------

    PC_normal = np.array([PC[1],-PC[0]])
    if np.dot(PC_normal,(g2.center-g1.center)/np.linalg.norm(g2.center-g1.center))<0 :
        PC_normal = np.array([-PC[1],PC[0]])

    #-----------------------------------------------------------------------------
    # Compute the overlap
    #-----------------------------------------------------------------------------

    d_b = np.dot(M-g2.l_border[:-1][ij_min[1]],PC_normal)
    d_a = np.dot(M-g1.l_border[:-1][ij_min[0]],PC_normal)
    overlap = d_b - d_a
    
    return overlap

#-------------------------------------------------------------------------------

def extract_vertices(g, angle_g_to_other_g) :
    """
    Extract a list of indices of vertices inside a angular window.

        Input :
            a grain (a grain)
            an angle (a float)
        Output :
            a list of indices (a list)
    """
    dtheta = 3*2*math.pi/len(g.l_border)
    angle_minus = angle_g_to_other_g - dtheta/2
    if angle_minus < 0 :
        angle_minus = angle_minus + 2*math.pi
    angle_plus  = angle_g_to_other_g + dtheta/2
    if 2*math.pi <= angle_plus :
        angle_plus = angle_plus - 2*math.pi
    i_minus = find_value_in_list(g.l_theta.copy(), angle_minus)
    i_plus = i_minus + find_value_in_list(g.l_theta[i_minus:].copy() + g.l_theta[:i_minus].copy(), angle_plus)

    L_i_vertices = []
    for i in range(i_minus, i_plus+1):
        if i < len(g.l_theta):
            L_i_vertices.append(i)
        else :
            L_i_vertices.append(i-len(g.l_theta))
    return L_i_vertices

#-------------------------------------------------------------------------------

def find_value_in_list(List_to_search, value_to_find) :
    """
    Extract the index of the nearest value from a target in a list.

        Input :
            a list of float (a list)
            a target (a float)
        Output :
            an index (an int)
    """
    L_search = list(abs(np.array(List_to_search)-value_to_find))
    return L_search.index(min(L_search))

#------------------------------------------------------------------------------
# User
#------------------------------------------------------------------------------

L_n_border = [10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]
radius = 1
L_overlap = [0.05, 0.1, 0.15, 0.2]
N_ite = 100

#------------------------------------------------------------------------------
# Main
#------------------------------------------------------------------------------

L_L_1q_2g = []
L_L_3q_2g = []
L_L_1q_w = []
L_L_3q_w = []

for i_overlap in range(len(L_overlap)):
    
    overlap = L_overlap[i_overlap]*radius
    
    M_overlap_2g = np.array(np.zeros((N_ite,len(L_n_border))))
    M_overlap_w = np.array(np.zeros((N_ite,len(L_n_border))))
    
    L_1q_2g = []
    L_3q_2g = []
    L_1q_w = []
    L_3q_w = []
    
    #prepare plot
    plt.close(i_overlap)
    plt.figure(i_overlap, figsize = (16,9))
    
    for i_n_border in range(len(L_n_border)) : 
    
        n_border = L_n_border[i_n_border] 
        L_overlap_computed = []
        L_overlap_w = []
        
        for i in range(N_ite) :
            
            #----------------------------------------------------------------------
            #2 grains case
            
            #Create grains
            grain_1 = Grain(np.array([-radius + overlap/2, 0]), radius, n_border)
            grain_1.random_rotation()
            grain_2 = Grain(np.array([ radius - overlap/2, 0]), radius, n_border)
            grain_2.random_rotation()    
            
            #Compute overlap
            overlap_computed = Compute_Overlap(grain_1, grain_2)
        
            #save
            M_overlap_2g[i][i_n_border] = overlap_computed
            L_overlap_computed.append(overlap_computed)
    
            #----------------------------------------------------------------------
            #grain - wall case
            
            #Create grains
            grain_w = Grain(np.array([radius - overlap, 0]), radius, n_border)
            grain_w.random_rotation()  
            
            #Compute overlap
            overlap_w = abs(min(grain_w.l_border_x))
        
            #save
            M_overlap_w[i][i_n_border] = overlap_w
            L_overlap_w.append(overlap_w)
    
        
        L_q_2g = np.quantile(L_overlap_computed, q = [0.25,0.75])
        L_1q_2g.append((L_q_2g[0]-overlap)/overlap)
        L_3q_2g.append((L_q_2g[1]-overlap)/overlap)
        
        L_q_w = np.quantile(L_overlap_w, q = [0.25,0.75])
        L_1q_w.append((L_q_w[0]-overlap)/overlap)
        L_3q_w.append((L_q_w[1]-overlap)/overlap)

    L_L_1q_2g.append(L_1q_2g)
    L_L_3q_2g.append(L_3q_2g)
    L_L_1q_w.append(L_1q_w)
    L_L_3q_w.append(L_3q_w)

    #------------------------------------------------------------------------------
    # Plot
    #------------------------------------------------------------------------------
    
    plt.subplot(1,2,1)
    plt.title('2 grains case', fontsize = 26)
    plt.boxplot(M_overlap_2g, positions = L_n_border)
    plt.xlabel('Number of vertices', fontsize = 20)
    plt.ylabel('Overlap computed', fontsize = 20)
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    
    plt.subplot(1,2,2)
    plt.title('Grain - wall case', fontsize = 26)
    plt.boxplot(M_overlap_w, positions = L_n_border)
    plt.xlabel('Number of vertices', fontsize = 20)
    plt.ylabel('Overlap computed', fontsize = 20)
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    
    plt.suptitle('Overlap = ' + str(int(100*L_overlap[i_overlap]))+' % of radius, '+str(N_ite)+' repetitions')
    plt.savefig('Distrib_Overlap_'+str(int(100*L_overlap[i_overlap]))+'.png')
    
    if i_overlap == 1 :

        plt.suptitle('')
        plt.savefig('Distrib_Overlap_'+str(int(100*L_overlap[i_overlap]))+'_correct.png')
        
    plt.close(i_overlap)        


plt.close(len(L_overlap))
plt.figure(len(L_overlap),figsize=(16,9))

L_color = ['b', 'r', 'g', 'darkorange', 'blueviolet']

plt.subplot(121)
plt.title('2 grains case',fontsize = 26)
L_patchs = []
for i_L in range(len(L_overlap)) :
    plt.plot(L_n_border, L_L_1q_2g[i_L], color = L_color[i_L],lw = 4)
    plt.plot(L_n_border, L_L_3q_2g[i_L], color = L_color[i_L],lw = 4)
    patch = mpatches.Patch(color=L_color[i_L], label = 'Overlap = '+str(int((100*L_overlap[i_L])))+' % of radius')
    L_patchs.append(patch)
plt.legend(handles=L_patchs, fontsize = 16, loc='lower right')
plt.xlabel('Number of vertices', fontsize = 20)
plt.xticks(fontsize = 14)
plt.ylabel(r'$\Delta$ overlap / real overlap', fontsize = 20)    
plt.yticks(fontsize = 14)

plt.subplot(122)
plt.title('Grain - wall case',fontsize = 26)
L_patchs = []
for i_L in range(len(L_overlap)) :
    plt.plot(L_n_border, L_L_1q_w[i_L], color = L_color[i_L],lw = 4)
    plt.plot(L_n_border, L_L_3q_w[i_L], color = L_color[i_L],lw = 4)
    patch = mpatches.Patch(color=L_color[i_L], label = 'Overlap = '+str(int((100*L_overlap[i_L])))+' % of radius')
    L_patchs.append(patch)
plt.legend(handles=L_patchs, fontsize = 16, loc='lower right')
plt.xlabel('Number of vertices', fontsize = 20)
plt.xticks(fontsize = 14)
plt.ylabel(r'$\Delta$ overlap / real overlap', fontsize = 20)
plt.yticks(fontsize = 14)

#plt.suptitle('1st & 3rd quartiles for '+str(N_ite)+' repetitions')
plt.savefig('1_3_q.png')
plt.close(len(L_overlap))


# extra...

if False :
    grain_1 = Grain(np.array([-radius + 0.2*radius/2, 0]), radius, 10)
    grain_1_real = Grain(np.array([-radius + 0.2*radius/2, 0]), radius, 360)
    grain_1.random_rotation()
    grain_2 = Grain(np.array([ radius - 0.2*radius/2, 0]), radius, 10)
    grain_2_real = Grain(np.array([ radius - 0.2*radius/2, 0]), radius, 360)
    grain_2.random_rotation()   
    
    plt.close(1)
    plt.figure(1,figsize=(16,9))
    
    plt.plot(grain_1.l_border_x, grain_1.l_border_y,'k', lw = 6)
    plt.plot(grain_1_real.l_border_x, grain_1_real.l_border_y,'k-.', lw = 6)
    plt.plot(grain_2.l_border_x, grain_2.l_border_y,'k', lw = 6)
    plt.plot(grain_2_real.l_border_x, grain_2_real.l_border_y,'k-.', lw = 6)
    
    plt.axis('equal')
    plt.tick_params(left = False, bottom = False, labelbottom = False, labelleft = False)
    
    plt.close(2)
    plt.figure(2,figsize=(16,9))
    
    plt.plot(grain_1.l_border_x, grain_1.l_border_y,'k', lw = 6)
    plt.plot(grain_1_real.l_border_x, grain_1_real.l_border_y,'k-.', lw = 6)
    plt.plot([min(grain_1_real.l_border_x)+overlap, min(grain_1_real.l_border_x)+overlap],
             [min(grain_1_real.l_border_y), max(grain_1_real.l_border_y)], 'k', lw = 6)
    
    plt.axis('equal')
    plt.tick_params(left = False, bottom = False, labelbottom = False, labelleft = False)
