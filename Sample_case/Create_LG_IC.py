# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

The goal of this file is to define an initial configuration.
Grains are considered as a square or a disk.
We have 2 temporary classes about grains and contact."""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import numpy as np
import math
import random
import matplotlib.pyplot as plt

#-------------------------------------------------------------------------------
#Classes definition
#-------------------------------------------------------------------------------

class Grain_Tempo:
  """
  A temporary grain used to generated an initial condition.
  """

#-------------------------------------------------------------------------------

  def __init__(self, ID, Center, Lenght, n_border, dict_material, Type):
    """Defining the grain.

        Input :
            itself (a grain_tempo)
            an id (a int)
            a center coordinate (a 1 x 2 numpy array)
            a dimension, radius if it is a disk or the size if it is a square (a float)
            a material dictionnary (a dict)
            a grain type, disk or square (a float)
        Output :
            Nothing, but a temporary grain is generated (a grain_tempo)
    """
    L_border = []
    L_border_x = []
    L_border_y = []
    L_r = []
    L_theta_r = []
    if Type == 'Disk':
        #Build the border
        for i in range(n_border):
            theta = 2*math.pi*i/n_border
            p = np.array(Center) + np.array([Lenght*math.cos(theta),Lenght*math.sin(theta)])
            L_border.append(p)
            L_border_x.append(p[0])
            L_border_y.append(p[1])
            L_r.append(Lenght)
            L_theta_r.append(theta)
        L_border.append(L_border[0])
        L_border_x.append(L_border_x[0])
        L_border_y.append(L_border_y[0])
        self.radius = Lenght
        self.theta = 0
        self.rho_surf = dict_material['rho_surf_disk']
        self.surface = math.pi*Lenght**2
        self.mass = math.pi*Lenght**2*dict_material['rho_surf_disk']
        self.inertia = self.mass*Lenght**2
    if Type == 'Square':
        if False:
            #initial random rotation
            Theta = random.uniform(0,math.pi/2)
            #Build the border
            P_rot = np.array([[math.cos(Theta), -math.sin(Theta)],
                              [math.sin(Theta),  math.cos(Theta)]]) #initial rotation
            for i in range(n_border):
                angle = 2*math.pi*i/n_border #angle of the vertex
                angle_to_determine_R = angle%(math.pi/2)
                if angle_to_determine_R > math.pi/4:
                    angle_to_determine_R = angle_to_determine_R - math.pi/2
                p = np.array([Lenght/2/math.cos(angle_to_determine_R),0]) #compute the radius, considering the grain as a square
                P_rot_p = np.array([[math.cos(angle), -math.sin(angle)],
                                    [math.sin(angle),  math.cos(angle)]])
                p = np.dot(P_rot_p,p)
                p = np.dot(P_rot,p)
                p = p + np.array(Center)
                L_border.append(p)
                L_border_x.append(p[0])
                L_border_y.append(p[1])
                L_r.append(Lenght/2/math.cos(angle_to_determine_R))
                angle_r = Theta + angle
                while angle_r >= 2*math.pi:
                    angle_r = angle_r - 2*math.pi
                while angle_r < 0:
                    angle_r = angle_r + 2*math.pi
                L_theta_r.append(angle_r)
            L_border.append(L_border[0])
            L_border_x.append(L_border_x[0])
            L_border_y.append(L_border_y[0])
            self.dimension = Lenght
            self.theta = Theta
            self.rho_surf = dict_material['rho_surf_square']
            self.surface = Lenght**2
            self.mass = Lenght**2*dict_material['rho_surf_square']
            self.inertia = 1/6*self.mass*Lenght**2
        else :
            #Build the border
            for i in range(n_border):
                theta = 2*math.pi*i/n_border
                p = np.array(Center) + np.array([Lenght/2*math.cos(theta),Lenght/2*math.sin(theta)])
                L_border.append(p)
                L_border_x.append(p[0])
                L_border_y.append(p[1])
                L_r.append(Lenght/2)
                L_theta_r.append(theta)
            L_border.append(L_border[0])
            L_border_x.append(L_border_x[0])
            L_border_y.append(L_border_y[0])
            self.dimension = Lenght/2
            self.theta = 0
            self.rho_surf = dict_material['rho_surf_disk']
            self.surface = math.pi*(Lenght/2)**2
            self.mass = math.pi*(Lenght/2)**2*dict_material['rho_surf_disk']
            self.inertia = self.mass*(Lenght/2)**2

    #save
    self.id = ID
    self.type = Type
    self.center = np.array(Center)
    self.l_border = L_border
    self.l_border_x = L_border_x
    self.l_border_y = L_border_y
    self.l_theta_r = L_theta_r
    self.l_r = L_r
    self.r_min = min(L_r)
    self.r_max = max(L_r)
    self.y = dict_material['Y']
    self.nu = dict_material['nu']
    self.g = dict_material['Y']/2/(1+dict_material['nu']) #shear modulus
    self.fx = 0
    self.fy = 0
    self.v = np.array([0,0])
    self.w = 0

#-------------------------------------------------------------------------------

  def add_F(self, F, p_application):
      """
      Add a force to the grain.

        Input :
            itself (a grain_tempo)
            a force applied (a 1 x 2 numpy array)
            a application point (a 1 x 2 numpy array)
        Output :
            Nothing, but attributes are updated (three floats)
      """
      self.fx = self.fx + F[0]
      self.fy = self.fy + F[1]
      v1 = np.array([p_application[0]-self.center[0], p_application[1]-self.center[1], 0])
      v2 = np.array([F[0], F[1], 0])
      self.mz = self.mz + np.cross(v1,v2)[2]

#-------------------------------------------------------------------------------

  def init_F_control(self,g):
      """
      Initialize the force applied to the grain.

      A gravity is assumed.

        Input :
            itself (a grain_tempo)
            a gravity value (a float)
        Output :
            Nothing, but attributes concerning the force applied are initialized (three floats)
      """
      self.fx = 0
      self.fy = -g*self.mass
      self.mz = 0

#-------------------------------------------------------------------------------

  def euler_semi_implicite(self,dt_DEM,factor):
    """
    Move the grain following a semi implicit euler scheme.

        Input :
            itself (a grain_tempo)
            a time step (a float)
            a factor to limite the displacement (a float)
        Output :
            Nothing, but the grain is moved
    """
    #translation
    a_i = np.array([self.fx,self.fy])/self.mass
    self.v = self.v + a_i*dt_DEM
    if self.type == 'Disk':
        if np.linalg.norm(self.v) > self.radius*factor/dt_DEM: #limitation of the speed
            self.v = self.v * self.radius*factor/dt_DEM/np.linalg.norm(self.v)
    elif self.type == 'Square':
        if np.linalg.norm(self.v) > self.dimension*factor/dt_DEM: #limitation of the speed
            self.v = self.v * self.dimension*factor/dt_DEM/np.linalg.norm(self.v)
    for i in range(len(self.l_border)):
        self.l_border[i] = self.l_border[i] + self.v*dt_DEM
        self.l_border_x[i] = self.l_border_x[i] + self.v[0]*dt_DEM
        self.l_border_y[i] = self.l_border_y[i] + self.v[1]*dt_DEM
    self.center = self.center + self.v*dt_DEM

    #rotation
    dw_i = self.mz/self.inertia
    self.w = self.w + dw_i*dt_DEM
    self.theta = self.theta + self.w*dt_DEM
    for i_theta_r in range(len(self.l_theta_r)) :
        theta_r = self.l_theta_r[i_theta_r]
        theta_r = theta_r + self.w*dt_DEM
        while theta_r >= 2*math.pi:
            theta_r = theta_r - 2*math.pi
        while theta_r < 0 :
            theta_r = theta_r + 2*math.pi
        self.l_theta_r[i_theta_r] = theta_r
    #rigib body rotation
    for i in range(len(self.l_border)):
        p = self.l_border[i] - self.center
        Rot_Matrix = np.array([[math.cos(self.w*dt_DEM), -math.sin(self.w*dt_DEM)],
                               [math.sin(self.w*dt_DEM),  math.cos(self.w*dt_DEM)]])
        p = np.dot(Rot_Matrix,p)
        self.l_border[i] = p + self.center
        self.l_border_x[i] = p[0] + self.center[0]
        self.l_border_y[i] = p[1] + self.center[1]

#-------------------------------------------------------------------------------

class Contact_Tempo:
  """
  A temporary contact grain - grain used to generated an initial condition.
  """

#-------------------------------------------------------------------------------

  def __init__(self, ID, G1, G2, dict_material):
    """
    Defining the contact grain-grain.

        Input :
            itself (a contact_tempo)
            an id (a int)
            two grains (two grain_tempos)
            a material dictionnary (a dict)
        Output :
            Nothing, but the contact grain - grain is generated (a contact_tempo)
    """
    self.id = ID
    self.g1 = G1
    self.g2 = G2
    self.ft = 0
    self.mu = 0
    self.coeff_restitution = dict_material['coeff_restitution']
    self.tangential_old_statut = False
    self.overlap_tangential = 0

#-------------------------------------------------------------------------------

  def normal(self):
    """
    Compute the normal reaction of a contact grain-grain.

    Here a pontual spring is considered.

        Input :
            itself (a contact_tempo)
        Output :
            Nothing, but attributes are updated
    """
    #compute angle between grains
    g1_to_g2 = self.g2.center - self.g1.center
    if g1_to_g2[1] >= 0 :
        angle_g1_to_g2 = math.acos(g1_to_g2[0]/np.linalg.norm(g1_to_g2))
        angle_g2_to_g1 = angle_g1_to_g2 + math.pi
    else :
        angle_g1_to_g2 = math.pi + math.acos(-g1_to_g2[0]/np.linalg.norm(g1_to_g2))
        angle_g2_to_g1 = angle_g1_to_g2 - math.pi

    #extract
    L_i_vertices_1 = extract_vertices(self.g1, angle_g1_to_g2)
    L_i_vertices_2 = extract_vertices(self.g2, angle_g2_to_g1)

    #looking for the nearest nodes
    d_virtual = max(self.g1.r_max,self.g2.r_max)
    ij_min = [0,0]
    d_ij_min = 100*d_virtual #Large
    for i in L_i_vertices_1:
        for j in L_i_vertices_2:
            d_ij = np.linalg.norm(self.g2.l_border[:-1][j]-self.g1.l_border[:-1][i]+d_virtual*(self.g2.center-self.g1.center)/np.linalg.norm(self.g2.center-self.g1.center))
            if d_ij < d_ij_min :
                d_ij_min = d_ij
                ij_min = [i,j]
    self.ij_min = ij_min

    #-----------------------------------------------------------------------------
    #Computing CP
    #-----------------------------------------------------------------------------

    M = (self.g1.l_border[:-1][ij_min[0]]+self.g2.l_border[:-1][ij_min[1]])/2
    # 5 candidates for CP
    N = np.array([self.g1.l_border[:-1][ij_min[0]][0] - self.g2.l_border[:-1][ij_min[1]][0],
                  self.g1.l_border[:-1][ij_min[0]][1] - self.g2.l_border[:-1][ij_min[1]][1]])
    N = N/np.linalg.norm(N)
    PB = np.array([-N[1] ,N[0]])
    PB = PB/np.linalg.norm(PB)

    #candidats from grain 1
    if ij_min[0] <len(self.g1.l_border[:-1]) - 1:
        M1 = self.g1.l_border[:-1][ij_min[0]+1]-self.g1.l_border[:-1][ij_min[0]]
    else :
        M1 = self.g1.l_border[:-1][0]-self.g1.l_border[:-1][ij_min[0]]
    M1 = M1/np.linalg.norm(M1)
    M3 = self.g1.l_border[:-1][ij_min[0]-1]-self.g1.l_border[:-1][ij_min[0]]
    M3 = M3/np.linalg.norm(M3)
    #reorganize the candidats
    if np.dot(M1,PB) < 0:
        Mtempo = M1.copy()
        M1 = M3.copy()
        M3 = Mtempo.copy()

    #candidats from grain 2
    if ij_min[1] <len(self.g2.l_border[:-1]) - 1:
        M2 = self.g2.l_border[:-1][ij_min[1]+1]-self.g2.l_border[:-1][ij_min[1]]
    else :
        M2 = self.g2.l_border[:-1][0]-self.g2.l_border[:-1][ij_min[1]]
    M2 = M2/np.linalg.norm(M2)
    M4 = self.g2.l_border[:-1][ij_min[1]-1]-self.g2.l_border[:-1][ij_min[1]]
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
    if np.dot(PC_normal,(self.g2.center-self.g1.center)/np.linalg.norm(self.g2.center-self.g1.center))<0 :
        PC_normal = np.array([-PC[1],PC[0]])
    self.pc_normal = PC_normal #n12
    self.pc_tangential = np.array([-PC_normal[1],PC_normal[0]])

    #-----------------------------------------------------------------------------
    # Compute the overlap
    #-----------------------------------------------------------------------------

    d_b = np.dot(M-self.g2.l_border[:-1][ij_min[1]],PC_normal)
    d_a = np.dot(M-self.g1.l_border[:-1][ij_min[0]],PC_normal)
    overlap = d_b - d_a
    self.overlap_normal = overlap

    if overlap > 0:
    #-----------------------------------------------------------------------------
    # Compute the reaction
    #-----------------------------------------------------------------------------

        #Spring term
        Y_eq = 1/((1-self.g1.nu*self.g1.nu)/self.g1.y+(1-self.g2.nu*self.g2.nu)/self.g2.y)
        R_eq = 1/(1/np.linalg.norm(self.g1.center-self.g1.l_border[self.ij_min[0]])+1/np.linalg.norm(self.g2.center-self.g2.l_border[self.ij_min[1]]))
        k = 4/3*Y_eq*math.sqrt(R_eq)
        F_2_1_n = -k * overlap**(3/2)  #unlinear spring
        F_2_1 = F_2_1_n * PC_normal
        self.F_2_1_n = F_2_1_n
        self.Ep_n = 2/5 * k * overlap**(5/2) #-dEp/dx = F_2_1_n
        self.g1.add_F( F_2_1, self.g1.l_border[:-1][ij_min[0]])
        self.g2.add_F(-F_2_1, self.g2.l_border[:-1][ij_min[1]])

        #Damping term
        gamma = -math.log(self.coeff_restitution)/math.sqrt(math.pi**2+math.log(self.coeff_restitution)**2)
        mass_eq = self.g1.mass*self.g2.mass/(self.g1.mass+self.g2.mass)
        eta = 2 * gamma * math.sqrt(mass_eq*k)
        F_2_1_damp_n = np.dot(self.g2.v - self.g1.v,PC_normal)*eta
        F_2_1_damp = F_2_1_damp_n *PC_normal
        self.F_2_1_damp = F_2_1_damp_n
        self.g1.add_F( F_2_1_damp, self.g1.l_border[:-1][ij_min[0]])
        self.g2.add_F(-F_2_1_damp, self.g2.l_border[:-1][ij_min[1]])

    #no contact finally
    else :
        self.F_2_1_n = 0
        self.F_2_1_damp = 0
        self.Ep_n = 0

#-------------------------------------------------------------------------------

  def tangential(self,dt_DEM):
    """
    Compute the tangential reaction of a contact grain-grain.

    Here a pontual spring is considered

        Input :
            itself (a contact_tempo)
            a time step (a float)
        Output :
            Nothing, but attributes are updated
    """
    if self.overlap_normal > 0 and self.mu > 0:

        if self.tangential_old_statut:
          #if a reaction has been already computed
          #need to project the tangential reaction on the new tangential plane
          self.ft = self.ft*np.dot(self.tangential_old,self.pc_tangential)
        else:
          self.tangential_old_statut = True

        G_eq = 1/((1-self.g1.nu)/self.g1.g+(1-self.g2.nu)/self.g2.g)
        R_eq = 1/(1/np.linalg.norm(self.g1.center-self.g1.l_border[self.ij_min[0]])+1/np.linalg.norm(self.g2.center-self.g2.l_border[self.ij_min[1]]))
        kt0 = 8 * G_eq *math.sqrt(R_eq*abs(self.overlap_normal))
        kt = kt0*math.sqrt(max(1-2/3*kt0*abs(self.overlap_tangential)/self.mu/abs(self.F_2_1_n),0)) #not linear spring

        r1 = np.linalg.norm(self.g1.l_border[:-1][self.ij_min[0]] - self.g1.center) - self.overlap_normal/2
        r2 = np.linalg.norm(self.g2.l_border[:-1][self.ij_min[1]] - self.g2.center) - self.overlap_normal/2
        Delta_Us = (np.dot(self.g1.v-self.g2.v,self.pc_tangential) + r1*self.g1.w + r2*self.g2.w)*dt_DEM
        self.overlap_tangential = self.overlap_tangential + Delta_Us
        self.ft = self.ft - kt*Delta_Us
        self.tangential_old = self.pc_tangential
        if abs(self.ft) > abs(self.mu*self.F_2_1_n) or kt == 0: #Coulomb criteria
            self.ft = self.mu * abs(self.F_2_1_n) * np.sign(self.ft)

        self.g1.add_F( self.ft*self.pc_tangential, self.g1.l_border[:-1][self.ij_min[0]])
        self.g2.add_F(-self.ft*self.pc_tangential, self.g2.l_border[:-1][self.ij_min[1]])

        #Damping term
        gamma = -math.log(self.coeff_restitution)/math.sqrt(math.pi**2+math.log(self.coeff_restitution)**2)
        mass_eq = self.g1.mass*self.g2.mass/(self.g1.mass+self.g2.mass)
        eta = 2 * gamma * math.sqrt(mass_eq*kt)
        F_2_1_damp_t = -Delta_Us/dt_DEM*eta/2
        F_2_1_damp = F_2_1_damp_t *self.pc_tangential
        self.ft_damp = F_2_1_damp_t
        self.g1.add_F( F_2_1_damp, self.g1.l_border[:-1][self.ij_min[0]])
        self.g2.add_F(-F_2_1_damp, self.g2.l_border[:-1][self.ij_min[1]])

    #no contact finally
    else :
        tangential_old_statut = False
        self.overlap_tangential = 0
        self.ft = 0
        self.ft_damp = 0

#-------------------------------------------------------------------------------

class Contact_gw_Tempo:
  """
  A temporary contact grain - wall used to generated an initial condition.
  """

#-------------------------------------------------------------------------------

  def __init__(self, ID, G, dict_material, Nature, Limit, Overlap):
    """
    Defining the contact grain-wall.

        Input :
             itself (a contact_gw_tempo)
             an id (a int)
             a grain (a grain_tempo)
             a material dictionnary (a dict)
             the nature of the wall (a string)
             the coordinate of the wall (a float)
             an overlap (a float)
    """
    self.id = ID
    self.g = G
    factor = 5 #factor just to increase the stiffness
    self.k = factor*4/3*self.g.y/(1-self.g.nu*self.g.nu)*math.sqrt(self.g.r_max) #Hertz law
    self.kt = 0
    self.ft = 0
    self.limit = Limit
    self.nature = Nature
    self.mu = 0
    self.coeff_restitution = dict_material['coeff_restitution']
    self.overlap = Overlap
    self.tangential_old_statut = False
    self.overlap_tangential = 0

#-------------------------------------------------------------------------------

  def update_overlap(self,new_overlap):
    '''
    Update the overlap of a contact already created.

        Input :
            itself (a contact_gw_tempo)
            an overlap (a float)
        Output :
            Nothing, but the attribut concerning the overlap is updated (a float)
    '''
    self.overlap = new_overlap

#-------------------------------------------------------------------------------

  def  normal(self):
    """
    Compute the normal reaction of a contact grain-wall.

    Here a pontual spring is considered

        Input :
            itself (a contact_gw_tempo)
        Output :
            Nothing, but attributes are updated
    """
    #conditions "if" are defined and same for each wall nature
    if self.nature == 'gwy_min':
        #unlinear stiffness
        nwg = np.array([0,1])
        self.nwg = nwg
        Fwg_n = self.k*self.overlap**(3/2)
        Fwg = Fwg_n*nwg
        self.Fwg_n = Fwg_n
        self.g.add_F(Fwg, self.g.l_border[self.g.l_border_y.index(min(self.g.l_border_y))])
        #damping
        gamma = -math.log(self.coeff_restitution)/math.sqrt(math.pi**2+math.log(self.coeff_restitution)**2)
        mass_eq = self.g.mass
        eta = 2 * gamma * math.sqrt(mass_eq*self.k)
        Fwg_damp_n = -np.dot(self.g.v,nwg)*eta
        Fwg_damp = Fwg_damp_n*nwg
        self.Fwg_damp_n = Fwg_damp_n
        self.g.add_F(Fwg_damp, self.g.l_border[self.g.l_border_y.index(min(self.g.l_border_y))])

    elif self.nature == 'gwy_max':
        #unlinear stiffness
        nwg = np.array([0,-1])
        self.nwg = nwg
        Fwg_n = self.k*self.overlap**(3/2)
        Fwg = Fwg_n*nwg
        self.Fwg_n = Fwg_n
        self.g.add_F(Fwg, self.g.l_border[self.g.l_border_y.index(max(self.g.l_border_y))])
        #damping
        Fwg_damp_n = 0
        self.Fwg_damp_n = Fwg_damp_n

    elif self.nature == 'gwx_min':
        #unlinear stiffness
        nwg = np.array([1,0])
        self.nwg = nwg
        Fwg_n = self.k*self.overlap**(3/2)
        Fwg = Fwg_n*nwg
        self.Fwg_n = Fwg_n
        self.g.add_F(Fwg, self.g.l_border[self.g.l_border_x.index(min(self.g.l_border_x))])
        #damping
        gamma = -math.log(self.coeff_restitution)/math.sqrt(math.pi**2+math.log(self.coeff_restitution)**2)
        mass_eq = self.g.mass
        eta = 2 * gamma * math.sqrt(mass_eq*self.k)
        Fwg_damp_n = -np.dot(self.g.v,nwg)*eta
        Fwg_damp = Fwg_damp_n*nwg
        self.Fwg_damp_n = Fwg_damp_n
        self.g.add_F(Fwg_damp, self.g.l_border[self.g.l_border_x.index(min(self.g.l_border_x))])

    elif self.nature == 'gwx_max':
        #unlinear stiffness
        nwg = np.array([-1,0])
        self.nwg = nwg
        Fwg_n = self.k*self.overlap**(3/2)
        Fwg = Fwg_n*nwg
        self.Fwg_n = Fwg_n
        self.g.add_F(Fwg, self.g.l_border[self.g.l_border_x.index(max(self.g.l_border_x))])
        #damping
        gamma = -math.log(self.coeff_restitution)/math.sqrt(math.pi**2+math.log(self.coeff_restitution)**2)
        mass_eq = self.g.mass
        eta = 2 * gamma * math.sqrt(mass_eq*self.k)
        Fwg_damp_n = -np.dot(self.g.v,nwg)*eta
        Fwg_damp = Fwg_damp_n*nwg
        self.Fwg_damp_n = Fwg_damp_n
        self.g.add_F(Fwg_damp, self.g.l_border[self.g.l_border_x.index(max(self.g.l_border_x))])

#-------------------------------------------------------------------------------

  def tangential(self, dt_DEM):
   """
   Compute the tangential reaction of a contact grain-wall.

   Here a pontual spring is considered.

        Input :
            itself (a contact_gw_tempo)
            a time step (a float)
        Output :
            Nothing, but attributes are updated
   """
   #conditions "if" are defined and same for each wall nature
   if self.nature == 'gwy_min':
       #unlinear stiffness
       twg = np.array([-1, 0])
       self.twg = twg
       r = np.linalg.norm(self.g.l_border[:-1][self.g.l_border_y.index(min(self.g.l_border_y))] - self.g.center) - self.overlap
       Delta_Us = (np.dot(self.g.v,self.twg) - r*self.g.w) * dt_DEM
       self.overlap_tangential = self.overlap_tangential + Delta_Us
       self.ft = self.ft - self.kt*Delta_Us
       if abs(self.ft) > abs(self.mu*self.Fwg_n) :
           self.ft = self.mu * abs(self.Fwg_n) * np.sign(self.ft)
       Fwg = self.ft*twg
       self.g.add_F(Fwg, self.g.l_border[self.g.l_border_y.index(min(self.g.l_border_y))])

   elif self.nature == 'gwy_max':
       #unlinear stiffness
       twg = np.array([1, 0])
       self.twg = twg
       r = np.linalg.norm(self.g.l_border[:-1][self.g.l_border_y.index(max(self.g.l_border_y))] - self.g.center) - self.overlap
       Delta_Us = (np.dot(self.g.v,self.twg) - r*self.g.w) * dt_DEM
       self.overlap_tangential = self.overlap_tangential + Delta_Us
       self.ft = self.ft - self.kt*Delta_Us
       if abs(self.ft) > abs(self.mu*self.Fwg_n) :
           self.ft = self.mu * abs(self.Fwg_n) * np.sign(self.ft)
       Fwg = self.ft*twg
       self.g.add_F(Fwg, self.g.l_border[self.g.l_border_y.index(max(self.g.l_border_y))])

   elif self.nature == 'gwx_min':
       #unlinear stiffness
       twg = np.array([0, 1])
       self.twg = twg
       r = np.linalg.norm(self.g.l_border[:-1][self.g.l_border_x.index(min(self.g.l_border_x))] - self.g.center) - self.overlap
       Delta_Us = (np.dot(self.g.v,self.twg) - r*self.g.w) * dt_DEM
       self.overlap_tangential = self.overlap_tangential + Delta_Us
       self.ft = self.ft - self.kt*Delta_Us
       if abs(self.ft) > abs(self.mu*self.Fwg_n) :
           self.ft = self.mu * abs(self.Fwg_n) * np.sign(self.ft)
       Fwg = self.ft*twg
       self.g.add_F(Fwg, self.g.l_border[self.g.l_border_x.index(min(self.g.l_border_x))])

   elif self.nature == 'gwx_max':
       #linear stiffness
       twg = np.array([0, -1])
       self.twg = twg
       r = np.linalg.norm(self.g.l_border[:-1][self.g.l_border_x.index(max(self.g.l_border_x))] - self.g.center) - self.overlap
       Delta_Us = (np.dot(self.g.v,self.twg) - r*self.g.w) * dt_DEM
       self.overlap_tangential = self.overlap_tangential + Delta_Us
       self.ft = self.ft - self.kt*Delta_Us
       if abs(self.ft) > abs(self.mu*self.Fwg_n) :
           self.ft = self.mu * abs(self.Fwg_n) * np.sign(self.ft)
       Fwg = self.ft*twg
       self.g.add_F(Fwg, self.g.l_border[self.g.l_border_x.index(max(self.g.l_border_x))])

#-------------------------------------------------------------------------------
#Function Definition
#-------------------------------------------------------------------------------

def LG_tempo(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report):
    """
    Create an initial condition

        Input :
            an algorithm dictionnary (a dict)
            a geometry dictionnary (a dict)
            an initial condition dictionnary (a dict)
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
            a sollicitations dictionnary (a dict)
            a simulation report (a report)
        Output :
            Nothing, but dictionnaries are updated
    """
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    #load data needed
    n_generation = dict_ic['n_generation']
    if n_generation != 2:
        simulation_report.write('n_generation must be equal to 2 !')
        raise ValueError('n_generation must be equal to 2 !')
    factor = dict_ic['factor_ymax_box']
    N_grain_disk = dict_geometry['N_grain_disk']
    N_grain_square = dict_geometry['N_grain_square']
    N_grain = dict_geometry['N_grain_disk'] + dict_geometry['N_grain_square']
    L_radius = dict_geometry['L_R']
    L_percentage_radius = dict_geometry['L_percentage_R']
    L_dimension = dict_geometry['L_Dimension']
    L_percentage_dimension = dict_geometry['L_percentage_Dimension']
    x_min = dict_sample['x_box_min']
    x_max = dict_sample['x_box_max']
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

    #define the y_max for the grains generation
    radius_mean = 0
    for i in range(len(L_radius)):
        radius_mean = radius_mean + L_radius[i]*L_percentage_radius[i]
    dimension_mean = 0
    for i in range(len(L_dimension)):
        dimension_mean = dimension_mean + L_dimension[i]*L_percentage_dimension[i]
    dy_creation = N_grain_disk/n_generation*factor*(2*radius_mean)**2/(x_max-x_min) + N_grain_square/n_generation*factor*(dimension_mean)**2/(x_max-x_min)

    #plan the grains generation
    L_n_grain_radius_try_one = []
    L_n_grain_radius = []
    L_n_grain_radius_done = []
    for percentage in L_percentage_radius:
        L_n_grain_radius_try_one.append(int(N_grain_disk*percentage/n_generation))
        L_n_grain_radius.append(int(N_grain_disk*percentage))
        L_n_grain_radius_done.append(0)
    L_n_grain_dimension_try_one = []
    L_n_grain_dimension = []
    L_n_grain_dimension_done = []
    for percentage in L_percentage_dimension:
        L_n_grain_dimension_try_one.append(int(N_grain_square*percentage/n_generation))
        L_n_grain_dimension.append(int(N_grain_square*percentage))
        L_n_grain_dimension_done.append(0)

    #Creation of grains
    #grains generation is decomposed in several steps (creation of grain then settlement)
 
    #---------------------------------------------------------------------------

    print('First generation of grains') 

    #add elements in dicts
    dict_ic['L_g_tempo'] = []
    dict_ic['L_L_g_tempo'] = []
    dict_ic['i_DEM_IC'] = 0
    dict_ic['L_n_grain_radius_try_one'] = L_n_grain_radius_try_one
    dict_ic['L_n_grain_radius'] = L_n_grain_radius
    dict_ic['L_n_grain_radius_done'] = L_n_grain_radius_done
    dict_ic['L_n_grain_dimension_try_one'] = L_n_grain_dimension_try_one
    dict_ic['L_n_grain_dimension'] = L_n_grain_dimension
    dict_ic['L_n_grain_dimension_done'] = L_n_grain_dimension_done
    dict_ic['last_id'] = 0
    dict_sample['y_box_min_ic'] = dict_sample['y_box_min']
    dict_sample['dy_creation'] = dy_creation

    Create_grains(dict_ic, dict_geometry, dict_sample, dict_material, 1, simulation_report)

    #DEM to find the steady-state configuration after loading
    #find the maximum y (center+radius)
    y_max = dict_sample['y_box_min_ic']
    for grain in dict_ic['L_g_tempo']:
        if max(grain.l_border_y) > y_max:
            y_max = max(grain.l_border_y)

    #add element in dict
    dict_sample['y_box_max'] = y_max

    DEM_loading(dict_ic, dict_material, dict_sample, dict_sollicitations, False, simulation_report)

    #---------------------------------------------------------------------------

    print('Second generation of grains')

    #Update dict
    dict_sample['y_box_min_ic'] = dict_sample['y_box_max']
    dict_ic['L_g_tempo'] = []

    Create_grains(dict_ic, dict_geometry, dict_sample, dict_material, 2, simulation_report)

    #DEM to find the steady-state configuration after loading
    #find the maximum y (center+radius)
    y_max = dict_sample['y_box_min_ic']
    for grain in dict_ic['L_g_tempo']:
        if max(grain.l_border_y) > y_max:
            y_max = max(grain.l_border_y)

    #add element in dict
    dict_sample['y_box_max'] = y_max

    DEM_loading(dict_ic, dict_material, dict_sample, dict_sollicitations, False, simulation_report)

    #update element in dict
    dict_sample['y_box_min_ic'] = dict_sample['y_box_max']

    #---------------------------------------------------------------------------

    print('Combine generations of grains')

    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    #load data needed
    L_L_g_tempo = dict_ic['L_L_g_tempo']
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

    L_g = []
    for L_g_tempo in L_L_g_tempo:
        for g_tempo in L_g_tempo:
            L_g.append(g_tempo)

    #update element in dict
    dict_ic['L_g_tempo'] = L_g

    DEM_loading(dict_ic, dict_material, dict_sample, dict_sollicitations, True, simulation_report)

    #update element in dict
    dict_sample['y_box_max'] = dict_ic['y_box_max']

    simulation_report.write_and_print(str(len(dict_ic['L_g_tempo']))+' / '+str(N_grain)+' grains have been created\n','\n'+str(len(dict_ic['L_g_tempo']))+' / '+str(N_grain)+' grains have been created\n')

#-------------------------------------------------------------------------------

def DEM_loading(dict_ic, dict_material, dict_sample, dict_sollicitations, multi_generation, simulation_report):
    """
    Loading the granular system.

        Input :
            an initial condition dictionnary (a dict)
            a material dictionnary (a dict)
            a smaple dictionnary (a dict)
            a sollicitations dictionnary (a dict)
            a Boolean to determine if multiple generations are considered (a Boolean)
            a simultion report (a report)
        Output :
            Nothing, but initial condition dictionnary is updated
    """
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    #load data needed
    dt_DEM = dict_ic['dt_DEM_IC']
    i_DEM_stop = dict_ic['i_DEM_stop_IC']
    i_DEM = dict_ic['i_DEM_IC']
    Ecin_ratio_IC = dict_ic['Ecin_ratio_IC']
    i_print_plot_IC = dict_ic['i_print_plot_IC']
    factor_neighborhood_IC = dict_ic['factor_neighborhood_IC']
    if multi_generation :
        i_update_neighbouroods = dict_ic['i_update_neighborhoods_com']
        y_min = dict_sample['y_box_min']
    else :
        i_update_neighbouroods = dict_ic['i_update_neighborhoods_gen']
        y_min = dict_sample['y_box_min_ic']
    x_min = dict_sample['x_box_min']
    x_max = dict_sample['x_box_max']
    y_max = dict_sample['y_box_max']
    Forcev_target = dict_sollicitations['Vertical_Confinement_Force']
    gravity = dict_sollicitations['gravity']
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

    i_DEM_0 = i_DEM
    DEM_loop_statut = True

    #Initialisation
    dict_ic['L_contact'] = []
    dict_ic['L_contact_ij'] = []
    dict_ic['L_contact_gw'] = []
    dict_ic['L_contact_gw_ij'] = []
    dict_ic['id_contact'] = 0

    #trackers and stop conditions
    Force_tracker = []
    Force_stop = 0
    Ecin_tracker = []
    Ecin_stop = 0
    Ymax_tracker = []
    k0_tracker = []
    Fv_tracker = []
    for grain in dict_ic['L_g_tempo']:
        Force_stop = Force_stop + 0.5*grain.mass*gravity
        Ecin_stop = Ecin_stop + 0.5*grain.mass*(Ecin_ratio_IC*grain.r_max/dt_DEM)**2
        grain.v = np.array([0, 0])

    while DEM_loop_statut :

        #calm down particle
        if (i_DEM - i_DEM_0) == int(i_DEM_stop*0.4) or (i_DEM - i_DEM_0) == int(i_DEM_stop*0.5) or (i_DEM - i_DEM_0) == int(i_DEM_stop*0.6) : 
            for grain in dict_ic['L_g_tempo'] : 
                grain.v = np.array([0, 0])
        
        i_DEM = i_DEM + 1

        #compute neighborhoods 
        if (i_DEM - i_DEM_0) < i_DEM_stop*0.1 :
            if (i_DEM-i_DEM_0-1) % int(i_update_neighbouroods/2) == 0:
                Update_Neighbouroods(dict_ic)
            if (i_DEM-i_DEM_0-1) % int(i_update_neighbouroods/4) == 0:
                wall_neighborhood = Update_wall_Neighborhoods(dict_ic['L_g_tempo'],factor_neighborhood_IC,x_min,x_max,y_min,y_max)
        else :
            if (i_DEM-i_DEM_0-1) % int(i_update_neighbouroods) == 0:
                Update_Neighbouroods(dict_ic)
            if (i_DEM-i_DEM_0-1) % int(i_update_neighbouroods) == 0:
                wall_neighborhood = Update_wall_Neighborhoods(dict_ic['L_g_tempo'],factor_neighborhood_IC,x_min,x_max,y_min,y_max)
        
        #Contact detection
        Grains_Polyhedral_contact_Neighbouroods(dict_ic,dict_material)
        Grains_Polyhedral_Wall_contact_Neighborhood(wall_neighborhood,x_min,x_max,y_min,y_max, dict_ic, dict_material)

        #Sollicitation computation
        for grain in dict_ic['L_g_tempo']:
             grain.init_F_control(gravity)
        for contact in  dict_ic['L_contact']+dict_ic['L_contact_gw']:
            contact.normal()
            contact.tangential(dt_DEM)

        #Move grains
        for grain in dict_ic['L_g_tempo']:
            grain.euler_semi_implicite(dt_DEM,25*Ecin_ratio_IC)

        #check if some grains are outside of the study box
        L_ig_to_delete = []
        for id_grain in range(len(dict_ic['L_g_tempo'])):
            if dict_ic['L_g_tempo'][id_grain].center[0] < x_min :
                L_ig_to_delete.append(id_grain)
            elif dict_ic['L_g_tempo'][id_grain].center[0] > x_max :
                L_ig_to_delete.append(id_grain)
            elif dict_ic['L_g_tempo'][id_grain].center[1] < y_min :
                L_ig_to_delete.append(id_grain)
            elif dict_ic['L_g_tempo'][id_grain].center[1] > y_max :
                L_ig_to_delete.append(id_grain)
        L_ig_to_delete.reverse()
        for id_grain in L_ig_to_delete:
            simulation_report.write_and_print('Grain '+str(dict_ic['L_g_tempo'][id_grain].id)+' has been deleted because it is out of the box\n','Grain '+str(dict_ic['L_g_tempo'][id_grain].id)+' has been deleted because it is out of the box')
            dict_ic['L_g_tempo'].pop(id_grain)

        #Control the y_max to have the pressure target
        y_max, Fv = Control_y_max_NR(y_max,Forcev_target,dict_ic['L_contact_gw'],dict_ic['L_g_tempo'])
        Fv_tracker.append(Fv)

        #Tracker
        F = F_total(dict_ic['L_g_tempo'])
        Ecin = E_cin_total(dict_ic['L_g_tempo'])
        Force_tracker.append(F)
        Ecin_tracker.append(Ecin)
        Ymax_tracker.append(y_max)
        F_xmin = 0
        F_xmax = 0
        for contact in dict_ic['L_contact_gw']:
            if contact.nature == 'gwx_min' :
                F_xmin = F_xmin + contact.Fwg_n
            if contact.nature == 'gwx_max' :
                F_xmax = F_xmax + contact.Fwg_n
        if Fv != 0: 
            k0 = np.mean([F_xmin/(y_max- y_min)*(x_max- x_min)/Fv,
                          F_xmax/(y_max- y_min)*(x_max- x_min)/Fv])
        else :
            k0 = 0
        k0_tracker.append(k0)

        if i_DEM % i_print_plot_IC ==0:
            if gravity > 0 :
                print('i_DEM',i_DEM,'and Ecin',int(100*Ecin/Ecin_stop),'% and Force',int(100*F/Force_stop),'% and Confinement',int(100*Fv/Forcev_target),'%')
            else :
                print('i_DEM',i_DEM,'and Ecin',int(100*Ecin/Ecin_stop),'% and Confinement',int(100*Fv/Forcev_target),'%')
            if dict_ic['Debug_DEM']:
                Plot_Config_Loaded(dict_ic['L_g_tempo'],x_min,x_max,y_min,y_max,i_DEM)

        #Check stop conditions for DEM
        if i_DEM >= i_DEM_stop + i_DEM_0:
             DEM_loop_statut = False
        if gravity > 0:
            if Ecin < Ecin_stop and F < Force_stop and (0.95*Forcev_target<Fv and Fv<1.05*Forcev_target):
                  DEM_loop_statut = False
        else:
            if Ecin < Ecin_stop and i_DEM >= i_DEM_stop*0.7 + i_DEM_0 and (0.95*Forcev_target<Fv and Fv<1.05*Forcev_target):
                DEM_loop_statut = False
        if dict_ic['L_g_tempo'] == []:
            DEM_loop_statut = False

    #Update dict
    dict_ic['L_L_g_tempo'].append(dict_ic['L_g_tempo'].copy())
    dict_ic['y_box_max'] = y_max
    dict_ic['i_DEM_IC'] = i_DEM
    dict_ic['Ecin_tracker'] = Ecin_tracker
    dict_ic['Ymax_tracker'] = Ymax_tracker
    dict_ic['k0_tracker'] = k0_tracker
    dict_ic['Fv_tracker'] = Fv_tracker

#-------------------------------------------------------------------------------

def Create_grains(dict_ic, dict_geometry, dict_sample, dict_material, id_generation, simulation_report):
    """
    Generate the grains.

    A position is tried. This position must not create overlap with already creared temporary grain. If there is no overlap, a temporary grai nis created.

        Input :
            an initial condition dictionnary (a dict)
            a geometry dictionnary (a dict)
            a sample dictionnary (a dict)
            a material dictionnary (a dict)
            a generation id (a int)
            a simulation report (a report)
        Output :
            Nothing, but initial configuration dictionnary is updated
    """
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    #load data needed
    L_g_tempo = dict_ic['L_g_tempo']
    N_test_max = dict_ic['N_test_max']
    if id_generation == 1:
        L_n_grain_radius = dict_ic['L_n_grain_radius_try_one']
        L_n_grain_dimension = dict_ic['L_n_grain_dimension_try_one']
    else :
        L_n_grain_radius = dict_ic['L_n_grain_radius']
        L_n_grain_dimension = dict_ic['L_n_grain_dimension']
    L_n_grain_radius_done = dict_ic['L_n_grain_radius_done']
    L_radius = dict_geometry['L_R']
    L_n_grain_dimension_done = dict_ic['L_n_grain_dimension_done']
    L_dimension = dict_geometry['L_Dimension']
    x_min = dict_sample['x_box_min']
    x_max = dict_sample['x_box_max']
    y_min = dict_sample['y_box_min_ic']
    dy_creation = dict_sample['dy_creation']
    y_max = y_min + dy_creation
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

    #Parameters for the method
    n_not_created = 0

    for i in range(len(L_radius)):
        radius = L_radius[i]
        n_grain = L_n_grain_radius[i]
        n_grain_done = L_n_grain_radius_done[i]
        last_id_grain_created = dict_ic['last_id']
        for id_grain in range(last_id_grain_created, last_id_grain_created + n_grain - n_grain_done):
            i_test = 0
            grain_created = False
            while (not grain_created) and i_test < N_test_max:
                i_test = i_test + 1
                center = np.array([random.uniform(x_min+1.1*radius,x_max-1.1*radius),random.uniform(y_min+1.1*radius,y_max)])
                g_tempo = Grain_Tempo(id_grain-n_not_created,center,radius,dict_geometry['grain_discretisation'],dict_material,'Disk')
                grain_created = True
                for grain in L_g_tempo:
                    if Grains_Polyhedral_contact_f(g_tempo,grain):
                        grain_created = False
            if i_test == N_test_max and not grain_created:
                n_not_created = n_not_created + 1
                simulation_report.write_and_print('Grain '+str(id_grain)+' has not been created after '+str(i_test)+' tries\n','Grain '+str(id_grain)+' has not been created after '+str(i_test)+' tries')
            else :
                L_g_tempo.append(g_tempo)
                L_n_grain_radius_done[i] = L_n_grain_radius_done[i] + 1
                dict_ic['last_id'] = dict_ic['last_id'] + 1
    for i in range(len(L_dimension)):
        dimension = L_dimension[i]
        n_grain = L_n_grain_dimension[i]
        n_grain_done = L_n_grain_dimension_done[i]
        last_id_grain_created = dict_ic['last_id']
        for id_grain in range(last_id_grain_created, last_id_grain_created + n_grain - n_grain_done):
            i_test = 0
            grain_created = False
            while (not grain_created) and i_test < N_test_max:
                i_test = i_test + 1
                center = np.array([random.uniform(x_min+1.1*dimension/2,x_max-1.1*dimension/2),random.uniform(y_min+1.1*dimension/2,y_max)])
                g_tempo = Grain_Tempo(id_grain-n_not_created,center,dimension,dict_geometry['grain_discretisation'],dict_material,'Square')
                grain_created = True
                for grain in L_g_tempo:
                    if Grains_Polyhedral_contact_f(g_tempo,grain):
                        grain_created = False
            if i_test == N_test_max and not grain_created:
                n_not_created = n_not_created + 1
                simulation_report.write_and_print('Grain '+str(id_grain)+' has not been created after '+str(i_test)+' tries\n','Grain '+str(id_grain)+' has not been created after '+str(i_test)+' tries')
            else :
                L_g_tempo.append(g_tempo)
                L_n_grain_dimension_done[i] = L_n_grain_dimension_done[i] + 1
                dict_ic['last_id'] = dict_ic['last_id'] + 1

    #Update dict
    dict_ic['L_g_tempo'] = L_g_tempo
    dict_ic['L_n_grain_dimension_done'] = L_n_grain_dimension_done
    dict_ic['L_n_grain_radius_done'] = L_n_grain_radius_done

#-------------------------------------------------------------------------------

def Grains_Polyhedral_contact_f(g1,g2):
  """
  Detect the contact grain-grain.

    Input :
        two temporary grains (two grain_tempos)
    Output :
        a Boolean, True if there is contaxct between the twwo grains (a Boolean)
  """
  if np.linalg.norm(g1.center-g2.center) < 1.5*(g1.r_max+g2.r_max):

      #looking for the nearest nodes
      d_virtual = max(g1.r_max,g2.r_max)
      ij_min = [0,0]
      d_ij_min = 100*d_virtual #Large
      for i in range(len(g1.l_border[:-1])):
        for j in range(len(g2.l_border[:-1])):
            d_ij = np.linalg.norm(g2.l_border[:-1][j]-g1.l_border[:-1][i]+d_virtual*(g2.center-g1.center)/np.linalg.norm(g2.center-g1.center))
            if d_ij < d_ij_min :
                d_ij_min = d_ij
                ij_min = [i,j]

      d_ij_min = np.dot(g2.l_border[:-1][ij_min[1]]-g1.l_border[:-1][ij_min[0]],-(g2.center-g1.center)/np.linalg.norm(g2.center-g1.center))
      return d_ij_min > 0

  else:
    return False

#-------------------------------------------------------------------------------

def E_cin_total(L_g):
    """
    Compute total kinetic energy.

        Input :
            a list of temporary grains (a list)
        Output :
            the total kinetic energy (a float)
    """
    Ecin = 0
    for grain in L_g:
        Ecin = Ecin + 1/2*grain.mass*np.dot(grain.v,grain.v)
    return Ecin

#-------------------------------------------------------------------------------

def F_total(L_g):
    """
    Compute total force applied on grains in the sample.

        Input :
            a list of temporary grains (a list)
        Output :
            the total force applied (a float)
    """
    F = 0
    for grain in L_g:
        F = F + np.linalg.norm([grain.fx, grain.fy])
    return F

#-------------------------------------------------------------------------------

def Control_y_max_NR(y_max,Force_target,L_contact_gw,L_g):
    """
    Control the upper wall to apply force.

    A Newton-Raphson method is applied to verify the confinement.
        Input :
            a coordinate of the upper wall (a float)
            a confinement value (a float)
            a list of contact grain - wall (a list)
            a list of temporary grain (a list)
        Output :
            the coordinate of the upper wall (a float)
            a force applied on the upper wall before control (a float)
    """
    F = 0
    overlap_L = []
    k_L = []
    for contact in L_contact_gw:
        if contact.nature == 'gwy_max':
            F = F + contact.Fwg_n
            overlap_L.append(contact.overlap)
            k_L.append(contact.k)
            #compute force applied, save contact overlap and spring

    if overlap_L != []:
        i_NR = 0
        dy = 0
        ite_criteria = True
        #control the upper wall
        if -0.01*Force_target<error_on_ymax_f(dy,overlap_L,k_L,Force_target) and error_on_ymax_f(dy,overlap_L,k_L,Force_target)<0.01*Force_target:
            ite_criteria = False
        while ite_criteria :
            i_NR = i_NR + 1
            dy = dy - error_on_ymax_f(dy,overlap_L,k_L,Force_target)/error_on_ymax_df(dy,overlap_L,k_L)
            if i_NR > 100: #Maximum try
                ite_criteria = False
            if -0.01*Force_target<error_on_ymax_f(dy,overlap_L,k_L,Force_target) and error_on_ymax_f(dy,overlap_L,k_L,Force_target)<0.01*Force_target:
                ite_criteria = False
        y_max = y_max + dy

    else :
        #if there is no contact with the upper wall, the wall is reset
        y_max = Reset_y_max(L_g,Force_target)

    for contact in L_contact_gw:
        if contact.nature == 'gwy_max':
            #reactualisation
            contact.limit = y_max

    return y_max, F

#-------------------------------------------------------------------------------

def error_on_ymax_f(dy,overlap_L,k_L,Force_target) :
    """
    Compute the function f to control the upper wall. It is the difference between the force applied and the target value.

        Input :
            an increment of the upper wall position (a float)
            a list of overlap for contact between temporary grain and upper wall (a list)
            a list of spring for contact between temporary grain and upper wall (a list)
            a confinement force (a float)
        Output :
            the difference between the force applied and the confinement (a float)
    """
    f = Force_target
    for i in range(len(overlap_L)):
        f = f - k_L[i]*(max(overlap_L[i]-dy,0))**(3/2)
    return f

#-------------------------------------------------------------------------------

def error_on_ymax_df(dy,overlap_L,k_L) :
    """
    Compute the derivative function df to control the upper wall (error_on_ymax_f()).

        Input :
            an increment of the upper wall position (a float)
            a list of overlap for contact between temporary grain and upper wall (a list)
            a list of spring for contact between temporary grain and upper wall (a list)
        Output :
            the derivative of error_on_ymax_f() (a float)
    """
    df = 0
    for i in range(len(overlap_L)):
        df = df + 3/2*k_L[i]*(max(overlap_L[i]-dy,0))**(1/2)
    return df

#-------------------------------------------------------------------------------

def Reset_y_max(L_g,Force):
    """
    The upper wall is located as a single contact verify the target value.

        Input :
            the list of temporary grains (a list)
            the confinement force (a float)
        Output :
            the upper wall position (a float)
    """
    y_max = None
    id_grain_max = None
    for id_grain in range(len(L_g)):
        grain = L_g[id_grain]
        y_max_grain = max(grain.l_border_y)

        if y_max != None and y_max_grain > y_max:
            y_max = y_max_grain
            id_grain_max = id_grain
        elif y_max == None:
            y_max = y_max_grain
            id_grain_max = id_grain

    factor = 5
    k = factor*4/3*L_g[id_grain_max].y/(1-L_g[id_grain_max].nu*L_g[id_grain_max].nu)*math.sqrt(L_g[id_grain_max].r_max)
    y_max = y_max - (Force/k)**(2/3)

    return y_max

#-------------------------------------------------------------------------------

def Grains_Polyhedral_contact_Neighborhoods_f(g1,g2):
  """
  Detect contact grain - grain.

    Input :
        two temporary grains (two grain_tempo)
    Output :
        a Boolean, True if there is a contact (a Boolean)
  """
  #looking for the nearest nodes
  d_virtual = max(g1.r_max,g2.r_max)
  ij_min = [0,0]
  d_ij_min = 100*d_virtual #Large
  for i in range(len(g1.l_border[:-1])):
    for j in range(len(g2.l_border[:-1])):
        d_ij = np.linalg.norm(g2.l_border[:-1][j]-g1.l_border[:-1][i]+d_virtual*(g2.center-g1.center)/np.linalg.norm(g2.center-g1.center))
        if d_ij < d_ij_min :
            d_ij_min = d_ij
            ij_min = [i,j]

  d_ij_min = np.dot(g2.l_border[:-1][ij_min[1]]-g1.l_border[:-1][ij_min[0]],-(g2.center-g1.center)/np.linalg.norm(g2.center-g1.center))
  return d_ij_min > 0

#-------------------------------------------------------------------------------

def Update_Neighbouroods(dict_ic):
    """
    Determine a neighborhood for each grain.

    This function is called every x time step. The contact is determined by Grains_Polyhedral_contact_Neighbouroods().
    Notice that if there is a potential contact between grain_i and grain_j, grain_i is not in the neighborhood of grain_j.
    Whereas grain_j is in the neighbourood of grain_i. With i_grain < j_grain.

        Input :
            an initial condition dictionnary (a dict)
        Output :
            Nothing, but the neighborhood of the temporary grains is updated
    """
    for i_grain in range(len(dict_ic['L_g_tempo'])-1) :
        neighbourood = []
        for j_grain in range(i_grain+1,len(dict_ic['L_g_tempo'])):
            if np.linalg.norm(dict_ic['L_g_tempo'][i_grain].center-dict_ic['L_g_tempo'][j_grain].center) < dict_ic['factor_neighborhood_IC']*(dict_ic['L_g_tempo'][i_grain].r_max+dict_ic['L_g_tempo'][j_grain].r_max):
                neighbourood.append(dict_ic['L_g_tempo'][j_grain])
        dict_ic['L_g_tempo'][i_grain].neighbourood = neighbourood

#-------------------------------------------------------------------------------

def Grains_Polyhedral_contact_Neighbouroods(dict_ic,dict_material):
    """
    Detect contact between a grain and grains from its neighbourood.

    The neighbourood is updated with Update_Neighbouroods().

        Input :
            an initial condition dictionnary (a dict)
            a material dictionnary (a dict)
        Output :
            Nothing, but the initial condition dictionnary is updated with grain - grain contacts
    """
    for i_grain in range(len(dict_ic['L_g_tempo'])-1) :
        grain_i = dict_ic['L_g_tempo'][i_grain]
        for neighbour in dict_ic['L_g_tempo'][i_grain].neighbourood:
            j_grain = neighbour.id
            grain_j = neighbour
            if Grains_Polyhedral_contact_Neighborhoods_f(grain_i,grain_j):
                if (i_grain,j_grain) not in dict_ic['L_contact_ij']:  #contact not detected previously
                   #creation of contact
                   dict_ic['L_contact_ij'].append((i_grain,j_grain))
                   dict_ic['L_contact'].append(Contact_Tempo(dict_ic['id_contact'], grain_i, grain_j, dict_material))
                   dict_ic['id_contact'] = dict_ic['id_contact'] + 1

            else :
                if (i_grain,j_grain) in dict_ic['L_contact_ij'] : #contact detected previously is not anymore
                       dict_ic['L_contact'].pop(dict_ic['L_contact_ij'].index((i_grain,j_grain)))
                       dict_ic['L_contact_ij'].remove((i_grain,j_grain))

#-------------------------------------------------------------------------------

def Update_wall_Neighborhoods(L_g_tempo,factor_neighborhood_IC,x_min,x_max,y_min,y_max):
    """
    Determine a neighborhood for wall.

    This function is called every x time step. The grain - wall contact is determined by Grains_Polyhedral_Wall_contact_Neighborhood().
    A factor determines the size of the neighborhood window.

        Input :
            a list of temporary grains (a list)
            a factor to determine the neighborhood window (a float)
            the coordinates of the left, right, lower, upper walls (four floats)
        Output :
            a list of temporary grains in the neighborhood of the walls (a list)
    """
    wall_neighborhood = []
    for grain in L_g_tempo:

        p_x_min = min(grain.l_border_x)
        p_x_max = max(grain.l_border_x)
        p_y_min = min(grain.l_border_y)
        p_y_max = max(grain.l_border_y)

        #grain-wall x_min
        if abs(p_x_min-x_min) < factor_neighborhood_IC*grain.r_max :
            wall_neighborhood.append(grain)
        #grain-wall x_max
        if abs(p_x_max-x_max) < factor_neighborhood_IC*grain.r_max :
            wall_neighborhood.append(grain)
        #grain-wall y_min
        if abs(p_y_min-y_min) < factor_neighborhood_IC*grain.r_max :
            wall_neighborhood.append(grain)
        #grain-wall y_max
        if abs(p_y_max-y_max) < factor_neighborhood_IC*grain.r_max :
            wall_neighborhood.append(grain)

    return wall_neighborhood

#-------------------------------------------------------------------------------

def Grains_Polyhedral_Wall_contact_Neighborhood(wall_neighborhood,x_box_min,x_box_max,y_box_min,y_box_max, dict_ic, dict_material):
  """
  Detect contact grain in the neighborhood of the wall and the wall.

  The neighbourood is updated with Update_wall_Neighborhoods(). An iteration over the grains in the wall neighborhood is done. A comparison is done with the coordinates of the wall to determine if there is a contact.

        Input :
            a walls neighborhood (a list)
            the coordinates of the walls (four floats)
            an initial condition dictionnary (a dict)
            a material dictionnary (a dict)
        Output :
            Nothing, but the initial condition dictionnary is updated with the contact grain - walls.
  """
  #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
  #load data needed
  L_ij_contact_gw = dict_ic['L_contact_gw_ij']
  L_contact_gw = dict_ic['L_contact_gw']
  id_contact = dict_ic['id_contact']
  #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

  for grain in wall_neighborhood:

      p_x_min = min(grain.l_border_x)
      p_x_max = max(grain.l_border_x)
      p_y_min = min(grain.l_border_y)
      p_y_max = max(grain.l_border_y)

      #grain-wall x_min
      if p_x_min < x_box_min and (grain.id,-1) not in L_ij_contact_gw:
          overlap = x_box_min - p_x_min
          L_contact_gw.append(Contact_gw_Tempo(id_contact, grain, dict_material, 'gwx_min', x_box_min, overlap))
          L_ij_contact_gw.append((grain.id,-1))
          id_contact = id_contact + 1
      elif p_x_min < x_box_min and (grain.id,-1) in L_ij_contact_gw:
          overlap = x_box_min - p_x_min
          L_contact_gw[L_ij_contact_gw.index((grain.id,-1))].update_overlap(overlap)
      elif p_x_min > x_box_min and (grain.id,-1) in L_ij_contact_gw:
          i_contact = L_ij_contact_gw.index((grain.id,-1))
          L_contact_gw.pop(i_contact)
          L_ij_contact_gw.pop(i_contact)
      #grain-wall x_max
      if p_x_max > x_box_max and (grain.id,-2) not in L_ij_contact_gw:
          overlap = p_x_max - x_box_max
          L_contact_gw.append(Contact_gw_Tempo(id_contact, grain, dict_material, 'gwx_max', x_box_max, overlap))
          L_ij_contact_gw.append((grain.id,-2))
          id_contact = id_contact + 1
      elif p_x_max > x_box_max and (grain.id,-2) in L_ij_contact_gw:
          overlap = p_x_max - x_box_max
          L_contact_gw[L_ij_contact_gw.index((grain.id,-2))].update_overlap(overlap)
      elif p_x_max < x_box_max and (grain.id,-2) in L_ij_contact_gw:
          i_contact = L_ij_contact_gw.index((grain.id,-2))
          L_contact_gw.pop(i_contact)
          L_ij_contact_gw.pop(i_contact)
      #grain-wall y_min
      if p_y_min < y_box_min and (grain.id,-3) not in L_ij_contact_gw:
          overlap = y_box_min - p_y_min
          L_contact_gw.append(Contact_gw_Tempo(id_contact, grain, dict_material, 'gwy_min', y_box_min, overlap))
          L_ij_contact_gw.append((grain.id,-3))
          id_contact = id_contact + 1
      elif p_y_min < y_box_min and (grain.id,-3) in L_ij_contact_gw:
          overlap = y_box_min - p_y_min
          L_contact_gw[L_ij_contact_gw.index((grain.id,-3))].update_overlap(overlap)
      elif p_y_min > y_box_min and (grain.id,-3) in L_ij_contact_gw:
          i_contact = L_ij_contact_gw.index((grain.id,-3))
          L_contact_gw.pop(i_contact)
          L_ij_contact_gw.pop(i_contact)
      #grain-wall y_max
      if p_y_max > y_box_max and (grain.id,-4) not in L_ij_contact_gw:
          overlap = p_y_max - y_box_max
          L_contact_gw.append(Contact_gw_Tempo(id_contact, grain, dict_material, 'gwy_max', y_box_max, overlap))
          L_ij_contact_gw.append((grain.id,-4))
          id_contact = id_contact + 1
      elif p_y_max > y_box_max and (grain.id,-4) in L_ij_contact_gw:
          overlap = p_y_max - y_box_max
          L_contact_gw[L_ij_contact_gw.index((grain.id,-4))].update_overlap(overlap)
      elif p_y_max < y_box_max and (grain.id,-4) in L_ij_contact_gw:
          i_contact = L_ij_contact_gw.index((grain.id,-4))
          L_contact_gw.pop(i_contact)
          L_ij_contact_gw.pop(i_contact)

      #Update dict
      dict_ic['L_contact_gw_ij'] = L_ij_contact_gw
      dict_ic['L_contact_gw'] = L_contact_gw
      dict_ic['id_contact'] = id_contact

#-------------------------------------------------------------------------------

def Plot_Config_Loaded(L_g,x_min,x_max,y_min,y_max,i):
    """
    Plot loaded configuration.

        Input :
            a list of temporary grain (a list)
            the coordinates of the walls (four floats)
            an iteration (a int)
        Output :
            Nothing, but a .png file is generated (a file)
    """
    plt.figure(1,figsize=(16,9))
    L_x = []
    L_y = []
    L_u = []
    L_v = []
    for grain in L_g:
        plt.plot(grain.l_border_x,grain.l_border_y,'k')
        plt.plot(grain.center[0],grain.center[1],'xk')
        L_x.append(grain.center[0])
        L_y.append(grain.center[1])
        L_u.append(grain.fx)
        L_v.append(grain.fy)
    plt.plot([x_min,x_min,x_max,x_max,x_min],[y_max,y_min,y_min,y_max,y_max],'k')
    plt.axis('equal')
    plt.savefig('Debug/Config_Loaded_'+str(i)+'.png')
    plt.close(1)

#-------------------------------------------------------------------------------

def Plot_Config_Loaded_End(L_g,x_min,x_max,y_min,y_max):
    """
    Plot loaded configuration at the end of the initial configuration.

        Input :
            a list of temporary grain (a list)
            the coordinates of the walls (four floats)
            an iteration (a int)
        Output :
            Nothing, but a .png file is generated (a file)
    """
    plt.figure(1,figsize=(16,9))
    L_x = []
    L_y = []
    L_u = []
    L_v = []
    for grain in L_g:
        plt.plot(grain.l_border_x,grain.l_border_y,'k')
        plt.plot(grain.center[0],grain.center[1],'xk')
        L_x.append(grain.center[0])
        L_y.append(grain.center[1])
        L_u.append(grain.fx)
        L_v.append(grain.fy)
    plt.plot([x_min,x_min,x_max,x_max,x_min],[y_max,y_min,y_min,y_max,y_max],'k')
    plt.axis('equal')
    plt.savefig('Debug/ConfigLoaded.png')
    plt.close(1)

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
    i_minus = find_value_in_list(g.l_theta_r.copy(), angle_minus)
    i_plus = i_minus + find_value_in_list(g.l_theta_r[i_minus:].copy() + g.l_theta_r[:i_minus].copy(), angle_plus)

    L_i_vertices = []
    for i in range(i_minus, i_plus+1):
        if i < len(g.l_theta_r):
            L_i_vertices.append(i)
        else :
            L_i_vertices.append(i-len(g.l_theta_r))
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

#-------------------------------------------------------------------------------