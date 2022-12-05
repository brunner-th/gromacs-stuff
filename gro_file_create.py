# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 16:40:17 2022

"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d



molecular_density = 33 # nm^-1

width_reservoir = 5
height_reservoir = 5 
lenght_reservoir = 1.5
lenght_channel = 5
radius_channel = 1

delta_wall = 0.3

m = 18.02 #g/mol
rho = 997 #kg/m^3
avogadro = 6.022*10**23 #1/mol



m_si = m/10**3
n_mol = rho/m_si # mol/m^3
v_mol = m_si/rho # m^3/mol
vol_per_molecule = v_mol/avogadro #m^3
vol_per_molecule_nm3 = vol_per_molecule*10**27 #nm^3
molecule_per_nm3 = 1/vol_per_molecule_nm3 #nm^-3
molecules_per_lenght = (molecule_per_nm3)**(1/3)




############## reservoir ######################################################


deltax = 1/molecules_per_lenght

k1 = 5 # adding const

x_vec = np.linspace(0,width_reservoir,int(width_reservoir*molecules_per_lenght)+k1)

y_vec = np.linspace(0,height_reservoir,int(height_reservoir*molecules_per_lenght))

z_vec = np.linspace(0,lenght_reservoir-delta_wall,int((lenght_reservoir-delta_wall)*molecules_per_lenght))

molecules_grid = len(x_vec)*len(y_vec)*len(z_vec)

molecule_per_nm3_grid = molecules_grid/(width_reservoir*height_reservoir*(lenght_reservoir-delta_wall))

print("reservoir density: %f, change k1" % (molecule_per_nm3_grid))

reserv1_vec = np.zeros((6,molecules_grid))

for idz, z in enumerate(z_vec):
    reserv1_vec[2,idz*len(x_vec)*len(y_vec):(idz+1)*len(x_vec)*len(y_vec)].fill(z)
    for idy, y in enumerate(y_vec):
        reserv1_vec[1,(idy)*len(x_vec)+idz*len(x_vec)*len(y_vec):(idy+1)*len(x_vec)+idz*len(x_vec)*len(y_vec)].fill(y)
        reserv1_vec[0,(idy)*len(x_vec)+idz*len(x_vec)*len(y_vec):(idy+1)*len(x_vec)+idz*len(x_vec)*len(y_vec)] = x_vec
        
        
reserv2_vec = np.copy(reserv1_vec)

reserv2_vec[2,:] = reserv2_vec[2,:]+lenght_channel+lenght_reservoir+delta_wall




################# wall ########################################################










################# channel #####################################################



k2 = 0 # adding const !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

x_vec = np.linspace(0,width_reservoir,int(width_reservoir*molecules_per_lenght)+k2)

y_vec = np.linspace(0,height_reservoir,int(height_reservoir*molecules_per_lenght))

z_vec = np.linspace(0,lenght_channel,int(lenght_channel*molecules_per_lenght))

molecules_grid = len(x_vec)*len(y_vec)*len(z_vec)

molecule_per_nm3_grid = molecules_grid/(width_reservoir*height_reservoir*lenght_channel)

print("channel density: %f, change k2" % (molecule_per_nm3_grid))


wall_vec = np.zeros((6,molecules_grid))

for idz, z in enumerate(z_vec):
    wall_vec[2,idz*len(x_vec)*len(y_vec):(idz+1)*len(x_vec)*len(y_vec)].fill(z+1.5)
    for idy, y in enumerate(y_vec):
        wall_vec[1,(idy)*len(x_vec)+idz*len(x_vec)*len(y_vec):(idy+1)*len(x_vec)+idz*len(x_vec)*len(y_vec)].fill(y)
        wall_vec[0,(idy)*len(x_vec)+idz*len(x_vec)*len(y_vec):(idy+1)*len(x_vec)+idz*len(x_vec)*len(y_vec)] = x_vec


vol_channel = lenght_channel*radius_channel**2*np.pi

num_of_molecules = molecule_per_nm3*vol_channel


channel_vec = np.zeros((6,molecules_grid))

delete_list = []
for idz, z in enumerate(wall_vec[2,:]):
    if (wall_vec[0,idz]-width_reservoir/2)**2 + (wall_vec[1,idz]-height_reservoir/2)**2 <= radius_channel**2:
        channel_vec[:,idz] = wall_vec[:,idz]
        wall_vec[:,idz] = np.zeros((6,))
        delete_list.append(idz)
        

    



# some ugly stuff

wall_vec = np.delete(wall_vec, delete_list, axis = 1)
idx_list = range(np.shape(channel_vec)[1])
z = list(set(idx_list) - set(delete_list))
channel_vec = np.delete(channel_vec, z, axis = 1)






################################################################################

fig = plt.figure()
ax = plt.axes(projection='3d')

    
ax.scatter(reserv1_vec[0,:], reserv1_vec[1,:], reserv1_vec[2,:], c='royalblue', marker='o', s = 2)
ax.scatter(reserv2_vec[0,:], reserv2_vec[1,:], reserv2_vec[2,:], c='firebrick', marker='o', s = 2)
#ax.scatter(wall_vec[0,:], wall_vec[1,:], wall_vec[2,:], c='gray', marker='o', s = 2)
ax.scatter(channel_vec[0,:], channel_vec[1,:], channel_vec[2,:], c='k', marker='o', s = 2)

ax.view_init(85, 50)
ax.view_init(20, 50)


