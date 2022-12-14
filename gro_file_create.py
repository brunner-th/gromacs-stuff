# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 16:40:17 2022

"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from math import sin, cos, pi

plt.rcParams['figure.figsize'] = [8, 6]
plt.rcParams['figure.dpi'] = 400
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

#################### input ####################################################

molecular_density = 33 # nm^-1

# dimensions:
width_reservoir = 5
height_reservoir = 5 
lenght_reservoir = 1.5
lenght_channel = 5
radius_channel = 0.5
delta_wall = 0.3

m = 18.02 #g/mol
rho = 997 #kg/m^3
avogadro = 6.022*10**23 #1/mol

################ basic stuff ##################################################

m_si = m/10**3
n_mol = rho/m_si # mol/m^3
v_mol = m_si/rho # m^3/mol
vol_per_molecule = v_mol/avogadro #m^3
vol_per_molecule_nm3 = vol_per_molecule*10**27 #nm^3
molecule_per_nm3 = 1/vol_per_molecule_nm3 #nm^-3
molecules_per_lenght = (molecule_per_nm3)**(1/3) #nm^-1


############## reservoir ######################################################


nz_planes = 5
k1 = 4 # adding const
z_vec = np.linspace(0,lenght_reservoir-delta_wall,nz_planes)
modified_molecules_per_lenght = (molecule_per_nm3/nz_planes)**(1/2)
x_vec = np.linspace(0,width_reservoir,int(width_reservoir*modified_molecules_per_lenght+k1/2))
y_vec = np.linspace(0,height_reservoir,int(height_reservoir*modified_molecules_per_lenght+k1/2))

# maybe add n random waters to get right density?

################ old version #######################

#k1 = 5 # adding const
#x_vec = np.linspace(0,width_reservoir,int(width_reservoir*molecules_per_lenght)+k1)
#y_vec = np.linspace(0,height_reservoir,int(height_reservoir*molecules_per_lenght))
#z_vec = np.linspace(0,lenght_reservoir-delta_wall,int((lenght_reservoir-delta_wall)*molecules_per_lenght))

####################################################

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
#print(np.shape(reserv1_vec))
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
    wall_vec[2,idz*len(x_vec)*len(y_vec):(idz+1)*len(x_vec)*len(y_vec)].fill(z+lenght_reservoir)
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



########################### water shift vectors ###############################

bond_lenght = 0.9572 # Angstroem
bond_angle = 104.52
alpha = 180-104.52
deltax = cos(alpha*pi/180)*bond_lenght/10
deltay = sin(alpha*pi/180)*bond_lenght/10

h1_shift_vec = np.array((0.952/10,0,0))
h2_shift_vec = np.array((-deltax,-deltay,0))



####################### XMol input ############################################

def ReadXMol(file):
    cnt_x = []
    cnt_y = []
    cnt_z = []
    with open(file) as f:
        for line in f.readlines()[2:-1]:
            #print(line)
            words = line.split(" ")
            #print(words)
            num_list = []
            for word in words:
                try:
                    if word[:-1] == "n" :
                        num_list.append(float(word[:-2]))
                    else:
                        num_list.append(float(word))
                except:
                    print("fail")
                    
            #print(num_list)
            cnt_x.append(num_list[0])
            cnt_y.append(num_list[1])
            cnt_z.append(num_list[2])
    return cnt_x, cnt_y, cnt_z


########################## cnt ################################################

cnt = ReadXMol('Nanotube_1nm.xmol') # Nanotube_CC.xmol

cnt_x = np.array(cnt[0])
cnt_y = np.array(cnt[1])
cnt_z = np.array(cnt[2])
        
cnt_x = cnt_x/10.0
cnt_y = cnt_y/10.0
cnt_z = cnt_z/10.0

cnt_x += width_reservoir/2
cnt_y += width_reservoir/2
cnt_z += lenght_reservoir 

cnt_vec = np.stack((cnt_x, cnt_y, cnt_z))

########################## graphene sheet #####################################

sheet = ReadXMol('NanotubeSheet_CC.xmol')

sheet_x = np.array(sheet[0])
sheet_y = np.array(sheet[1])
sheet_z = np.array(sheet[2])
        
sheet_x = sheet_x/10.0
sheet_y = sheet_y/10.0
sheet_z = sheet_z/10.0


delete_list = []



for idz, z in enumerate(sheet_x):
    if (sheet_x[idz]-2.5)**2 + (sheet_y[idz]-2.5)**2 <= radius_channel**2:
        delete_list.append(idz)
            

# some ugly stuff

sheet_x = np.delete(sheet_x, delete_list, axis = 0)
sheet_y = np.delete(sheet_y, delete_list, axis = 0)
sheet_z = np.delete(sheet_z, delete_list, axis = 0)

sheet1_z = sheet_z+lenght_reservoir 
sheet2_z = sheet_z+lenght_reservoir+lenght_channel


sheet1_vec = np.stack((sheet_x, sheet_y, sheet1_z))
sheet2_vec = np.stack((sheet_x, sheet_y, sheet2_z))


atom_name_list = ["C"]*(2*len(sheet_z)+len(cnt_z))

carbon_vec = np.concatenate((cnt_vec,sheet1_vec,sheet2_vec), axis = 1)

zeros = np.zeros_like(carbon_vec)

carbon_vec = np.concatenate((carbon_vec,zeros), axis = 0)



########################### rotate waters #####################################


oxygen_vec = np.concatenate((channel_vec, reserv1_vec, reserv2_vec), axis = 1)

hydrogen1_vec = np.zeros_like(oxygen_vec) 
hydrogen2_vec = np.zeros_like(oxygen_vec)

result = np.zeros_like(m) 
for i in range(np.shape(oxygen_vec)[1]-1):
    
    alph = (np.random.random((1))*2*np.pi-np.pi)[0]
    beta = (np.random.random((1))*2*np.pi-np.pi)[0]
    gamm = (np.random.random((1))*2*np.pi-np.pi)[0]
    
    #alph = (np.random.random((1))*np.pi)[0]
    #beta = (np.random.random((1))*np.pi)[0]
    #gamm = (np.random.random((1))*np.pi)[0]
    
    rotation_mat = np.array([[np.cos(alph)*np.cos(beta),
                             np.cos(alph)*np.sin(beta)*np.sin(gamm)-np.sin(alph)*np.cos(gamm),
                             np.cos(alph)*np.sin(beta)*np.cos(gamm)+np.sin(alph)*np.sin(gamm)
                             ],
                            [np.sin(alph)*np.cos(beta),
                             np.sin(alph)*np.sin(beta)*np.sin(gamm)+np.cos(alph)*np.cos(gamm),
                             np.sin(alph)*np.sin(beta)*np.cos(gamm)-np.cos(alph)*np.sin(gamm)
                             ],
                            [-np.sin(beta),
                             np.cos(beta)*np.sin(gamm),
                             np.cos(beta)*np.cos(gamm)
                             ]])
    
    h1_shift_vec = np.dot(rotation_mat.T,h1_shift_vec)
    h2_shift_vec = np.dot(rotation_mat.T,h2_shift_vec)
    
    hydrogen1_vec[:3, i] = oxygen_vec[:3, i] + h1_shift_vec
    hydrogen2_vec[:3, i] = oxygen_vec[:3, i] + h2_shift_vec
  


# number of steps & number of atoms
n_steps = 1
n_mol = np.shape(oxygen_vec)[1]
n_carb = np.shape(carbon_vec)[1]
n_atoms =  3*n_mol+n_carb # change to right 


atom_names = ["OW","HW1","HW2"]
atom_name_list = atom_names*n_atoms

molecule1_name = "WAT"
molecule2_name = "CAR"


def MoleculeDataToStringList(molecule_data, molecule_name, atom_name):
    num_mol = np.shape(molecule_data)[1]
    molecule_list = list(np.linspace(1,num_mol, num_mol))
    molecule_name_list = []
    atom_name_list = []
    idx_list = []
    
    for mol in molecule_list:
        molecule_name_list.append(molecule_name)
        atom_name_list.append(atom_name)
        idx_list.append((int(mol)))
    return (molecule_name_list, atom_name_list, idx_list)


OXY = (MoleculeDataToStringList(oxygen_vec, molecule1_name, atom_names[0]))
H1 = (MoleculeDataToStringList(hydrogen1_vec, molecule1_name, atom_names[1]))
H2 = (MoleculeDataToStringList(hydrogen2_vec, molecule1_name, atom_names[2]))
Carbon = (MoleculeDataToStringList(carbon_vec, molecule2_name, "C"))

sorted_vec = np.zeros_like(np.concatenate((oxygen_vec, hydrogen1_vec, hydrogen2_vec), axis = 1))

full_molname_vec = []
full_atomname_vec = []
full_idx_vec = []
full_vec = np.concatenate((oxygen_vec, hydrogen1_vec, hydrogen2_vec, carbon_vec), axis = 1)

for ind, molecule in enumerate(oxygen_vec[0,:]):
    
    sorted_vec[:,ind*3] = oxygen_vec[:,ind]
    sorted_vec[:,ind*3+1] = hydrogen1_vec[:,ind]
    sorted_vec[:,ind*3+2] = hydrogen2_vec[:,ind]
    
    full_molname_vec.append(molecule1_name)
    full_molname_vec.append(molecule1_name)
    full_molname_vec.append(molecule1_name)
    full_atomname_vec.append(atom_names[0])
    full_atomname_vec.append(atom_names[1])
    full_atomname_vec.append(atom_names[2])
    full_idx_vec.append((ind))
    full_idx_vec.append((ind))
    full_idx_vec.append((ind))

sorted_vec = np.concatenate((sorted_vec, carbon_vec), axis = 1)

full_molname_vec = full_molname_vec+Carbon[0]
full_atomname_vec = full_atomname_vec+Carbon[1]
full_idx_vec = full_idx_vec+[1]*len(Carbon[2])


print("full_vec lenght:")
print(np.shape(full_vec)[1])
print("full_molname_vec lenght:")
print(len(full_molname_vec))
print("full_atomname_vec lenght:")
print(len(full_atomname_vec))


########################## plotting #########################################

fig = plt.figure()
ax = plt.axes(projection='3d')

#ax.scatter(reserv1_vec[0,:], reserv1_vec[1,:], reserv1_vec[2,:], c='royalblue', marker='o', s = 2)
#ax.scatter(reserv2_vec[0,:], reserv2_vec[1,:], reserv2_vec[2,:], c='firebrick', marker='o', s = 2)
#ax.scatter(wall_vec[0,:], wall_vec[1,:], wall_vec[2,:], c='gray', marker='o', s = 2)
#ax.scatter(channel_vec[0,:], channel_vec[1,:], channel_vec[2,:], c='k', marker='o', s = 2)
#ax.scatter(full_vec[0,:], full_vec[1,:], full_vec[2,:], c='k', marker='o', s = 2)

ax.scatter(full_vec[0,:2528], full_vec[1,:2528], full_vec[2,:2528], c='royalblue', marker='o', s = 2)
ax.scatter(full_vec[0,2528:2528*2], full_vec[1,2528:2528*2], full_vec[2,2528:2528*2], c='firebrick', marker='o', s = 1)
ax.scatter(full_vec[0,2528*2:2528*2+len(channel_vec)], full_vec[1,2528*2+len(channel_vec)], full_vec[2,2528*2:2528*2+len(channel_vec)], c='firebrick', marker='o', s = 1)
ax.scatter(full_vec[0,2528*2+len(channel_vec):], full_vec[1,2528*2+len(channel_vec):], full_vec[2,2528*2+len(channel_vec):], c='green', marker='+', s = 2)
#ax.view_init(85, 50)
ax.view_init(20, 50)

###############################################################################


x = sorted_vec[:3,:]
atom_name = full_atomname_vec
molecule_name = full_molname_vec
molecule_number = full_idx_vec

box = np.zeros(3)

# assign arbitrary values for now
box[0] = 10e0
box[1] = 10e0
box[2] = 10e0

title = "test"


########################### output ############################################
###############################################################################


########################### gro file create ###################################

f = open("demofile2.gro", "w")

for t in range(n_steps):

  f.write("%-10s%2s%12.4f\n" % (title,"t=",t))
  f.write("%5d\n" % (n_atoms))

  for i in range(n_atoms):
    #print(i)
    # format: (i5,2a5,i5,3f8.3,3f8.4) (the final 3 values can be used for velocities)
    f.write("%5d%-5s%-5s%5d%8.3f%8.3f%8.3f\n" % (molecule_number[i],molecule_name[i],atom_name[i],i,x[0,i],x[1,i],x[2,i]))

  f.write("%10.5f%10.5f%10.5f\n" % (box[0], box[1], box[2]))

f.close()




############################# ndx file create #################################

name_list = ["System", "Water", "OW", "HW1_HW2", "Carbon"]


sys_lenght = len(full_atomname_vec)
sys_ind_array = np.linspace(1,sys_lenght, sys_lenght)

water_num = np.shape(oxygen_vec)[1]+np.shape(hydrogen1_vec)[1]+np.shape(hydrogen2_vec)[1]
water_ind_array = np.linspace(1,water_num, water_num)

ow_num = int((water_num)/3)
ow_ind_array = np.linspace(1,ow_num, ow_num)

h12_ind_array = np.linspace(ow_num+1,water_num, ow_num*2)

ind_list = [sys_ind_array, water_ind_array, ow_ind_array, h12_ind_array]



#f = open("test.ndx", "w")
#
#for ind, name in enumerate(name_list):
#
#    f.write("[ "+name+" ]\n")
#
#    for i in ind_list[ind]:
#        f.write(" "+str(int(i))+" ")
#
#    f.write("\n")
  
  
#f.close()



########################### itp file create ###################################



residue = ["PHPC"]*np.shape(sorted_vec)[1]




