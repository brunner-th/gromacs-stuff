# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 16:44:48 2022

@author: tomal
"""
import re
import numpy as np

cnt_name_list = []
cnt_x = []
cnt_y = []
cnt_z = []














with open('NanotubeSheet_CC.xmol') as f:
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
        
    
                

cnt_x = np.array(cnt_x)
cnt_y = np.array(cnt_y)
cnt_z = np.array(cnt_z)
        
cnt_x = cnt_x/10.0
cnt_y = cnt_y/10.0
cnt_z = cnt_z/10.0

cnt_x += 0
cnt_y += 0
cnt_z += 1.5+5

        
####################### plane ##################################################

delete_list = []


for idz, z in enumerate(cnt_x):
    if (cnt_x[idz]-2.5)**2 + (cnt_y[idz]-2.5)**2 <= 1**2:
        delete_list.append(idz)
            

# some ugly stuff



cnt_x = np.delete(cnt_x, delete_list, axis = 0)
cnt_y = np.delete(cnt_y, delete_list, axis = 0)
cnt_z = np.delete(cnt_z, delete_list, axis = 0)



###############################################################################


n_steps = 1

n_atoms = np.shape(cnt_x)[0]

molecule_number = np.linspace(1,n_atoms, n_atoms)

box = [10,10,10]
title = "cnt"
atom_name = "C"
molecule_name = "CNT"

f = open("sheet2.gro", "w")

for t in range(n_steps):

  f.write("%-10s%2s%12.4f\n" % (title,"t=",t))
  f.write("%5d\n" % (n_atoms))

  for i in range(n_atoms):

    
    
    # format: (i5,2a5,i5,3f8.3,3f8.4) (the final 3 values can be used for velocities)
    f.write("%5d%-5s%-5s%5d%8.3f%8.3f%8.3f\n" % (molecule_number[i],molecule_name,atom_name,i,cnt_x[i],cnt_y[i],cnt_z[i]))

  f.write("%10.5f%10.5f%10.5f\n" % (box[0], box[1], box[2]))

f.close()