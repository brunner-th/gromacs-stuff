import numpy as np

# number of steps & number of atoms
n_steps = 1
n_atoms = 6

# position, molecule names, atom names & simulation box size
x = np.zeros((3,n_atoms))
atom_name = np.empty(n_atoms,dtype=object)
molecule_name = np.empty(n_atoms,dtype=object)
molecule_number = np.arange(n_atoms)
box = np.zeros(3)

# assign arbitrary values for now
box[0] = 1e0
box[1] = 1e0
box[2] = 1e0
molecule_name[0:6] = "Mol"
atom_name[0:3] = "A"
atom_name[3:6] = "B"
title = "title"

f = open("demofile.gro", "w")

for t in range(n_steps):

  f.write("%-10s%2s%12.4f\n" % (title,"t=",t))
  f.write("%5d\n" % (n_atoms))

  for i in range(n_atoms):

    x[:,i] = np.random.random(3)

    # format: (i5,2a5,i5,3f8.3,3f8.4) (the final 3 values can be used for velocities)
    f.write("%5d%-5s%-5s%5d%8.3f%8.3f%8.3f\n" % (molecule_number[i],molecule_name[i],atom_name[i],i,x[0,i],x[1,i],x[2,i]))

  f.write("%10.5f%10.5f%10.5f\n" % (box[0], box[1], box[2]))

f.close()
