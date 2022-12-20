#!/bin/bash
gmx grompp -c water.gro -f water.mdp -p water.top -o water.tpr -n water.ndx
gmx mdrun -nt 4 -v -s water.tpr -c water.gro -g water.log -o water.trr

