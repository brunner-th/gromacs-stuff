# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 16:41:50 2023

@author: tomal
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from math import sin, cos, pi
from scipy.optimize import curve_fit
from scipy.interpolate import interp2d


plt.rcParams['figure.figsize'] = [8, 6]
plt.rcParams['figure.dpi'] = 400
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

#################### input ####################################################

n_oxy = 2000
n_steps = 10
nbins = 50
deltabin = 2

def scatter_hist(x, y, ax, ax_histx, ax_histy):
    
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    ax.scatter(x, y,marker = ".", s = 1)
    
    xmax = np.max(x)+deltabin
    ymax = np.max(y)+deltabin
    xmin = np.min(x)-deltabin
    ymin = np.min(y)-deltabin
    
    outx = ax_histx.hist(x, nbins, (xmin,xmax))
    outy = ax_histy.hist(y, nbins, (ymin,ymax) ,orientation='horizontal')
    
    nx = outx[0]
    ny = outy[0]
    binx = outx[1]
    biny = outy[1]
    
    
    return nx, ny, binx, biny
    



def ReadXTC(file):
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





#xyz_oxy = np.zeros((3,n_oxy, n_steps))
xyz_oxy = np.random.rand(3,n_oxy, n_steps)+ np.random.default_rng().normal(loc=0.5, scale=5, size=(3,n_oxy,n_steps))


uxyz_oxy = xyz_oxy[:,:,:-1]-xyz_oxy[:,:,1:]
u_oxy = np.sqrt(uxyz_oxy[0,:,:]**2+uxyz_oxy[1,:,:]**2+uxyz_oxy[2,:,:]**2)
uz_oxy = uxyz_oxy[2,:,:]+10



######################### hist ################################################



fig = plt.figure(constrained_layout=True)
ax = fig.add_gridspec(top=0.75, right=0.75).subplots()
ax.set(aspect=0.8)
ax_histx = ax.inset_axes([0, 1.05, 1, 0.25], sharex=ax)
ax_histy = ax.inset_axes([1.05, 0, 0.25, 1], sharey=ax)





###############################################################################


x_dist = [[],[]]
y_dist = [[],[]]

for ind, xyz in enumerate(u_oxy[0,:]):
    
    out = scatter_hist(xyz_oxy[0,:,ind], uz_oxy[:,ind], ax, ax_histx, ax_histy)
    
    x_dist[0].append(out[0])
    x_dist[1].append(out[2])
    y_dist[0].append(out[1])
    y_dist[1].append(out[3])
    
    #scatter_hist(xyz_oxy[0,:,ind], u_oxy[:,ind], ax, ax_histx, ax_histy)
    #plt.plot(xyz_oxy[0,:,ind], uz_oxy[:,ind], marker = ".", linestyle = " ", markersize = 1)

plt.show()


for ind, dist in enumerate(x_dist[0]):
    xdata = (x_dist[1][ind][:-1]-x_dist[1][ind][1:])/2+x_dist[1][ind][:-1]
    plt.plot(xdata, dist, marker = "+",  linestyle = " ")



##################### fitting #################################################


def Func(data, a, b, c):
    return a*data**2+b*data+c


initial_guess = np.array([-100,0, 200])
boun = [(-1000,-1000,-1000), (1000, 1000,1000)]
xdata = xdata
ydata = dist

params, pcov = curve_fit(Func, xdata, ydata, initial_guess, bounds=boun)

print("fit params:")
print(params)


plt.plot(xdata, Func(xdata, params[0], params[1], params[2]), color = "k")





