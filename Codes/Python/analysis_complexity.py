#---------Import---------
import sys
sys.path.append('../Codes/')
sys.path.append('../Codes/Python/')
from Immuno_models import*
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.ticker as mtick
from matplotlib.ticker import PercentFormatter
from matplotlib.lines import Line2D


#------------------------------------

Text_files_path = '../../../../../Dropbox/Research/Evolution_Immune_System/Text_files/Complexity/'

colors = plt.cm.plasma(np.linspace(0,1,10))
T = 8
dt = 1
Ls = [10, 15, 20, 50]
Types = ["MM", "MJ"]
M_epitopes = 10

for Type in Types:
    if(Type == "MJ"):
        T = 150
    elif(Type == "MM"):
        T = 8
    for L in Ls:
        fig, ax = plt.subplots(figsize = (12,10), gridspec_kw={'top':.95, 'left':.18    })
        for i in range(M_epitopes):
            data = np.loadtxt(Text_files_path + 'state_L-%d_L_alphabet-20_N_epitopes-%d_'%(L, i+1)+Type+'.txt')
            time = np.linspace(1*dt, T*dt, T)
            ax.plot(time, data, color = colors[i], linestyle = '--', label = "%d"%(i+1))
            #ax.plot(time, (1-(1-np.exp(-time))**(i+1)), color = colors[i])
        my_plot_layout(ax=ax, xlabel = 'Steps', ylabel = "Retained memory", x_fontsize = 24, y_fontsize = 24)
        ax.set_ylim(-.1, 1.1)
        ax.legend(loc = 5, fontsize = 24)
        fig.savefig('../../Figures/4_Complexity/memory_lost_L-%d_'%(L)+Type+'.png')