#---------Import---------
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.ticker as mtick
from matplotlib.ticker import PercentFormatter
from matplotlib.lines import Line2D

#---------Functions---------
def my_linear_func(x, a, b):
    return a + b*x
def my_quadratic_func(x, a, b, c):
    return np.log(a)+np.log(np.sqrt(-b)) + b*(x-c)**2
def my_plot_layout(ax, yscale = 'linear', xscale = 'linear', ticks_labelsize = 24,
                   xlabel = '', ylabel = '', title = '', x_fontsize=24, y_fontsize = 24,
                   t_fontsize = 24):
    ax.tick_params(labelsize = ticks_labelsize)
    ax.set_yscale(yscale)
    ax.set_xscale(xscale)
    ax.set_xlabel(xlabel, fontsize = x_fontsize)
    ax.set_ylabel(ylabel, fontsize = y_fontsize)
    ax.set_title(title, fontsize = y_fontsize)
#------------------------------------

Text_files_path = '../../../../../Dropbox/Research/Evolution_Immune_System/Text_files/Complexity/'

fig, ax = plt.subplots(figsize = (12,10), gridspec_kw={'top':.95, 'left':.18    })
colors = plt.cm.plasma(np.linspace(0,1,10))

T = 8
dt = 1
L = 50
Type = "MM"
M_epitopes = 10

for i in range(M_epitopes):
    data = np.loadtxt(Text_files_path + 'state_L-%d_L_alphabet-20_N_epitopes-%d_'%(L, i+1)+Type+'.txt')
    time = np.linspace(1*dt, T*dt, T)
    ax.plot(time, data, color = colors[i], linestyle = '--', label = "%d"%(i+1))
    #ax.plot(time, (1-(1-np.exp(-time))**(i+1)), color = colors[i])
my_plot_layout(ax=ax, xlabel = 'Steps', ylabel = "Memory lost", x_fontsize = 24, y_fontsize = 24)
ax.legend(loc = 5, fontsize = 24)

fig.savefig('../../Figures/4_Complexity/MJ_memory_lost_L-%d_'%(L)+Type+'.png')