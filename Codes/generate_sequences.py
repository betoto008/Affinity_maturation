import sys
import numpy as np
import matplotlib.pyplot as plt
from Immuno_models import*
from Bio import Phylo
from io import StringIO
from matplotlib.lines import Line2D
from datetime import datetime, timedelta
import scipy.special as sc
import os.path
import pickle
from matplotlib import style
from scipy.interpolate import interp1d

M2 = np.loadtxt('../Text_files/MJ2.txt', skiprows= 1, usecols=range(1,21)).tolist()

n_seq = input('How many sequences do you want to create? ')
Sequences = generate_Sequences(int(n_seq), Energy_Matrix = M2)
pickle.dump( Sequences, open( "../Text_files/Sequences-n_seq-%d.pkl"%(int(n_seq)), "wb" ) )