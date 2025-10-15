#
# clear
#
from IPython import get_ipython
get_ipython().magic('reset -sf');
import os
os.system('clear');
#
import pandas as pd;
import numpy as np;
from matplotlib import pyplot as plt    
#
L0 = 2.0*np.pi;
N = 512;
tot_row = 18;
istep = 30838;
istep_c = f'{istep:09d}'
pos = 1;
pos_c = f'{pos:03d}'
work_dir  = "/home/nn8802/Documents/simulation/basilisk_dotc_turb/windwave_turb/";
#ux_bin = np.fromfile(work_dir+'slices/ux_2d_'+istep_c+'_'+pos_c+'.bin');
uxA_bin = np.fromfile(work_dir+'field/ux_3d_'+istep_c+'.bin');
uxA_bin = np.transpose(uxA_bin.reshape([N, N, N]));
#
x_int = np.linspace(-L0/2,L0/2,N,endpoint=False)+L0/N/2;
y_int = np.linspace(0,L0,N,endpoint=False)+L0/N/2;
x_til, y_til = np.meshgrid(x_int,y_int);
#
fig = plt.figure(figsize=(1.1*10,5));
plt.contourf(x_til,y_til,np.transpose(np.mean(uxA_bin,axis=-1)),100);
#plt.ylim([0.90,1.20]);
