#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
import pandas as pd;
import numpy as np;
from matplotlib import pyplot as plt    
#
L0 = 2.0*np.pi;
N = 512;
tot_row = 18;
istep = 621442;
istep_c = f'{istep:09d}'
pos = 1;
pos_c = f'{pos:03d}'
work_dir  = "/home/nn8802/Documents/simulation/basilisk_dotc_turb/windwave_turb/";
ux_bin = np.fromfile(work_dir+'slices/ux_2d_'+istep_c+'_'+pos_c+'.bin');
ux_bin = np.transpose(ux_bin.reshape([N, N]));
fv_bin = np.fromfile(work_dir+'slices/fv_2d_'+istep_c+'_'+pos_c+'.bin');
fv_bin = np.transpose(fv_bin.reshape([N, N]));
#
x_int = np.linspace(-L0/2,L0/2,N,endpoint=False)+L0/N/2;
y_int = np.linspace(0,L0,N,endpoint=False)+L0/N/2;
x_til, y_til = np.meshgrid(x_int,y_int);
#
plt.contourf(x_til,y_til,ux_bin,100);
plt.contour(x_til,y_til,fv_bin,[0.5]);
plt.ylim([0.90,1.20]);