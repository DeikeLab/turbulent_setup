#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
import pandas as pd;
import numpy as np;
from matplotlib import pyplot as plt    
from my_funcs import get_amp, interp_2d;
#
tot_row = 18;
work_dir  = "/home/nn8802/Documents/simulation/basilisk_dotc_turb/windwave_turb/";
#
istep   = 30838;
istep_c = f'{istep:09d}';
etalo = np.fromfile(work_dir+'eta/eta_loc/eta_loc_t'+istep_c+'.bin')
size  = etalo.shape;
tot_row_i = int(size[0]/tot_row);
print(tot_row_i);
etalo = etalo.reshape([tot_row_i, tot_row]);
#
L0 = 2.0*np.pi;
k  = 4; N = 512;
eta = etalo[:, 12];
pre = etalo[:, 2];
xpo = etalo[:,0];  
zpo = etalo[:,1];
x_int = np.linspace(-L0/2,L0/2,N,endpoint=False)+L0/N/2;
z_int = np.linspace(-L0/2,L0/2,N,endpoint=False)+L0/N/2;
x_til, z_til = np.meshgrid(x_int,z_int);
eta_int = interp_2d(xpo, zpo, x_til, z_til, eta);
pre_int = interp_2d(xpo, zpo, x_til, z_til, pre);
#plt.contourf(x_til,z_til,eta_int)
plt.contourf(x_til,z_til,pre_int-np.average(pre_int))
#plt.plot(x_int,np.average(eta_int,axis=0));
#plt.plot(x_int,eta_int[1])
#


#
ak = get_amp(eta_int,k,L0);