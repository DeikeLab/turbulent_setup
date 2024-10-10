#
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from specgen import spectrum_PM, spectrum_gen_linear, eta_random, spectrum_PM_piece
#
# L = 200;    # domain size
# P = 0.02;        # energy level (estimated so that kpHs is reasonable)
#P = [0.0016/4,0.0016,0.0016*4,0.0016*8.6,0.0016*8.6*2];
P = 0.0016;          # energy level (estimated so that kpHs is reasonable)
Bo_1 = 200
Bo_2 = 1
N = 4
L = 2.0*np.pi    # domain size
kp = 2*np.pi/(L/N)    # peak wavenumber
N_mode =  64    # number of modes
N_power = 5      # directional spreading coeff
k_cf = (kp**2*Bo_1/Bo_2)**0.5
print(k_cf);
#
#
def shape(kmod):
    ''' Choose values here '''
    global P, kp
    #F_kmod = spectrum_PM (P, kp, kmod);
    F_kmod = spectrum_PM_piece (P, kp, k_cf, kmod);
    #F_kmod = spectrum_PM(P, kp, kmod)
    return F_kmod
#
kmod, F_kmod, kx, ky, F_kxky_tile = spectrum_gen_linear(
    shape, N_mode=N_mode, L=L, N_power=N_power)
#
''' Generate a grid in x-y to visualize random eta '''
N_grid = 256  # L = 200
x = np.linspace(-L/2, L/2, N_grid)
y = np.linspace(-L/2, L/2, N_grid)
x_tile, y_tile = np.meshgrid(x, y)
kx_tile, ky_tile = np.meshgrid(kx, ky)
t = 0
eta_tile, phase_tile = eta_random(
    t, kx_tile, ky_tile, F_kxky_tile, x_tile, y_tile)
kpHs = kp*np.std(eta_tile)*4;
print('kpHs = %g' % (kpHs))
#
print('Spectrum array shape:', F_kxky_tile.shape)
path = "";
#
file = "F_kxky_P%.5f_L%gm_N%g_kpHs%g_N_mode%g" % (P, L, N_power, kpHs, N_mode)
fF = open(path + file, 'bw')
print('Written to', path+file)
F_output = F_kxky_tile.astype('float32')
#F_output = F_kxky_tile.astype('double');
F_output.tofile(fF)
#
file = "kx_P%.5f_L%gm_N%g_kpHs%g_N_mode%g" % (P, L, N_power, kpHs, N_mode)
fF = open(path + file, 'bw')
print('Written to', path+file)
F_output = kx.astype('float32')
#F_output = kx.astype('double');
F_output.tofile(fF)
#
file = "ky_P%.5f_L%gm_N%g_kpHs%g_N_mode%g" % (P, L, N_power, kpHs, N_mode)
fF = open(path + file, 'bw')
print('Written to', path+file)
F_output = ky.astype('float32');
#F_output = ky.astype('double');
F_output.tofile(fF)
