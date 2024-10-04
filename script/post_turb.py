#
# clear
#
#from IPython import get_ipython
#get_ipython().magic('reset -sf');
import os
#os.system('clear');
#
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from my_funcs import interp_2d,phase_partion,return_file,mom_flux_p,mom_flux_v,ene_flux_p,ene_flux_v,cart_to_wf,get_amp
from my_funcs import mom_flux_p_alt, mom_flux_v_alt
#
L0    = 2.0*np.pi; h = 1; u_s = 0.25;
N     = 512; k_ = 4; 
rho_r = 1000/1.225; Re_tau = 720;
rho1  = 1; rho2 = rho1/rho_r;
mu2   = rho2*u_s*(L0-h)/Re_tau;
tot_row = 18;
#
os.system('rm -rf __pycache__ && rm -f *.out && rm -rf wave_coord && mkdir wave_coord');
#
work_dir = '/home/nn8802/Documents/document_postdoc/proj_wave_break/turbulent_results/uoc_0p5/';
#work_dir = '/scratch/gpfs/ns8802/multiphase_case/re720_bo0200_ak0p3_uoc0p5_L10/';
#
time_fld = pd.read_csv(work_dir+'field/log_field.out', header=None, sep=" ");
time_fld = time_fld.to_numpy();
time_eta = pd.read_csv(work_dir+'eta/global_int.out', header=None, sep=" ");
time_eta = time_eta.to_numpy();
#
x_int = np.linspace(-L0/2,L0/2,N,endpoint=False)+L0/N/2;
z_int = np.linspace(-L0/2,L0/2,N,endpoint=False)+L0/N/2;
x_til, z_til = np.meshgrid(x_int,z_int);
#
print(work_dir);
for i in range(len(time_fld)):
    #
    # define the time and row
    #
    time    = time_fld[i,0];
    istep   = int(time_fld[i,1]);
    istep_c = f'{istep:09d}';
    print("****************");
    print(istep, time, i);
    print("****************");
    #
    # load eta_loc
    #
    etalo = np.fromfile(work_dir+'eta/eta_loc/eta_loc_t'+istep_c+'.bin')
    size  = etalo.shape;
    tot_row_i = int(size[0]/tot_row);
    print(tot_row_i);
    etalo = etalo.reshape([tot_row_i, tot_row]);
    #
    # we remove
    #
    print("First pass of remove")
    eta_m0  = 1.0; cirp_th = 0.20;
    new_row = 0;
    for i in range(tot_row_i):
        if ( abs(etalo[i][12]-eta_m0) < cirp_th ):
           new_row += 1;
    #
    print("Second pass of remove")
    etal = np.zeros([new_row, 18]);
    for i in range(new_row):
        if ( abs(etalo[i][12]-eta_m0) < cirp_th ):
           etal[i][:] = etalo[i][:];
    #
    print("Assign array")
    xpo = etal[:,0];  zpo = etal[:,1];
    pre = etal[:,2]; 
    Sxx = etal[:,3];  Syy = etal[:,4];  Szz = etal[:,5]; 
    Sxy = etal[:,6];  Sxz = etal[:,7];  Syz = etal[:,8]; 
    uxi = etal[:,9];  uyi = etal[:,10]; uzi = etal[:,11]; 
    eta = etal[:,12]; eps = etal[:,13];
    n_x = etal[:,14]; n_y = etal[:,15]; n_z = etal[:,16]; 
    #
    # we interpolate eta on a 2D cartesian grid with equidistant spacing equal to the printing resolution (2**9)
    #
    print("Interpolation to a Cartesian grid")
    pre_int = interp_2d(xpo, zpo, x_til, z_til, pre);
    Sxx_int = interp_2d(xpo, zpo, x_til, z_til, Sxx);
    Syy_int = interp_2d(xpo, zpo, x_til, z_til, Syy);
    Szz_int = interp_2d(xpo, zpo, x_til, z_til, Szz);
    Sxy_int = interp_2d(xpo, zpo, x_til, z_til, Sxy);
    Sxz_int = interp_2d(xpo, zpo, x_til, z_til, Sxz);
    Syz_int = interp_2d(xpo, zpo, x_til, z_til, Syz);
    uxi_int = interp_2d(xpo, zpo, x_til, z_til, uxi);
    uyi_int = interp_2d(xpo, zpo, x_til, z_til, uyi);
    uzi_int = interp_2d(xpo, zpo, x_til, z_til, uzi);
    eta_int = interp_2d(xpo, zpo, x_til, z_til, eta);
    eps_int = interp_2d(xpo, zpo, x_til, z_til, eps);
    n_x_int = interp_2d(xpo, zpo, x_til, z_til, n_x);
    n_y_int = interp_2d(xpo, zpo, x_til, z_til, n_y);
    n_z_int = interp_2d(xpo, zpo, x_til, z_til, n_z);   
    #
    # compute momentum and energy fluxes
    #
    print("Compute momentum flux - pressure")
    [mf_px,mf_py,mf_pz] = mom_flux_p(pre_int,n_x_int,n_y_int,n_z_int);
    print("Compute momentum flux - viscous dissipation")
    [mf_vx,mf_vy,mf_vz] = mom_flux_v(Sxx_int,Sxy_int,Sxz_int,Syy_int,Syz_int,Szz_int,n_x_int,n_y_int,n_z_int,mu2);
    print("Energy flux - pressure")
    en_p = ene_flux_p(pre_int,uxi_int,uyi_int,uzi_int,n_x_int,n_y_int,n_z_int);
    print("Energy flux - viscous dissipation")
    en_v = ene_flux_v(Sxx_int,Sxy_int,Sxz_int,Syy_int,Syz_int,Szz_int,uxi_int,uyi_int,uzi_int,n_x_int,n_y_int,n_z_int,mu2);
    #
    # compute stress budget (pressure and viscous term)
    #
    eta_1d = np.average(eta_int,axis=0)-np.average(eta);
    pre_1d = np.average(pre_int,axis=0)-np.average(pre);
    Sxx_1d = np.average(Sxx_int,axis=0);
    Sxy_1d = np.average(Sxy_int,axis=0);
    mf_px_alt = mom_flux_p_alt(pre_1d,eta_1d,L0,N);
    mf_vx_alt = mom_flux_v_alt(Sxx_1d,Sxy_1d,eta_1d,L0,N,mu2);
    #
    # compute amplitude
    #
    ak = get_amp(eta_int,k_,L0);
    #
    # load the 2d span-averaged binary files
    #
    print("Load field and phase partition")
    fv_2d = return_file(work_dir,'fv',istep_c,N,1);
    ux_2d = return_file(work_dir,'ux',istep_c,N,1);
    uy_2d = return_file(work_dir,'uy',istep_c,N,1);
    pr_2d = return_file(work_dir,'pr',istep_c,N,1);
    di_2d = return_file(work_dir,'di',istep_c,N,1);
    rh_2d = return_file(work_dir,'ru',istep_c,N,1);
    md_2d = return_file(work_dir,'md',istep_c,N,1);
    [ux_2d_air,ux_2d_wat,uy_2d_air,uy_2d_wat] = phase_partion(ux_2d,uy_2d,fv_2d,1.0,1.0);
    [pr_2d_air,pr_2d_wat,di_2d_air,di_2d_wat] = phase_partion(pr_2d,di_2d,fv_2d,1.0,1.0);
    [rh_2d_air,rh_2d_wat,md_2d_air,md_2d_wat] = phase_partion(rh_2d,md_2d,fv_2d,1.0,1.0);
    #
    eta_m0 = 1;
    print("From cartesian to wf")
    [ux_air_2d_wf, ux_air_1d_wf, zplot_air, zeta_air] = cart_to_wf(ux_2d_air,eta_1d,N,L0,k_,eta_m0); # ux_air, x-vel
    [ux_wat_2d_wf, ux_wat_1d_wf, zplot_wat, zeta_wat] = cart_to_wf(ux_2d_wat,eta_1d,N,L0,k_,eta_m0); # ux_wat, x-vel
    [pr_air_2d_wf, pr_air_1d_wf, zplot_air, zeta_air] = cart_to_wf(pr_2d_air,eta_1d,N,L0,k_,eta_m0); # pr_air, pressure
    [pr_wat_2d_wf, pr_wat_1d_wf, zplot_wat, zeta_wat] = cart_to_wf(pr_2d_wat,eta_1d,N,L0,k_,eta_m0); # pr_wat, pressure
    [di_air_2d_wf, di_air_1d_wf, zplot_air, zeta_air] = cart_to_wf(di_2d_air,eta_1d,N,L0,k_,eta_m0); # di_air, dissipation
    [di_wat_2d_wf, di_wat_1d_wf, zplot_wat, zeta_wat] = cart_to_wf(di_2d_wat,eta_1d,N,L0,k_,eta_m0); # di_wat, dissipation
    [rh_air_2d_wf, rh_air_1d_wf, zplot_air, zeta_air] = cart_to_wf(rh_2d_air,eta_1d,N,L0,k_,eta_m0); # rh_air, rho*uv
    [rh_wat_2d_wf, rh_wat_1d_wf, zplot_wat, zeta_wat] = cart_to_wf(rh_2d_wat,eta_1d,N,L0,k_,eta_m0); # rh_wat, rho*uv
    [md_air_2d_wf, md_air_1d_wf, zplot_air, zeta_air] = cart_to_wf(md_2d_air,eta_1d,N,L0,k_,eta_m0); # md_air, dudy
    [md_wat_2d_wf, md_wat_1d_wf, zplot_wat, zeta_wat] = cart_to_wf(md_2d_wat,eta_1d,N,L0,k_,eta_m0); # md_wat, dudy
    #
    print("Final print of glo_obs")
    f = open("glo_obs_post.out","a");
    f.write("%.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f \n" %(1.0*istep, time, ak, mf_px, mf_py, mf_pz, en_p, en_v));
    f.flush();
    f.close();
    #
    print("Final print of glo_obs alt")
    f = open("glo_obs_post_alt.out","a");
    f.write("%.15f %.15f %.15f %.15f \n" %(1.0*istep, time, mf_px_alt, mf_vx_alt));
    f.flush();
    f.close();
    #
    print("Final print of wf")
    f = open("wave_coord/prof_wf_"+istep_c+".out","w");
    for i in range(N):
        f.write("%.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f\n" %(zeta_air[i], zeta_wat[i], ux_air_1d_wf[i], ux_wat_1d_wf[i], 
                                                                                                                        pr_air_1d_wf[i], pr_wat_1d_wf[i], 
                                                                                                                        di_air_1d_wf[i], di_wat_1d_wf[i],
                                                                                                                        rh_air_1d_wf[i], rh_wat_1d_wf[i],
                                                                                                                        md_air_1d_wf[i], md_wat_1d_wf[i]) );
    f.flush();
    f.close();
