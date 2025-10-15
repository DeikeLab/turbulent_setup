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
from my_funcs import interp_2d,return_field,split_profile,get_amp
from my_funcs import mapping_function_1,mapping_function_2,mapping_function_3,mapping_function_4,transform_to_wave_following
#
L0    = 2.0*np.pi; h = 1; u_s = 0.25;
N     = 512; r_L0lam = 4;
lam   = L0/r_L0lam;
k_    = 2.0*np.pi/lam; 
rho_r = 1000/1.225; Re_tau = 720;
rho1  = 1; rho2 = rho1/rho_r;
mu2   = rho2*u_s*(L0-h)/Re_tau;
tot_col = 18;
#
os.system('rm -rf __pycache__ && rm -f *.out && rm -rf wave_coord && mkdir wave_coord');
#
#work_dir = '/scratch/cimes/ns8802/multiphase_case_be/re720_bo0200_ak0p30_uoc0p9_L10/';
work_dir = '/home/nn8802/Documents/simulation/basilisk_dotc_turb/script/'
#
time_fld = pd.read_csv(work_dir+'field/log_field.out', header=None, sep=" ");
time_fld = time_fld.to_numpy();
time_eta = pd.read_csv(work_dir+'eta/global_int.out', header=None, sep=" ");
time_eta = time_eta.to_numpy();
#
# we follow the convection employed in the literature
#  x: streawise, y: spanwise, z: wall-normal.
#
x_int = np.linspace(-L0/2,L0/2,N,endpoint=False)+L0/N/2;
y_int = np.linspace(-L0/2,L0/2,N,endpoint=False)+L0/N/2;
z_int = np.linspace(    0,L0  ,N,endpoint=False)+L0/N/2;
x_til, y_til = np.meshgrid(x_int, y_int);
#
print(work_dir);
for i in range(len(time_fld)):
    #
    # define the time and row
    #
    pos = 1;
    time    = time_fld[i,0];
    istep   = int(time_fld[i,1]);
    istep_c = f'{istep:09d}';
    pos_c   = f'{pos:03d}';
    print("****************");
    print(istep, time, i);
    print("****************");
    #
    # load eta_loc
    #
    etalo = np.fromfile(work_dir+'eta/eta_loc/eta_loc_t'+istep_c+'.bin')
    size  = etalo.shape;
    tot_row_i = int(size[0]/tot_col);
    print(tot_row_i);
    etalo = etalo.reshape([tot_row_i, tot_col]);
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
    pre_int = interp_2d(xpo, zpo, x_til, y_til, pre);
    Sxx_int = interp_2d(xpo, zpo, x_til, y_til, Sxx);
    Syy_int = interp_2d(xpo, zpo, x_til, y_til, Syy);
    Szz_int = interp_2d(xpo, zpo, x_til, y_til, Szz);
    Sxy_int = interp_2d(xpo, zpo, x_til, y_til, Sxy);
    Sxz_int = interp_2d(xpo, zpo, x_til, y_til, Sxz);
    Syz_int = interp_2d(xpo, zpo, x_til, y_til, Syz);
    uxi_int = interp_2d(xpo, zpo, x_til, y_til, uxi);
    uyi_int = interp_2d(xpo, zpo, x_til, y_til, uyi);
    uzi_int = interp_2d(xpo, zpo, x_til, y_til, uzi);
    eta_int = interp_2d(xpo, zpo, x_til, y_til, eta);
    eps_int = interp_2d(xpo, zpo, x_til, y_til, eps);
    n_x_int = interp_2d(xpo, zpo, x_til, y_til, n_x);
    n_y_int = interp_2d(xpo, zpo, x_til, y_til, n_y);
    n_z_int = interp_2d(xpo, zpo, x_til, y_til, n_z);   
    #
    # compute momentum and energy fluxes
    #
    #print("Compute momentum flux - pressure")
    #[mf_px,mf_py,mf_pz] = mom_flux_p(pre_int,n_x_int,n_y_int,n_z_int);
    #print("Compute momentum flux - viscous dissipation")
    #[mf_vx,mf_vy,mf_vz] = mom_flux_v(Sxx_int,Sxy_int,Sxz_int,Syy_int,Syz_int,Szz_int,n_x_int,n_y_int,n_z_int,mu2);
    #print("Energy flux - pressure")
    #en_p = ene_flux_p(pre_int,uxi_int,uyi_int,uzi_int,n_x_int,n_y_int,n_z_int);
    #print("Energy flux - viscous dissipation")
    #en_v = ene_flux_v(Sxx_int,Sxy_int,Sxz_int,Syy_int,Syz_int,Szz_int,uxi_int,uyi_int,uzi_int,n_x_int,n_y_int,n_z_int,mu2);
    #
    # compute stress budget (pressure and viscous term)
    #
    eta_1d = np.average(eta_int,axis=0);#-np.average(eta);
    #pre_1d = np.average(pre_int,axis=0)-np.average(pre);
    #Sxx_1d = np.average(Sxx_int,axis=0);
    #Sxy_1d = np.average(Sxy_int,axis=0);
    #mf_px_alt = mom_flux_p_alt(pre_1d,eta_1d,L0,N);
    #mf_vx_alt = mom_flux_v_alt(Sxx_1d,Sxy_1d,eta_1d,L0,N,mu2);
    #
    # compute amplitude
    #
    [eta2,akrms,akpp] = get_amp(eta_int,k_,L0);
    #
    # load the 2d span-averaged binary files
    #
    print("Load field and phase partition")
    fv_2D = return_field(work_dir,'fv',istep_c,N,1,1,pos_c);
    ux_2D = return_field(work_dir,'ux',istep_c,N,1,1,pos_c);
    uy_2D = return_field(work_dir,'uy',istep_c,N,1,1,pos_c);
    uz_2D = return_field(work_dir,'uz',istep_c,N,1,1,pos_c);
    pr_2D = return_field(work_dir,'pr',istep_c,N,1,1,pos_c);
    uv_2D = return_field(work_dir,'uv',istep_c,N,1,1,pos_c);
    du_2D = return_field(work_dir,'du',istep_c,N,1,1,pos_c);
    ke_2D = return_field(work_dir,'ke',istep_c,N,1,1,pos_c);
    vk_2D = return_field(work_dir,'vk',istep_c,N,1,1,pos_c);
    di_2D = return_field(work_dir,'di',istep_c,N,1,1,pos_c);
    dk_2D = return_field(work_dir,'dk',istep_c,N,1,1,pos_c);
    #
    print("From cartesian to wf")
    chosen_mapping = mapping_function_3
    fv_2D_wf, zeta = transform_to_wave_following(fv_2D, eta_1d, x_int, z_int, chosen_mapping, N, L0, k_)
    ux_2D_wf, zeta = transform_to_wave_following(ux_2D, eta_1d, x_int, z_int, chosen_mapping, N, L0, k_)
    uy_2D_wf, zeta = transform_to_wave_following(uy_2D, eta_1d, x_int, z_int, chosen_mapping, N, L0, k_)
    uz_2D_wf, zeta = transform_to_wave_following(uz_2D, eta_1d, x_int, z_int, chosen_mapping, N, L0, k_)
    pr_2D_wf, zeta = transform_to_wave_following(pr_2D, eta_1d, x_int, z_int, chosen_mapping, N, L0, k_)
    uv_2D_wf, zeta = transform_to_wave_following(uv_2D, eta_1d, x_int, z_int, chosen_mapping, N, L0, k_)
    du_2D_wf, zeta = transform_to_wave_following(du_2D, eta_1d, x_int, z_int, chosen_mapping, N, L0, k_)
    ke_2D_wf, zeta = transform_to_wave_following(ke_2D, eta_1d, x_int, z_int, chosen_mapping, N, L0, k_)
    vk_2D_wf, zeta = transform_to_wave_following(vk_2D, eta_1d, x_int, z_int, chosen_mapping, N, L0, k_)
    di_2D_wf, zeta = transform_to_wave_following(di_2D, eta_1d, x_int, z_int, chosen_mapping, N, L0, k_)
    dk_2D_wf, zeta = transform_to_wave_following(dk_2D, eta_1d, x_int, z_int, chosen_mapping, N, L0, k_)
    #
    # Split in air and water phases
    #
    print("Split profiles")
    count_wat = 0;
    count_air = 0;
    for q in range(N):
      #
      if(zeta[q] <= np.average(eta_1d)):
        count_wat +=1;
      else:
        count_air +=1;
    #
    zeta_wat = zeta[0:count_wat];
    zeta_air = zeta[count_wat:N];
    #
    fv_wat_1D_wf, fv_air_1D_wf = split_profile(fv_2D_wf,count_wat,count_air,N,-1) 
    ux_wat_1D_wf, ux_air_1D_wf = split_profile(ux_2D_wf,count_wat,count_air,N,-1) 
    uy_wat_1D_wf, uy_air_1D_wf = split_profile(uy_2D_wf,count_wat,count_air,N,-1) 
    uz_wat_1D_wf, uz_air_1D_wf = split_profile(uz_2D_wf,count_wat,count_air,N,-1) 
    pr_wat_1D_wf, pr_air_1D_wf = split_profile(pr_2D_wf,count_wat,count_air,N,-1) 
    uv_wat_1D_wf, uv_air_1D_wf = split_profile(uv_2D_wf,count_wat,count_air,N,-1) 
    du_wat_1D_wf, du_air_1D_wf = split_profile(du_2D_wf,count_wat,count_air,N,-1) 
    ke_wat_1D_wf, ke_air_1D_wf = split_profile(ke_2D_wf,count_wat,count_air,N,-1) 
    vk_wat_1D_wf, vk_air_1D_wf = split_profile(vk_2D_wf,count_wat,count_air,N,-1) 
    di_wat_1D_wf, di_air_1D_wf = split_profile(di_2D_wf,count_wat,count_air,N,-1) 
    dk_wat_1D_wf, dk_air_1D_wf = split_profile(dk_2D_wf,count_wat,count_air,N,-1) 
    #
    #print("Final print of glo_obs")
    #f = open("glo_obs_post.out","a");
    #f.write("%.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f \n" %(1.0*istep, time, akrms, akpp, eta2, mf_px, mf_py, mf_pz, en_p, en_v));
    #f.flush();
    #f.close();
    #
    #print("Final print of glo_obs alt")
    #f = open("glo_obs_post_alt.out","a");
    #f.write("%.15f %.15f %.15f %.15f \n" %(1.0*istep, time, mf_px_alt, mf_vx_alt));
    #f.flush();
    #f.close();
    #
    plt.contour(x_int,z_int,fv_2D_wf,[0.5])
    plt.contourf(x_int,z_int,ux_2D)
    #
    print("Final print of wf")
    f = open("wave_coord/prof_wf_air_"+istep_c+".out","w");
    for i in range(count_air):
        f.write("%.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f\n" %(zeta_air[i], 
                                                                                        ux_air_1D_wf[i],   # 1 
                                                                                        uy_air_1D_wf[i],   # 2
                                                                                        uz_air_1D_wf[i],   # 3
                                                                                        pr_air_1D_wf[i],   # 4
                                                                                        uv_air_1D_wf[i],   # 5
                                                                                        du_air_1D_wf[i],   # 6
                                                                                        ke_air_1D_wf[i],   # 7
                                                                                        vk_air_1D_wf[i],   # 8
                                                                                        di_air_1D_wf[i],   # 9
                                                                                        dk_air_1D_wf[i])); # 10
    f.flush();
    f.close();
    #
    f = open("wave_coord/prof_wf_wat_"+istep_c+".out","w");
    for i in range(count_wat):
        f.write("%.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f\n" %(zeta_wat[i], 
                                                                                        ux_wat_1D_wf[i],   # 1 
                                                                                        uy_wat_1D_wf[i],   # 2
                                                                                        uz_wat_1D_wf[i],   # 3
                                                                                        pr_wat_1D_wf[i],   # 4
                                                                                        uv_wat_1D_wf[i],   # 5
                                                                                        du_wat_1D_wf[i],   # 6
                                                                                        ke_wat_1D_wf[i],   # 7
                                                                                        vk_wat_1D_wf[i],   # 8
                                                                                        di_wat_1D_wf[i],   # 9
                                                                                        dk_wat_1D_wf[i])); # 10
    f.flush();
    f.close();

