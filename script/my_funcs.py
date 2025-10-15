#
# module import
#
import numpy as np
from scipy.interpolate import griddata, interp1d
from scipy.signal import savgol_filter
#
def interp_2d(xdata,zdata,xtile,ztile,fld):
    #
    fld_int = griddata((xdata.ravel(), zdata.ravel()), fld.ravel(), (xtile, ztile), method='nearest');
    #
    return fld_int;
#
def load_bin(pfile,N,is_2d):
    #
    snapshot = np.fromfile(pfile, dtype=np.float64);
    if is_2d == 1:
      snapshot = snapshot.reshape([N,N]);
    else:
      snapshot = snapshot.reshape([N,N,N]);
    snapshot = np.transpose(snapshot);
    #
    return snapshot;
#
def return_field(common_path,name,istep,N,is_2d,span_avg,pos):
    #
    if is_2d == 1:
      if span_avg == 1:
        pfile = common_path+'field/'+name+'_2d_avg_'+istep+'.bin';
      else: 
        pfile = common_path+'slices/'+name+'_2d_'+pos+'_'+istep+'.bin';
    else:
        pfile = common_path+'field/'+name+'_3d_'+istep+'.bin';
    #
    field = load_bin(pfile,N,is_2d);
    #
    return field;
#
def split_profile(p_2D,count_wat,count_air,N,i_dir):
    #
    p_wat_2D_wf = p_2D[0:count_wat,:];
    p_air_2D_wf = p_2D[count_wat:N,:];
    p_wat_1D_wf = np.average(p_wat_2D_wf,i_dir);
    p_air_1D_wf = np.average(p_air_2D_wf,i_dir);
    #
    if(count_wat+count_air != N):
        print("Warning in split profile");
    #
    return p_wat_1D_wf,p_air_1D_wf;
#
# Define different mapping functions
#
def mapping_function_1(xi, eta_x, mean_interface_position, L0, k_):
    #
    z_mapped_water = xi[xi <= mean_interface_position] * eta_x  # Scale water phase from 0 to eta(x)
    z_mapped_air = eta_x + (xi[xi > mean_interface_position] - mean_interface_position) * (L0 - eta_x) / (L0 - mean_interface_position)
    #
    return np.concatenate((z_mapped_water, z_mapped_air))

def mapping_function_2(xi, eta_x, mean_interface_position, L0, k_):
    #
    z_mapped_water = xi[xi <= mean_interface_position] * eta_x**0.8  # Non-linear scaling for water
    z_mapped_air = eta_x + (xi[xi > mean_interface_position] - mean_interface_position) * (L0 - eta_x)**0.5 / ((L0 - mean_interface_position)**0.5)
    #
    return np.concatenate((z_mapped_water, z_mapped_air))

def mapping_function_3(xi, eta_x, mean_interface_position, L0, k_):
    #
    z_mapped_water = xi[xi <= mean_interface_position] * eta_x  # Scale water phase linearly from 0 to eta(x)
    z_mapped_air = eta_x + (xi[xi > mean_interface_position] - mean_interface_position) * (L0 - eta_x) / (L0 - mean_interface_position)
    #
    return np.concatenate((z_mapped_water, z_mapped_air))

def mapping_function_4(xi, eta_x, mean_interface_position, L0, k_):
    #
    sigma = 5e-3;
    #
    return xi + (eta_x-mean_interface_position)*np.exp(-sigma*k_*np.abs(xi));
#
# Function to perform the field transformation based on the chosen mapping
#
def transform_to_wave_following(p_2D_cart, eta, x, z, mapping_func, N, L0, k_):
    #
    p_2D_cart = np.transpose(p_2D_cart);
    xi = np.linspace(0, L0, N)  # Define xi over the full range
    p_2D_wave_following = np.zeros((N, N))  # Initialize transformed field
    mean_interface_position = np.mean(eta)  # Compute mean interface position
    #
    for i in range(N):
        # Use the eta array value at the current x position
        eta_x = eta[i]
       
        # Apply the chosen mapping function to get z_mapped
        z_mapped = mapping_func(xi, eta_x, mean_interface_position, L0, k_)
       
        # Interpolate the field using the mapped z values
        p_2D_wave_following[i, :] = np.interp(z_mapped, z, p_2D_cart[i, :])
    #
    return np.transpose(p_2D_wave_following), xi
#
def mom_flux_p(pre,nx,ny,nz):
    #
    # compute the mean pressure
    #
    pre_mean = np.average(pre);
    #
    # compute momentum flux due to pressure
    #
    mf_px = - np.average( (pre[:]-pre_mean)*nx[:] );
    mf_py = - np.average( (pre[:]-pre_mean)*ny[:] );
    mf_pz = - np.average( (pre[:]-pre_mean)*nz[:] );
    #
    return mf_px,mf_py,mf_pz
#
def mom_flux_p_alt(pre_1d,eta_1d,L0,N):
    #
    # compute momentum flux due to pressure
    #
    etahat = savgol_filter(eta_1d, 31, 4);
    eps    = np.gradient(etahat)/(L0/N);
    mf_px  = np.average( pre_1d[:]*eps[:] );
    #
    return mf_px
#
def mom_flux_v_alt(Sxx_1d,Sxy_1d,eta_1d,L0,N,mu2):
    #
    # compute momentum flux due to pressure
    #
    etahat = savgol_filter(eta_1d, 31, 4);
    eps    = np.gradient(etahat)/(L0/N);
    mf_vx  = mu2*np.average( Sxy_1d[:] - Sxx_1d[:]*eps[:] );
    #
    return mf_vx
#
def mom_flux_v(Sxx,Sxy,Sxz,Syy,Syz,Szz,nx,ny,nz,mu2):
    #
    # compute momentum flux due to viscous dissipation 
    #
    mf_vx = + mu2*np.average( (Sxx[:]*nx[:] + Sxy[:]*ny[:] + Sxz[:]*nz[:]) );
    mf_vy = + mu2*np.average( (Sxy[:]*nx[:] + Syy[:]*ny[:] + Syz[:]*nz[:]) );
    mf_vz = + mu2*np.average( (Sxz[:]*nx[:] + Syz[:]*ny[:] + Szz[:]*nz[:]) );
    #
    return mf_vx,mf_vy,mf_vz
#
def ene_flux_p(pre,u_x,u_y,u_z,nx,ny,nz):
    #
    # compute the mean pressure and velocity
    #
    pre_mean = np.average(pre);
    u_x_mean = np.average(u_x);
    u_y_mean = np.average(u_y);
    u_z_mean = np.average(u_z);
    #
    # compute energy flux due to pressure
    #
    en_p = - np.average( (pre[:]-pre_mean)*(u_x[:]-u_x_mean)*nx[:] +
                         (pre[:]-pre_mean)*(u_y[:]-u_y_mean)*ny[:] +
                         (pre[:]-pre_mean)*(u_z[:]-u_z_mean)*nz[:] );
    #
    return en_p
#
def ene_flux_v(Sxx,Sxy,Sxz,Syy,Syz,Szz,u_x,u_y,u_z,nx,ny,nz,mu2):
    #
    # compute the mean velocity
    #
    u_x_mean = np.average(u_x);
    u_y_mean = np.average(u_y);
    u_z_mean = np.average(u_z);
    #
    # compute energy flux due to viscous dissipation
    #
    en_v = + mu2*np.average( ( Sxx[:]*nx[:]+Sxy[:]*ny[:]+Sxz[:]*nz[:] )*(u_x[:]-u_x_mean) +
                             ( Sxy[:]*nx[:]+Syy[:]*ny[:]+Syz[:]*nz[:] )*(u_y[:]-u_y_mean) +
                             ( Sxz[:]*nx[:]+Syz[:]*ny[:]+Szz[:]*nz[:] )*(u_z[:]-u_z_mean) );
    #
    return en_v
#
def get_amp(eta,k,L0):
    #
    eta2  = np.sqrt(np.mean(eta[:])**2); # mean square slope 
    akrms = k*( (2/(L0**2))*np.std( (eta[:])**2) )**0.5; # amplitude rms
    akpp  = k*0.5*(eta[:].max()-eta[:].min()); # amplitude peak-to-peak
    #
    return eta2,akrms,akpp
