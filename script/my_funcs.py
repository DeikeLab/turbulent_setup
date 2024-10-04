#
# module import
#
import numpy as np
from scipy.interpolate import griddata, interp1d
from scipy.signal import savgol_filter
#
def phase_partion(ux,uy,f,s1_exp,s2_exp):
    ux_air   = ux*(1.0-f)**s1_exp
    ux_water = ux*f**(1/s2_exp)
    uy_air   = uy*(1.0-f)**s1_exp
    uy_water = uy*f**(1/s2_exp)
    return ux_air,ux_water,uy_air,uy_water
#
def interp_2d(xdata,zdata,xtile,ztile,fld):
    #
    fld_int = griddata((xdata.ravel(), zdata.ravel()), fld.ravel(), (xtile, ztile), method='nearest');
    #
    return fld_int;
#
def load_bin(pfile,N,is_2d):
    #
    if is_2d == 1:
      snapshot = np.fromfile(pfile, dtype=np.float64);
      snapshot = snapshot.reshape([N,N]);
      snapshot = np.transpose(snapshot);
    else:
      snapshot = np.fromfile(pfile, dtype=np.float64);
      snapshot = snapshot.reshape([N,N,N]);
      snapshot = np.transpose(snapshot);
    #
    return snapshot;
#
def return_file(common_path,name,istep,N,is_2d):
    #
    if is_2d == 1:
      pfile = common_path+'field/'+name+'_2d_avg_'+istep+'.bin';
      file  = load_bin(pfile,N,is_2d);
    else:
      pfile = common_path+'field/'+name+'_3d_'+istep+'.bin';
      file  = load_bin(pfile,N,is_2d);
    #
    return file;
#
def cart_to_wf(p_2d,eta_1d,N,L0,k_,eta_m0):
    #
    # note p is 2d (in the x, z), eta is 1d (only x)
    #
    eta_1d      = eta_1d;         # the mean has been already subtracted
    p_2d_interp = np.zeros([N,N])
    zplot       = np.zeros([N,N]) # To show in the original cartesian grid z where the interpolating grid z' is
    for i in range(N): # For each x
      #
      z        = np.linspace(-eta_m0,L0-eta_m0,N,endpoint=False);
      f        = interp1d(z,p_2d[i,:],kind='quadratic',fill_value = "extrapolate");
      zeta     = z;
      zplot[i] = zeta + eta_1d[i]*np.exp(-k_*np.abs(zeta));
      p_2d_interp[i] = f(zplot[i]);
      #
    p_1d_interp = np.average(p_2d_interp,axis=-1); # average along the x direction
    #
    return p_2d_interp, p_1d_interp, zplot, zeta
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
    ak = k*( (2/(L0**2))*np.std( (eta[:])**2) )**0.5;
    #
    return ak
