/** 
This is the wave wind interaction simulation with logarithmic wind profile and adaptive grid. This is a turbulent configuration. */

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"  // reduced gravity
#include "view.h"
#include "tag.h"
#include "lambda2.h"
#include "navier-stokes/perfs.h"
#include "sandbox/profile6.h"   // From Antoon
//#include "sandbox/frac-dist.h"  // From Antoon
//#include "sandbox/LS_reinit.h"  // From Alex

/**
   By fixing the thermophysical properties, 
   the primary parameters are i) the friction Reynolds number $RE_tau$, ii) the Bond, $Bo$, 
   iii) the maximum level of the simulation MAXLEVEL, iv) the initial steepness. 

   Note that we define the values here and then later we pass as arguments in <int main>, together
   with wave speed, RELEASETIME, uemaxRATIO, alterMU and end_sim 
   
   */

double RE_tau;      // Friction Reynolds number
double BO;          // Bond number
int MAXLEVEL;       // max level of the simulation
//double UstarRATIO;  // inverse of wave age
double c_;          // wave-speed
double ak;          // initial steepness
double RELEASETIME; // physical time before which the precursor simulation should be run
double uemaxRATIO;  // Ratio to control the maximum error on the water velocity
double alter_MU;    // Factor to alter the water dynamic viscosity.
double end_sim;     // end of the simulation (physical time unit)
int do_eta_loc;     // output or not eta_loc
int do_profile;     // output or not profiles
int do_fields;      // output or not fields
double tout_eta;    // output frequency of interfacial quantity
double tout_slc;    // output frequency of slices
double tout_pro;    // output frequency of 1d profile
double tout_res;    // output frequency of restart
double tout_mov;    // output frequency of the movie
int prt_res;        // printing resolution

/**
   We define these values: the wave number, fluid depth, wave period, acceleration, wave RE,
   friction velocity Ustar. */

double k_  = 1.0; // we change later 
double h_  = 1.0; // we change later 
double T0_ = 1.0; // we change later
double g_  = 1.0; // we change later
double RE  = 1.0; // we change later
double Ustar = 1.0; // we change later
//double c_  = 1.0; // we change later

/**
   Clean the VoF from bubbles/droplets. */

bool do_remove   = true;   // do it or not
int min_size     = 3;      // size of the bubble/droplet in grid cells
double threshold = 1.0e-4; // VoF and 1.0-VoF threshold
double cirp_th   = 0.20;   // minimum threshold for position
double eta_m0    = 1.0;    // initial mean eta0

/**
   Define some output frequencies. */

double tout_glo_my = 128.0;  // output frequency of global observables
double tout_eta_my = 32.00;  // output frequency of interfacial quantity
double tout_slc_my = 32.00;  // output frequency of slices
double tout_pro_my = 8.00;   // output frequency of 1d profile
double tout_res_my = 4.00;   // output frequency of restart
double tout_mov_my = 12.00;  // output frequency of the movie
double tout_img_my = 6.00;   // output frequency of the images
double tout_rbk_my = 0.50;   // output frequency of the backup restarting files

/**
   For the restarting step. */

int counter_max = 2;
int counter = 0;

/**
   The error on the components of the velocity field used for adaptive
   refinement. */

double uemax  = 0.01;
double femax  = 0.0001;
double uwemax = 0.001;

/**
   The density and viscosity ratios are those of air and water. 
   note: the ratio is air properties over water properties */

double RHO_RATIO = 1.225/1000.;
double MU_RATIO  = 18.31e-6/10.0e-4;

/** 
   We define a water velocity field for the refinement criteria */

vector u_water[]; 

/**
   We need to store the variable forcing field av[]. */

face vector av[];

/**
   Input paramters are passed in from the command line. */

int main(int argc, char *argv[]) {

  /**
     Provide the inputs */

  if (argc > 1)
    RE_tau  = atof (argv[1]);
  if (argc > 2)
    BO = atof(argv[2]);
  if (argc > 3)
    MAXLEVEL = atoi(argv[3]);
  if (argc > 4)
    //UstarRATIO = atof(argv[4]);
    c_ = atof(argv[4]);
  if (argc > 5)
    ak = atof(argv[5]);
  if (argc > 6)
    RELEASETIME = atof(argv[6]);
  if (argc > 7)
    uemaxRATIO = atof(argv[7]);
  if (argc > 8)
    alter_MU = atof(argv[8]);
  if (argc > 9)
    end_sim  = atof(argv[9]);
  if (argc > 10)
    do_eta_loc = atoi(argv[10]); // do_eta_loc = 1 --> true, else --> false
  if (argc > 11)
    do_profile = atoi(argv[11]); // do_profile = 1 --> true, else --> false
  if (argc > 12)
    do_fields  = atoi(argv[12]); // do_fields  = 1 --> true, else --> false
  if (argc > 13)
    tout_eta   = atof(argv[13]); 
  if (argc > 14)
    tout_slc   = atof(argv[14]); 
  if (argc > 15)
    tout_pro   = atof(argv[15]); 
  if (argc > 16)
    tout_res   = atof(argv[16]); 
  if (argc > 17)
    tout_mov   = atof(argv[17]); 
  if (argc > 18)
    prt_res    = atoi(argv[18]); 

  if (argc < 19) {

    fprintf(ferr, "Lack of command line arguments. Check! Need %d more arguments\n", 19-argc);
    return 1;

  }

  fprintf(stderr, "************************\n"), fflush (stderr);
  fprintf(stderr, "Check of input parameters only\n"), fflush (stderr);
  fprintf(stderr, " RE_tau = %8E\n ", RE_tau), fflush (stderr);
  fprintf(stderr, " Bo     = %8E\n ", BO), fflush (stderr);
  fprintf(stderr, " LEVEL  = %d\n ", MAXLEVEL), fflush (stderr);
  //fprintf(stderr, " UstarRATIO = %8E\n ", Ustar/c_), fflush (stderr);
  fprintf(stderr, " Wave speed = %8E\n ", c_), fflush (stderr);
  fprintf(stderr, " a_0k  = %8E\n ", ak), fflush (stderr);
  fprintf(stderr, " RELEASETIME = %8E\n ", RELEASETIME), fflush (stderr);
  fprintf(stderr, " uemaxRatio = %8E\n ", uemaxRATIO), fflush (stderr);
  fprintf(stderr, " alter_MU = %8E\n ", alter_MU), fflush (stderr);
  fprintf(stderr, " tend  = %8E\n ", end_sim), fflush (stderr);
  fprintf(stderr, " do_eta_loc = %d\n ", do_eta_loc), fflush (stderr);
  fprintf(stderr, " do_profile = %d\n ", do_profile), fflush (stderr);
  fprintf(stderr, " do_fields  = %d\n ", do_fields), fflush (stderr);
  fprintf(stderr, " Output frequency for eta_loc  = %8E\n ", tout_eta), fflush (stderr);
  fprintf(stderr, " Output frequency for slices   = %8E\n ", tout_slc), fflush (stderr);
  fprintf(stderr, " Output frequency for profiles = %8E\n ", tout_pro), fflush (stderr);
  fprintf(stderr, " Output frequency for restart  = %8E\n ", tout_res), fflush (stderr);
  fprintf(stderr, " Output frequency for movies   = %8E\n ", tout_mov), fflush (stderr);
  fprintf(stderr, " Printing resolution = %d\n ", prt_res), fflush (stderr);
  fprintf(stderr, "************************\n"), fflush (stderr);

  /**
     The domain is a cubic box centered on the origin and of length
     $L_0=2*\pi$, periodic in the x- and z-directions. */

  L0 = 2*pi;
  h_ = 1; // Water depth set L0/(2*pi) following Wu, Popinet & Deike JFM2022
  k_ = 4; // Four waves per box
  origin (-L0/2., 0, -L0/2.);

  // According to http://basilisk.fr/Basilisk%20C#boundary-conditions
  // for top, u.n = u.y, u.t = u.z, u.r = u.x
  u.r[top]    = neumann(0);
  u.r[bottom] = neumann(0);
  u.n[top]    = dirichlet(0);  
  u.n[bottom] = dirichlet(0);
  u.t[top]    = neumann(0);
  u.t[bottom] = neumann(0);
  periodic (right);
  periodic (front);

  /**
     Here we set the densities and viscosities corresponding to the
     parameters above. Note that these variables are defined in two-phase.h already. 
     To set the dimensional parameters we can use two strategies:
       1. If we want to sweep the space $a_0k$ and $U_\star/c$, we set $Ustar$ to a fixed value, say 0.25, and we modify
          $c$. This strategy allows to keep the same $Re_\tau$ and if we also alter the water viscosity, we can keep fixed Re_wave (strategy of Wu, Popinet and Deike, JFM2022)
       2. More general choice: keep fix $c$ and modify $U_\star/c$, but we cannot keep Re_tau and Re_wave simultaneously fixed (unless to alter the properties ratio)
     */

  //Ustar   = c_*UstarRATIO; // Obsolete
  Ustar   = 0.25; // Pick a fixed value
  rho1    = 1.;
  rho2    = rho1*RHO_RATIO;
  g_      = k_*sq(c_)/( 1.0+(rho1-rho2)/(rho1*BO) ); // tune g_ 
  f.sigma = g_*(rho1-rho2)/(BO*sq(k_));
  G.y     = -g_;
  mu2     = Ustar*rho2*(L0-h_)/RE_tau;
  mu1     = mu2/(MU_RATIO)/alter_MU;
  RE      = rho1*c_*(2*pi/k_)/mu1; // RE now becomes a dependent Non-dim number on c
  T0_     = 2.*pi/sqrt(g_*k_ + f.sigma*k_*k_*k_/rho1);
    
  /*
  rho1    = 1.;
  rho2    = rho1*RHO_RATIO;
  g_      = k_*sq(c_)/( 1.0+(rho1-rho2)/(rho1*BO) ); // tune g_ 
  f.sigma = g_*(rho1-rho2)/(BO*sq(k_));
  c_      = sqrt(g_/k_ + f.sigma*k_/rho0);
  Ustar   = c_*UstarRATIO;
  G.y     = -g_;
  mu2     = Ustar*rho2*(L0-h_)/RE_tau;
  //mu1     = mu2/(MU_RATIO)/alter_MU;
  mu1     = mu2/MU_RATIO;
  RE      = rho1*c_*(2*pi/k_)/mu1; // RE now becomes a dependent Non-dim number on c
  T0_     = 2.*pi/sqrt(g_*k_ + f.sigma*k_*k_*k_/rho1);
  */

  /**
     Give the address of av[] to a[] so that acceleration can be changed */

  a = av;
  init_grid (1 << 7);
  // Refine according to 
  uemax  = uemaxRATIO*Ustar;
  uwemax = 0.001*c_;
  fprintf (stderr, "RELEASETIME = %8E, uemax = %8E \n", RELEASETIME, uemax);
  run();

}

/** 
   Specify the interface shape using a third-order Stokes wave. */

double WaveProfile (double x, double y) {

  double a_ = ak/k_;
  double eta1 = a_*cos(k_*x);
  double alpa = 1./tanh(k_*h_);
  double eta2 = 1./4.*alpa*(3.*sq(alpa) - 1.)*sq(a_)*k_*cos(2.*k_*x);
  double eta3 = -3./8.*(cube(alpa)*alpa - 
			3.*sq(alpa) + 3.)*cube(a_)*sq(k_)*cos(k_*x) + 
    3./64.*(8.*cube(alpa)*cube(alpa) + 
	    (sq(alpa) - 1.)*(sq(alpa) - 1.))*cube(a_)*sq(k_)*cos(3.*k_*x);
  return eta1 + eta2 + eta3 + h_;

}

/**
   Random noise gets killed by adaptive mesh refinement. Therefore, instead of initializing with a mean log profile plus random noise, we initialize with only the mean, and we wait for the instability to naturally develop (because of the shear profile). */

event init (i = 0) {

  /**
     Print relevant simulation parameters */

  fprintf(stderr, "************************\n"), fflush (stderr);
  fprintf(stderr, "A-posteriori check of simulation parameters\n"), fflush (stderr);
  fprintf(stderr, " RE_tau = %8E\n ", Ustar*rho2*(L0-h_)/mu2), fflush (stderr);
  fprintf(stderr, " Bo    = %8E\n ", g_*(rho1-rho2)/(f.sigma*sq(k_))), fflush (stderr);
  fprintf(stderr, " LEVEL = %d\n ", MAXLEVEL), fflush (stderr);
  fprintf(stderr, " wave speed  = %8E\n ", c_), fflush (stderr);
  fprintf(stderr, " a_0k  = %8E\n ", ak), fflush (stderr);
  fprintf(stderr, " RELEASETIME  = %8E\n ", RELEASETIME), fflush (stderr);
  fprintf(stderr, " uemaxRATIO  = %8E\n ", uemaxRATIO), fflush (stderr);
  fprintf(stderr, " alter_MU  = %8E\n ", alter_MU), fflush (stderr);
  fprintf(stderr, " end_sim  = %8E\n ", end_sim), fflush (stderr);
  fprintf(stderr, " do_eta_loc = %d\n ", do_eta_loc), fflush (stderr);
  fprintf(stderr, " do_profile = %d\n ", do_profile), fflush (stderr);
  fprintf(stderr, " do_fields  = %d\n ", do_fields), fflush (stderr);
  fprintf(stderr, " gravity = %8E\n ", g_), fflush (stderr);
  fprintf(stderr, " RE_wave = %8E\n ", rho1*c_*(2*pi/k_)/mu1), fflush (stderr);
  fprintf(stderr, " rho_r = %8E\n ", (rho1/rho2)), fflush (stderr);
  fprintf(stderr, " mu_r  = %8E\n ", (mu1/mu2)), fflush (stderr);
  fprintf(stderr, " UstarRATIO = %8E\n ", Ustar/c_), fflush (stderr);
  fprintf(stderr, " T0 = %8E\n ", T0_), fflush (stderr);
  //fprintf(stderr, " Output frequency for eta_loc  = %8E\n ", tout_eta), fflush (stderr);
  //fprintf(stderr, " Output frequency for slices   = %8E\n ", tout_slc), fflush (stderr);
  //fprintf(stderr, " Output frequency for profiles = %8E\n ", tout_pro), fflush (stderr);
  //fprintf(stderr, " Output frequency for restart  = %8E\n ", tout_res), fflush (stderr);
  //fprintf(stderr, " Output frequency for movies   = %8E\n ", tout_mov), fflush (stderr);
  fprintf(stderr, " Output frequency for eta_loc  = %8E\n ", tout_eta_my), fflush (stderr);
  fprintf(stderr, " Output frequency for slices   = %8E\n ", tout_slc_my), fflush (stderr);
  fprintf(stderr, " Output frequency for profiles = %8E\n ", tout_pro_my), fflush (stderr);
  fprintf(stderr, " Output frequency for restart  = %8E\n ", tout_res_my), fflush (stderr);
  fprintf(stderr, " Output frequency for movies   = %8E\n ", tout_mov_my), fflush (stderr);

  fprintf(stderr, " We_w = %8E\n ", rho1*sq(c_)*(2.0*pi/k_)/f.sigma ), fflush (stderr);
  fprintf(stderr, " Fr_w = %8E\n ", sq(c_)/(g_*2.0*pi/k_) ), fflush (stderr);
  fprintf(stderr, " We_a = %8E\n ", rho2*sq(Ustar)*(L0-h_)/f.sigma ), fflush (stderr);
  fprintf(stderr, " Fr_a = %8E\n ", sq(Ustar)/(g_*(L0-h_)) ), fflush (stderr);

  fprintf(stderr, " Printing resolution = %d\n ", prt_res), fflush (stderr);
  fprintf(stderr, "************************\n"), fflush (stderr);

  /**
     Create directories for better organization of the output files */

  fprintf(stderr, "Create directories (if needed) \n"), fflush (stderr);

  char comm[80];
  sprintf (comm, "mkdir -p profiles");
  system(comm);
  sprintf (comm, "mkdir -p eta");
  system(comm);
  sprintf (comm, "mkdir -p eta/eta_loc");
  system(comm);
  sprintf (comm, "mkdir -p field");
  system(comm);
  sprintf (comm, "mkdir -p images");
  system(comm);
  sprintf (comm, "mkdir -p restart_bk");
  system(comm);
  sprintf (comm, "mkdir -p tagging_f1");
  system(comm);
  sprintf (comm, "mkdir -p tagging_f2");
  system(comm);

  /**
     Restart from the latest dump files or initialize a new simulation */

  if (restore ("restart.bin")) {
    fprintf(stderr, "Simulation restarts from a dumped file\n"), fflush (stderr);
  }
  else {
    fprintf(stderr, "Initialization of all the variables\n"), fflush (stderr);
    double ytau = (mu2/rho2)/Ustar;
    do {
      fraction (f, WaveProfile(x,z)-y);
      foreach() {
	// Initialize with a fairly accurate profile
	if ((y-WaveProfile(x,z)) > 0.05) {
	  u.x[] = (1-f[])*(log((y-WaveProfile(x,z))/ytau)*Ustar/0.41);
	}
	else {
	  u.x[] = 0.;
	}
	u.y[] = 0.;
	u.z[] = 0.;
      }
      boundary ((scalar *){u});
    }
#if TREE
    // No need for adaptation when starting 
    while (0);
#else
    while (0);
#endif
  }

}

/**
## Outputs

   We are interested in the viscous dissipation rate in both water and air. */

int dissipation_rate (double* rates) {
	
  double rateWater = 0.0;
  double rateAir = 0.0;
  foreach (reduction (+:rateWater) reduction (+:rateAir)) {
    double dudx = (u.x[1]     - u.x[-1]    )/(2.*Delta);
    double dudy = (u.x[0,1]   - u.x[0,-1]  )/(2.*Delta);
    double dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
    double dvdx = (u.y[1]     - u.y[-1]    )/(2.*Delta);
    double dvdy = (u.y[0,1]   - u.y[0,-1]  )/(2.*Delta);
    double dvdz = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
    double dwdx = (u.z[1]     - u.z[-1]    )/(2.*Delta);
    double dwdy = (u.z[0,1]   - u.z[0,-1]  )/(2.*Delta);
    double dwdz = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);
    double SDeformxx = dudx;
    double SDeformxy = 0.5*(dudy + dvdx);
    double SDeformxz = 0.5*(dudz + dwdx);
    double SDeformyx = SDeformxy;
    double SDeformyy = dvdy;
    double SDeformyz = 0.5*(dvdz + dwdy);
    double SDeformzx = SDeformxz;
    double SDeformzy = SDeformyz;
    double SDeformzz = dwdz; 
    double sqterm = 2.*dv()*(sq(SDeformxx) + sq(SDeformxy) + sq(SDeformxz) +
			     sq(SDeformyx) + sq(SDeformyy) + sq(SDeformyz) +
			     sq(SDeformzx) + sq(SDeformzy) + sq(SDeformzz)) ;
    rateWater += mu1/rho[]*f[]*sqterm; //water
    rateAir   += mu2/rho[]*(1. - f[])*sqterm; //air
  }
  rates[0] = rateWater;
  rates[1] = rateAir;
  return 0;

}

/**
   We log the physical time, the number of time-steps and dt of the simulation. */

event log_simulation (i += 10) {
  
  if (pid() == 0) { 

    fflush(stderr);
    char name_1[80];
    sprintf (name_1,"log_simulation.out");
    FILE * log_sim = fopen(name_1,"a");
    fprintf (log_sim, "%8E %8E %8E\n", t, 1.0*i, dt);
    fclose(log_sim);

  }

}

/**
  Keep the wave profile and set the wave velocity to 0, until the RELEASETIME. */
  
event set_wave(i=0; i++; t<RELEASETIME) {

  fraction (f, WaveProfile(x,z)-y);
  foreach() {
    u.x[] = (1.0 - f[])*u.x[];
    u.y[] = (1.0 - f[])*u.y[];
    u.z[] = (1.0 - f[])*u.z[];
  }
  boundary ((scalar *){u});

}

/**
   Release the wave at RELEASETIME. We don't do any adaptation at this step. 
   */

// TO-DO: maybe include g correction by Bond number. High Bond should be fine.

double u_x (double x, double y) {

  double alpa = 1./tanh(k_*h_);
  double a_ = ak/k_;
  double sgma = sqrt(g_*k_*tanh(k_*h_)*
		     (1. + k_*k_*a_*a_*(9./8.*(sq(alpa) - 1.)*
					(sq(alpa) - 1.) + sq(alpa))));
  double A_ = a_*g_/sgma;
  return A_*cosh(k_*(y + h_))/cosh(k_*h_)*k_*cos(k_*x) +
    3.*ak*A_/(8.*alpa)*(sq(alpa) - 1.)*(sq(alpa) - 1.)*
    cosh(2.0*k_*(y + h_))*2.*k_*cos(2.0*k_*x)/cosh(2.0*k_*h_) +
    1./64.*(sq(alpa) - 1.)*(sq(alpa) + 3.)*
    (9.*sq(alpa) - 13.)*
    cosh(3.*k_*(y + h_))/cosh(3.*k_*h_)*ak*ak*A_*3.*k_*cos(3.*k_*x);

}

double u_y (double x, double y) {

  double alpa = 1./tanh(k_*h_);
  double a_ = ak/k_;
  double sgma = sqrt(g_*k_*tanh(k_*h_)*
		     (1. + k_*k_*a_*a_*(9./8.*(sq(alpa) - 1.)*
					(sq(alpa) - 1.) + sq(alpa))));
  double A_ = a_*g_/sgma;
  return A_*k_*sinh(k_*(y + h_))/cosh(k_*h_)*sin(k_*x) +
    3.*ak*A_/(8.*alpa)*(sq(alpa) - 1.)*(sq(alpa) - 1.)*
    2.*k_*sinh(2.0*k_*(y + h_))*sin(2.0*k_*x)/cosh(2.0*k_*h_) +
    1./64.*(sq(alpa) - 1.)*(sq(alpa) + 3.)*
    (9.*sq(alpa) - 13.)*
    3.*k_*sinh(3.*k_*(y + h_))/cosh(3.*k_*h_)*ak*ak*A_*sin(3.*k_*x);

}

event start(t = RELEASETIME) {
  // A slightly changed version of the test/stokes.h wave function as y = 0 at the bottom now so y+h -> y
  
  fprintf(stderr, "t = RELEASETIME, we initialize the liquid velocity\n"), fflush (stderr);
  fraction (f, WaveProfile(x,z)-y);
  foreach () {
    u.x[] += u_x(x, y-h_)*f[];
    u.y[] += u_y(x, y-h_)*f[];
  }
  boundary ((scalar *){u});

}

/**
   Forcing term to sustain the flow equivalent to the pressure gradient in x direction. */

event acceleration (i++) {

  double ampl = sq(Ustar)/(L0-h_);
  foreach_face(x)
    av.x[] += ampl*(1.-f[]);

  double u_air_mean = 0.;
  double u_wat_mean = 0.;
  double voll_mean  = 0.;
  double volg_mean  = 0.;
  foreach(reduction(+:u_air_mean) reduction(+:u_wat_mean)
          reduction(+:volg_mean ) reduction(+:voll_mean )) {
    u_air_mean += u.x[]*(1.0-f[])*dv();
    u_wat_mean += u.x[]*(0.0+f[])*dv();
    volg_mean  += (1.0-f[])*dv();
    voll_mean  += (0.0+f[])*dv();
  }

  if (pid() == 0) { 

    fflush(stderr);
    char name_1[80];
    sprintf (name_1,"log_forcing.out");
    FILE * log_sim = fopen(name_1,"a");
    fprintf (log_sim, "%8E %8E %8E %8E %8E %8E\n", t, 1.0*i, volg_mean, u_air_mean/volg_mean,
		                                             voll_mean, u_wat_mean/voll_mean);
    fclose(log_sim);

  }

}

double my_interpolation(scalar s, double xp, double yp, double zp) {

  //boundary({s});

#if dimension == 2
  Point point = locate(xp,yp);
  if(point.level < 0) {
    return nodata;
  }
  else {
    double xpo = (xp - x)/Delta - s.d.x/2.0;
    double ypo = (yp - y)/Delta - s.d.y/2.0;
      
    int i = sign(xpo), j = sign(ypo);
    xpo = fabs(xpo), ypo = fabs(ypo);
    
    return ((s[]*(1. - xpo) + s[i]*xpo)*(1. - ypo) +
	    (s[0,j]*(1. - xpo) + s[i,j]*xpo)*ypo);
  }
#else
  Point point = locate(xp,yp,zp);
  if(point.level < 0) {
    return nodata;
  }
  else {
    double xpo = (xp - x)/Delta - s.d.x/2.0;
    double ypo = (yp - y)/Delta - s.d.y/2.0;
    double zpo = (zp - z)/Delta - s.d.z/2.0;
      
    int i = sign(xpo), j = sign(ypo), k = sign(zpo);
    xpo = fabs(xpo), ypo = fabs(ypo), zpo = fabs(zpo);
    
    return (((s[]*(1. - xpo) + s[i]*xpo)*(1. - ypo) + 
	    (s[0,j]*(1. - xpo) + s[i,j]*xpo)*ypo)*(1. - zpo) +
	    ((s[0,0,k]*(1. - xpo) + s[i,0,k]*xpo)*(1. - ypo) + 
	    (s[0,j,k]*(1. - xpo) + s[i,j,k]*xpo)*ypo)*zpo); 
  }
#endif

}

/** 
   Output video and field in uncompressed format 
   (we can compress later using convert <FILE>.ppm to <FILE>.mpg and open with mplayer) */

# define POPEN(name, mode) fopen (name ".ppm", mode)

//event movies (t += T0_/tout_mov) {
event movies (t += T0_/tout_mov_my) {
//event movies (i += 2) {

  //fprintf(stderr, "I output the movies every t=%8E\n", T0_/tout_mov); fflush (stderr);
  fprintf(stderr, "I output the movies every t=%8E\n", T0_/tout_mov_my); fflush (stderr);
  
  scalar l2[];
  lambda2 (u, l2);
  scalar omega[];
  vorticity (u, omega);

  char stg[80];
 
  clear();
  view (fov = 27.5,  
        theta = pi/4.0, phi = pi/50.0, psi = 0.0, ty = -0.5,
	//width = 1300, height = 1300, bg = {1,1,1}, samples = 4);
	width = 900, height = 900, bg = {1,1,1}, samples = 4);
  box(false, lc = {1,1,1}, lw = 0.1);
  draw_vof ("f", color = "u.x");
  squares ("u.x", linear = true, n = {0,0,1}, alpha = -L0/2.0);
  squares ("omega", linear = true, n = {1,0,0}, alpha = -L0/2.0);
  cells (n = {1,0,0}, alpha = -L0/2.0);
  sprintf (stg, "t = %0.3f", 2.0*pi*(t-RELEASETIME)/T0_);
  draw_string (stg, size = 30); 
  {
    static FILE * fp = POPEN ("3D", "a");
    save (fp = fp);
    fflush (fp);
  }

  clear();
  view (fov = 40, camera = "iso", ty = -0.25,
        //width = 1000, height = 1000, bg = {1,1,1}, samples = 4);
        width = 900, height = 900, bg = {1,1,1}, samples = 4);
  box(false, lc = {1,1,1}, lw = 0.1);
  draw_vof ("f", color = "u.x");
  squares ("u.y", linear = true, n = {0,0,1}, alpha = -L0/2.0);
  squares ("u.x", linear = true, n = {1,0,0}, alpha = -L0/2.0);
  cells (n = {1,0,0}, alpha = -L0/2.0);
  isosurface ("l2", -1);
  sprintf (stg, "t = %0.3f", 2.0*pi*(t-RELEASETIME)/T0_);
  draw_string (stg, size = 30); 
  {
    static FILE * fp = POPEN ("vortex", "a");
    save (fp = fp);
    fflush (fp);
  }

  clear();
  view (fov = 22.5, camera = "front", ty = -0.50,
        //theta = 0.0, phi = -pi/240.0, psi = 0.0, ty = -0.5,
        //width = 1300, height = 1300, bg = {1,1,1}, samples = 4);
        width = 900, height = 900, bg = {1,1,1}, samples = 4);
  box(false, lc = {1,1,1}, lw = 0.1);
  draw_vof ("f", color = "omega");
  squares ("u.x", linear = true, n = {0,0,1}, alpha = -L0/2.0);
  sprintf (stg, "t = %0.3f", 2.0*pi*(t-RELEASETIME)/T0_);
  draw_string (stg, size = 30, pos = 3); 
  {
    static FILE * fp = POPEN ("int", "a");
    save (fp = fp);
    fflush (fp);
  }

}

/*
//event images (t += T0_/tout_img) {
event images (t += T0_/tout_img_my) {
//event images (i += 2) {

  //fprintf(stderr, "I output the movies every t=%8E\n", T0_/tout_mov); fflush (stderr);
  fprintf(stderr, "I output the images every t=%8E\n", T0_/tout_mov_my); fflush (stderr);
  
  scalar l2[];
  lambda2 (u, l2);
  scalar omega[];
  vorticity (u, omega);
  
  char stg[80];
  char file[99];
 
  clear();
  view (fov = 27.5,  
        theta = pi/4.0, phi = pi/50.0, psi = 0.0, ty = -0.5,
	width = 1300, height = 1300, bg = {1,1,1}, samples = 4);
  box(false, lc = {1,1,1}, lw = 0.1);
  draw_vof ("f", color = "u.x");
  squares ("u.x", linear = true, n = {0,0,1}, alpha = -L0/2.0);
  squares ("omega", linear = true, n = {1,0,0}, alpha = -L0/2.0);
  cells (n = {1,0,0}, alpha = -L0/2.0);
  sprintf (stg, "t = %0.3f", 2.0*pi*(t-RELEASETIME)/T0_);
  draw_string (stg, size = 30); 
  sprintf (file, "./images/3D_%09d.ppm", i);
  save(file);

  clear();
  view (fov = 40, camera = "iso", ty = -0.25,
        width = 1000, height = 1000, bg = {1,1,1}, samples = 4);
  box(false, lc = {1,1,1}, lw = 0.1);
  draw_vof ("f", color = "u.x");
  squares ("u.y", linear = true, n = {0,0,1}, alpha = -L0/2.0);
  squares ("u.x", linear = true, n = {1,0,0}, alpha = -L0/2.0);
  cells (n = {1,0,0}, alpha = -L0/2.0);
  isosurface ("l2", -1);
  sprintf (stg, "t = %0.3f", 2.0*pi*(t-RELEASETIME)/T0_);
  draw_string (stg, size = 30); 
  sprintf (file, "./images/vortex_%09d.ppm", i);
  save(file);

  clear();
  view (fov = 22.5, camera = "front", ty = -0.50,
        //theta = 0.0, phi = -pi/240.0, psi = 0.0, ty = -0.5,
        width = 1300, height = 1300, bg = {1,1,1}, samples = 4);
  box(false, lc = {1,1,1}, lw = 0.1);
  draw_vof ("f", color = "omega");
  squares ("u.x", linear = true, n = {0,0,1}, alpha = -L0/2.0);
  sprintf (stg, "t = %0.3f", 2.0*pi*(t-RELEASETIME)/T0_);
  draw_string (stg, size = 30, pos = 3); 
  sprintf (file, "./images/int_pic_%09d.ppm", i);
  save(file);

}
*/


void profile_output (int istep) {

  /**
     Compute the vorticity. */

  scalar omega[];
  vorticity (u, omega);

  char file[99];
  sprintf (file, "./profiles/prof_%09d.out", istep);
  scalar uxuy[],uxux[],uyuy[],uzuz[];
  foreach () {
    uxuy[] = u.x[]*u.y[];
    uxux[] = u.x[]*u.x[];
    uyuy[] = u.y[]*u.y[];
    uzuz[] = u.z[]*u.z[];
  }

  /*
  scalar phi[], phi0[];
  foreach() {
    phi[]  = y;
  }
  boundary({phi,phi0});

  int nit = LS_reinit(phi);
  */
  scalar phi[];
  foreach_vertex() {
    phi[] = y;
  }

  profiles ({u.x, u.y, u.z, omega, uxuy, uxux, uyuy, uzuz}, phi, rf = 1.0, fname = file, n = 1 << MAXLEVEL);

}

void cmpt_stress (int istep) {

  // Name of the file
  char file[99];
  sprintf (file, "./profiles/budget_%09d.out", istep);
  
  // We want to compute the streamwise stress budget 
  // using a Favre average. The budget is computed
  // in the overall domain and per phase
  
  vector u_til[] , rhou[];
  scalar uxuy[]  , uxux[], uyuy[], dudy[];
  scalar int_ft[], rhouv[];
#if dimension > 2
  scalar uzuz[], uzuz_w[], uzuz_a[];
#endif

  // Compute phi_tmp = sigma*kappa
  scalar phi_tmp[];
  curvature (f, phi_tmp, f.sigma, add = false);

  foreach () {

    // For the global budget
    foreach_dimension() {
      u_til.x[] = rho[]*u.x[];
      rhou.x[]  = rho[]*u.x[];
    }
    rhouv[] = rho[]*u.x[]*u.y[];
    uxuy[]  = u.x[]*u.y[];
    uxux[]  = u.x[]*u.x[];
    uyuy[]  = u.y[]*u.y[];
    double mu_tmp = mu1*f[]+mu2*(1.0-f[]);
    dudy[]  = mu_tmp*(u.x[0,1]-u.x[0,-1])/(2.*Delta);
#if dimension > 2
    uzuz[]  = u.z[]*u.z[];
#endif

    // Interface
    if(interfacial(point,f) && phi_tmp[] != nodata) {
      int_ft[] = phi_tmp[]*(f[1]-f[-1])/(2.*Delta);
    }
    else {
      int_ft[] = 0;
    }

  }

  // vertical position
  vertex scalar pos_v[];
  foreach_vertex () {
    pos_v[] = y;
  }
  
  /*
#if dimension > 2
  scalar * list_fld = {u_til.x,u_til.y,u_til.z,
		       rhou.x,rhou.y,rhou.z,
		       rhouv,
		       u.x,u.y,u.z,
	               uxux,uyuy,uzuz,uxuy,
		       dudy,
		       f,int_ft};
#else
  scalar * list_fld = {u_til.x,u_til.y,
		       rhou.x,rhou.y,
		       rhouv,
		       u.x,u.y,
	               uxux,uyuy,uxuy,
		       dudz,
		       f,int_ft};
#endif
  */
  profiles ({u_til.x,u_til.y,u_til.z,
             rhou.x,rhou.y,rhou.z,rhouv,u.x,u.y,u.z,
             uxux,uyuy,uzuz,uxuy,
             dudy,
             f,int_ft}, pos_v, rf = 1.0, fname = file, n = 1 << MAXLEVEL);

}

//event out_pro (t += T0_/tout_pro) {
event out_pro (t += T0_/tout_pro_my) {
//event out_pro (i += 2) {

  if (do_profile == 1) {

    //fprintf(stderr, "I output the vertical profiles file every t=%8E\n", T0_/tout_pro);
    fprintf(stderr, "I output the vertical profiles file every t=%8E\n", T0_/tout_pro_my);
    profile_output (i);
    //cmpt_stress (i);
    fflush (stderr);
   
    if (pid() == 0) {
    
      char name_1[80];
      sprintf (name_1,"./profiles/log_pro.out");
      FILE * log_sim = fopen(name_1,"a");
      fprintf (log_sim, "%8E %8E\n", t, 1.0*i);
      fclose(log_sim);
    
    }

  }

}

void * matrix_new_3d (int nx, int ny, int nz, size_t size) {
  void ** m = qmalloc(nx, void *);  //Define a pointer that points to every x coordinate.
  char * a  = qmalloc(nx*ny*nz*size, char);
  for (int i=0; i<nx; i++) {
    m[i] = a+i*ny*nz*size;
  }
  return m;
}

void output_3d(char * fname, scalar s, int maxlevel, bool do_linear, bool print_bin) {

  FILE * fpver = fopen (fname,"w");
  int nn = (1<<maxlevel);
  double stp = 0.999999*(L0+X0-X0)/(double)nn; // to avoid interpolated point coincides with fine grid boundary
  double ** field = (double **) matrix_new_3d (nn, nn, nn, sizeof(double));
  for (int i = 0; i < nn; i++) {
    double xp = stp*i + X0 + stp/2.;
    for (int j = 0; j < nn; j++) {
      double yp = stp*j + Y0 + stp/2.;
      for (int k = 0; k < nn; k++) {
        double zp = stp*k + Z0 + stp/2.;
        if(do_linear) {
          field[i][j*nn+k] = interpolate(s,xp,yp,zp);
	}
        else {
          Point point = locate (xp,yp,zp);
          field[i][j*nn+k] = point.level >= 0 ? s[] : nodata;
	}
      }
    }
  }
  if (pid() == 0) { // master
#if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], cube(nn), MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
#endif
    if(print_bin) { // print in binary format
      for (int i = 0; i < nn; i++) {
        for (int j = 0; j < nn; j++) {
          for (int k = 0; k < nn; k++) {
            fwrite ( &field[i][j*nn+k], sizeof(double), 1, fpver );
          }
        }
      }
    }
    else { // print in ascii format
      for (int i = 0; i < nn; i++) {
        for (int j = 0; j < nn; j++) {
          for(int k = 0; k < nn; k++) {
            fprintf(fpver, "%8E", field[i][j*nn+k]);
            fputc('\n', fpver);  // not double quotation""
          }
        }
      }
    }
    fflush(fpver);
  }
#if _MPI
  else // slave
    MPI_Reduce (field[0], NULL, cube(nn), MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
#endif
  matrix_free (field);
  fclose (fpver); // we close at the end
}

/*
void output_2d_span_avg(char * fname, scalar s, int maxlevel, bool do_linear, bool print_bin) {

  // spanwise averaging
  int nn = (1<<maxlevel);
  double stp = 0.999999*(L0+X0-X0)/(double)nn; // to avoid interpolated point coincides with fine grid boundary
  double ** field = (double **) matrix_new (nn, nn, sizeof(double)); // spanwise averaged
  for (int i = 0; i < nn; i++) {
    double xp = stp*i + X0 + stp/2.;
    for (int j = 0; j < nn; j++) {
      double yp = stp*j + Y0 + stp/2.;
      double add_z = 0;
      field[i][j]  = nodata;
      for (int k = 0; k < nn; k++) {
        double zp = stp*k + Z0 + stp/2.;
        if(do_linear) {
          double avg = interpolate(s,xp,yp,zp);
          add_z += avg;
	}
        else {
          Point point = locate (xp,yp,zp);
          double avg = point.level >= 0 ? s[] : nodata;
          add_z += avg;
	}
      }
      field[i][j] = add_z/(double)nn;
    }
  }

  // parallel reduction and print in a file
  FILE * fpver = fopen (fname,"w");
  if (pid() == 0) { // master
#if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], sq(nn), MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
#endif
    if(print_bin) { // print in binary format
      for (int i = 0; i < nn; i++) {
        for (int j = 0; j < nn; j++) {
	  fwrite ( &field[i][j], sizeof(double), 1, fpver );
        }
      }
    }
    else { // print in ascii format
      for (int i = 0; i < nn; i++) {
        for (int j = 0; j < nn; j++) {
          fprintf(fpver, "%8E", field[i][j]);
          fputc('\n', fpver);  // not double quotation""
        }
      }
    }
    fflush(fpver);
  }
#if _MPI
  else // slave
    MPI_Reduce (field[0], NULL, sq(nn), MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
#endif
  matrix_free (field);
  fclose (fpver); // we close at the end

}
*/

void output_2d_span_avg(char * fname, scalar s, int maxlevel, bool do_linear, bool print_bin) {

  int nn = (1<<maxlevel);
  double stp = 0.999999*(L0+X0-X0)/(double)nn; // to avoid interpolated point coincides with fine grid boundary

  // first create a 3D field
  double ** field = (double **) matrix_new_3d (nn, nn, nn, sizeof(double));
  for (int i = 0; i < nn; i++) {
    double xp = stp*i + X0 + stp/2.;
    for (int j = 0; j < nn; j++) {
      double yp = stp*j + Y0 + stp/2.;
      for (int k = 0; k < nn; k++) {
        double zp = stp*k + Z0 + stp/2.;
        if(do_linear) {
          field[i][j*nn+k] = interpolate(s,xp,yp,zp);
	}
        else {
          Point point = locate (xp,yp,zp);
          field[i][j*nn+k] = point.level >= 0 ? s[] : nodata;
	}
      }
    }
  }

  FILE * fpver = fopen (fname,"w");
  if (pid() == 0) { // master

    // reduce it to pid()
#if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], cube(nn), MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
#endif

    // spanwise averaging
    double ** field_avg_z = (double **) matrix_new (nn, nn, sizeof(double));
    for (int i = 0; i < nn; i++) {
      for (int j = 0; j < nn; j++) {
	field_avg_z[i][j] = 0.0;
        for (int k = 0; k < nn; k++) {
          field_avg_z[i][j] += field[i][j*nn+k]/(double)nn;
        }
      }
    }

    // print
    if(print_bin) { // print in binary format
      for (int i = 0; i < nn; i++) {
        for (int j = 0; j < nn; j++) {
	  fwrite ( &field_avg_z[i][j], sizeof(double), 1, fpver );
        }
      }
    }
    else { // print in ascii format
      for (int i = 0; i < nn; i++) {
        for (int j = 0; j < nn; j++) {
          fprintf(fpver, "%8E", field_avg_z[i][j]);
          fputc('\n', fpver);  // not double quotation""
        }
      }
    }
    fflush(fpver);
    matrix_free (field_avg_z);

  }
#if _MPI
  else // slave
    MPI_Reduce (field[0], NULL, cube(nn), MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
#endif

  matrix_free (field);
  fclose (fpver); // we close at the end

}


/*
//event slice_stat (t += T0_/tout_slc) {
event slice_stat (t += T0_/tout_slc_my) {
//event slice_stat (i += 2) {

  if (do_fields == 1) {

    //fprintf(stderr, "I output the slices every t=%8E\n", T0_/tout_slc), fflush (stderr);
    fprintf(stderr, "I output the slices every t=%8E\n", T0_/tout_slc_my), fflush (stderr);
    
    // Compute the vorticity and the pressure 
    
    scalar omega[];
    vorticity (u, omega);
    boundary({omega});
    
    char filename[100];
    int res        = prt_res;
    //int Nslice    = 0.5*pow(2.0,1.0*res);
    //int Nslice    = 1.0;
    bool do_linear = true;
    bool print_bin = true;
    //double zslice = -L0/2+L0/2./Nslice;
    fprintf(stderr, "I print 3d\n"), fflush (stderr);

    sprintf (filename, "./field/ux_3d_%09d.bin", i); // x-vel
    output_3d (filename,u.x  ,res, do_linear, print_bin);
    sprintf (filename, "./field/uy_3d_%09d.bin", i); // y-vel
    output_3d (filename,u.y  ,res, do_linear, print_bin);
    //sprintf (filename, "./field/uz_3d_%09d.bin", i); // z-vel
    //output_3d (filename,u.z  ,res, do_linear, print_bin);
    //sprintf (filename, "./field/om_3d_%09d.bin", i); // vorticity
    //output_3d (filename,omega,res, do_linear, print_bin);
    sprintf (filename, "./field/pa_3d_%09d.bin", i); // pressure
    output_3d (filename,p    ,res, do_linear, print_bin);
    sprintf (filename, "./field/fv_3d_%09d.bin", i); // VoF
    output_3d (filename,f    ,res, do_linear, print_bin);
    
    //sprintf (filename, "./field/fv_3d_%09d.out", i); // VoF
    //FILE * fpver = fopen (filename,"w");
    //output_field_3d({f},fpver,res,do_linear);

    if (pid() == 0) {
    
      char name_1[80];
      sprintf (name_1,"./field/log_field.out");
      FILE * log_sim = fopen(name_1,"a");
      fprintf (log_sim, "%8E %8E\n", t, 1.0*i);
      fclose(log_sim);
    
    }

  }

}
*/

//event slice_stat (t += T0_/tout_slc) {
event slice_stat (t += T0_/tout_slc_my) {
//event slice_stat (i += 2) {

  if (do_fields == 1) {

    //fprintf(stderr, "I output the slices every t=%8E\n", T0_/tout_slc), fflush (stderr);
    fprintf(stderr, "I output the slices every t=%8E\n", T0_/tout_slc_my), fflush (stderr);
    
    // Compute some quantities
    scalar diss[], rh_uv[], mdudy[];
    foreach () {

      double dudx = (u.x[1]     - u.x[-1]    )/(2.*Delta);
      double dudy = (u.x[0,1]   - u.x[0,-1]  )/(2.*Delta);
      double dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
      double dvdx = (u.y[1]     - u.y[-1]    )/(2.*Delta);
      double dvdy = (u.y[0,1]   - u.y[0,-1]  )/(2.*Delta);
      double dvdz = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
      double dwdx = (u.z[1]     - u.z[-1]    )/(2.*Delta);
      double dwdy = (u.z[0,1]   - u.z[0,-1]  )/(2.*Delta);
      double dwdz = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);
      double SDeformxx = dudx;
      double SDeformxy = 0.5*(dudy + dvdx);
      double SDeformxz = 0.5*(dudz + dwdx);
      double SDeformyx = SDeformxy;
      double SDeformyy = dvdy;
      double SDeformyz = 0.5*(dvdz + dwdy);
      double SDeformzx = SDeformxz;
      double SDeformzy = SDeformyz;
      double SDeformzz = dwdz; 
      double sqterm = 2.*dv()*( sq(SDeformxx) + sq(SDeformxy) + sq(SDeformxz) +
          		        sq(SDeformyx) + sq(SDeformyy) + sq(SDeformyz) +
          		        sq(SDeformzx) + sq(SDeformzy) + sq(SDeformzz) );
      double mu_loc = mu1*f[]+mu2*(1.0-f[]);
      diss[]  = mu_loc/rho[]*sqterm; 

      rh_uv[] = rho[]*u.x[]*u.y[];
      mdudy[] = mu_loc*dudy; 

    }
    boundary({diss,rh_uv,mdudy});
    
    char filename[100];
    int res        = prt_res;
    bool do_linear = true;
    bool print_bin = true;

    sprintf (filename, "./field/ux_2d_avg_%09d.bin", i); // x-vel
    output_2d_span_avg (filename,u.x  ,res, do_linear, print_bin);
    sprintf (filename, "./field/uy_2d_avg_%09d.bin", i); // y-vel
    output_2d_span_avg (filename,u.y  ,res, do_linear, print_bin);
    sprintf (filename, "./field/uz_2d_avg_%09d.bin", i); // z-vel
    output_2d_span_avg (filename,u.z  ,res, do_linear, print_bin);
    sprintf (filename, "./field/fv_2d_avg_%09d.bin", i); // VoF
    output_2d_span_avg (filename,f    ,res, do_linear, print_bin);
    sprintf (filename, "./field/pr_2d_avg_%09d.bin", i); // pressure
    output_2d_span_avg (filename,p    ,res, do_linear, print_bin);
    sprintf (filename, "./field/di_2d_avg_%09d.bin", i); // dissipation
    output_2d_span_avg (filename,diss ,res, do_linear, print_bin);
    sprintf (filename, "./field/ru_2d_avg_%09d.bin", i); // rh_uv
    output_2d_span_avg (filename,rh_uv,res, do_linear, print_bin);
    sprintf (filename, "./field/md_2d_avg_%09d.bin", i); // mdudy
    output_2d_span_avg (filename,mdudy,res, do_linear, print_bin);

    /*
    sprintf (filename, "./field/ux_3d_%09d.bin", i); // x-vel
    output_3d (filename,u.x  ,res, do_linear, print_bin);
    sprintf (filename, "./field/uy_3d_%09d.bin", i); // y-vel
    output_3d (filename,u.y  ,res, do_linear, print_bin);
    sprintf (filename, "./field/uz_3d_%09d.bin", i); // z-vel
    output_3d (filename,u.z  ,res, do_linear, print_bin);
    sprintf (filename, "./field/fv_3d_%09d.bin", i); // VoF
    output_3d (filename,f    ,res, do_linear, print_bin);
    sprintf (filename, "./field/pr_3d_%09d.bin", i); // pressure
    output_3d (filename,p    ,res, do_linear, print_bin);
    sprintf (filename, "./field/om_3d_%09d.bin", i); // omega
    output_3d (filename,omega,res, do_linear, print_bin);
    sprintf (filename, "./field/di_3d_%09d.bin", i); // dissipation
    output_3d (filename,diss ,res, do_linear, print_bin);
    */

    if (pid() == 0) {
    
      char name_1[80];
      sprintf (name_1,"./field/log_field.out");
      FILE * log_sim = fopen(name_1,"a");
      fprintf (log_sim, "%8E %8E\n", t, 1.0*i);
      fclose(log_sim);
    
    }

  }
  //return 1;

}

/**
   We want to compute some global observables - second way */

void output_global_obs_1 (char * fname, int istep, scalar c, double stp_eta) {

  // Preliminary calculations
 
  // --> position
  scalar pos[];
#if dimension > 2
  coord G = {0.,1.,0.}, Z = {0.,0.,0.};
#else
  coord G = {0.,1.}, Z = {0.,0.};
#endif
  position (c, pos, G, Z);
  boundary ({pos});

  // --> my diagnosis of eta_my, area_my, amp_my
  double area_my = 0; double eta_my = 0;
  foreach(reduction(+:area_my), reduction(+:eta_my)) {
    if (interfacial (point, c)) {
      //if (point.level == MAXLEVEL) {
      if (point.level == MAXLEVEL && abs(pos[]-eta_m0) < cirp_th) {
        coord n;
	n = interface_normal(point, c);
	double deltaD = sqrt(sq(n.x) + sq(n.y) + sq(n.z));
	double ar     = deltaD*dv();
	area_my += ar;
	eta_my  += ar*pos[]; // already defined at the interface
      }
    }
  }
  area_my += 1.0e-12; // to prevent arithmetic failure if area = 0
  eta_my   = eta_my/area_my;
  fprintf(stderr, "area %8E, eta %8E\n", area_my, eta_my);

  // --> my diagnosis of eta_my, area_my, amp_my (we want to just do it at the interface)
  double amp2_my = 0;
  foreach(reduction(+:amp2_my)) {
    if (interfacial (point, c)) {
      //if (point.level == MAXLEVEL) {
      if (point.level == MAXLEVEL && abs(pos[]-eta_m0) < cirp_th) {
        coord n;
	n = interface_normal(point, c);
	double deltaD = sqrt(sq(n.x) + sq(n.y) + sq(n.z));
	double ar     = deltaD*dv();
        amp2_my += ar*2.0*sq(pos[]-eta_my);
      }
    }
  }
  amp2_my = amp2_my/area_my;

  // --> compute some mean quantities to be used later
  double area = 0; double eta_m = 0; double pr_m = 0;
  double u_xi = 0; double u_yi  = 0; double u_zi = 0;
  foreach(reduction(+:area), reduction(+:eta_m), reduction(+:pr_m),
	  reduction(+:u_xi), reduction(+:u_yi) , reduction(+:u_zi)) {
    if (interfacial (point, c)) {
      //if (point.level == MAXLEVEL) {
      if (point.level == MAXLEVEL && abs(pos[]-eta_m0) < cirp_th) {
        coord n, pp;
	n = interface_normal(point, c);
        double alpha1 = plane_alpha (c[], n);
        plane_area_center(n, alpha1, &pp);
	double deltaD = sqrt(sq(n.x) + sq(n.y) + sq(n.z));
	double ar     = deltaD*dv();
	double eta = pos[]; // must be here since here it is defined

	double yp = y + stp_eta;
	
#if dimension > 2
	Point point = locate (x, yp, z);
#else
	Point point = locate (x, yp);
#endif
	if (point.level > 0) {
        
	  POINT_VARIABLES;
	
	  area += ar;
	  eta_m += ar*eta; // already defined at the interface

	  double pr_m_t1 = 0.0, xc = 0.0, yc = 0.0, zc = 0.0;
	  xc = x + Delta*pp.x, yc = yp + Delta*pp.y, zc = z + Delta*pp.z;
          pr_m_t1 = my_interpolation(p, xc, yc, zc);

	  if(pr_m_t1 != nodata) { // best case: we can locate the centroid location
            pr_m += ar*pr_m_t1;
            u_xi += ar*my_interpolation(u.x, xc, yc, zc);
            u_yi += ar*my_interpolation(u.y, xc, yc, zc);
#if dimension > 2
            u_zi += ar*my_interpolation(u.z, xc, yc, zc);
#else
	    u_zi = 0.;
#endif
	  }
	  else { // we avoid interpolation
	    pr_m += ar*p[];
	    u_xi += ar*u.x[];
	    u_yi += ar*u.y[];
	    u_zi += ar*u.z[];
#if dimension > 2
            u_zi += ar*u.z[];
#else
	    u_zi += 0.;
#endif
	  }

	}
      }
    }
  }
  area  += 1.0e-12; // to prevent arithmetic failure if area = 0
  eta_m  = eta_m/area;
  pr_m   = pr_m/area;
  u_xi   = u_xi/area; u_yi = u_yi/area; u_zi = u_zi/area;
  //fprintf(stderr, "area %8E, eta %8E, pr_m %8E\n", area, eta_m, pr_m);
  //fprintf(stderr, "u.x %8E, u.y %8E, u.z %8E\n", u_xi, u_yi, u_zi);

  // --> quantities to be probed at the interface
  // n.b.: for the stress we can use u since it contains derivatives
  scalar Sxx[]; scalar Syy[]; scalar Szz[];
  scalar Sxy[]; scalar Sxz[]; scalar Syz[];
  foreach() {
    double dudx = (u.x[1]     - u.x[-1]    )/(2.*Delta);
    double dudy = (u.x[0,1]   - u.x[0,-1]  )/(2.*Delta);
    double dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
    double dvdx = (u.y[1]     - u.y[-1]    )/(2.*Delta);
    double dvdy = (u.y[0,1]   - u.y[0,-1]  )/(2.*Delta);
    double dvdz = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
    double dwdx = (u.z[1]     - u.z[-1]    )/(2.*Delta);
    double dwdy = (u.z[0,1]   - u.z[0,-1]  )/(2.*Delta);
    double dwdz = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);
    Sxx[] = dudx+dudx; 
    Syy[] = dvdy+dvdy;
    Szz[] = dwdz+dwdz;
    Sxy[] = dudy+dvdx;
    Sxz[] = dudz+dwdx;
    Syz[] = dvdz+dwdy;
  }
  boundary({Sxx,Syy,Szz,Sxy,Sxz,Syz});

  // --> compute the amplitude, stress, velocity at the interface, energy fluxes
  double amp2  = 0.0;
  double mf_px = 0.0; double mf_py = 0.0; double mf_pz = 0.0;
  double mf_vx = 0.0; double mf_vy = 0.0; double mf_vz = 0.0;
  double ef_px = 0.0; double ef_py = 0.0; double ef_pz = 0.0;
  double ef_vx = 0.0; double ef_vy = 0.0; double ef_vz = 0.0;
  double ef_pp = 0.0; double ef_vp = 0.0;
  foreach(reduction(+:amp2), // amplitude 
	  reduction(+:mf_px), reduction(+:mf_py), reduction(+:mf_pz), // pressure - momentum flux (x,y,z)
	  reduction(+:mf_vx), reduction(+:mf_vy), reduction(+:mf_vz), // viscous  - momentum flux (x,y,z)
	  reduction(+:ef_px), reduction(+:ef_py), reduction(+:ef_pz), // pressure - en flux (x,y,z)
	  reduction(+:ef_vx), reduction(+:ef_vy), reduction(+:ef_vz), // viscous  - en flux (x,y,z)
	  reduction(+:ef_pp), reduction(+:ef_vp)) { // total en flux
    if (interfacial (point, c)) {
      //if (point.level == MAXLEVEL) {
      if (point.level == MAXLEVEL && abs(pos[]-eta_m0) < cirp_th) {
        coord n, pp;
	n = interface_normal(point, c);
        double alpha1 = plane_alpha (c[], n);
        plane_area_center(n, alpha1, &pp);
	double deltaD = sqrt(sq(n.x) + sq(n.y) + sq(n.z));
	double ar     = deltaD*dv();
	double eta    = pos[]; // must be here since here it is defined
	normalize (&n); // |n| = 1, we should normalize since we use the normals for online calculations

	double yp = y + stp_eta;
	
#if dimension > 2
	Point point = locate (x, yp, z);
#else
	Point point = locate (x, yp);
#endif
	if (point.level > 0) {
        
	  POINT_VARIABLES;
	
	  amp2 += ar*2.0*sq(eta-eta_m); // amplitude
 
          double pr_int = 0;   
	  double Sxx_int = 0; double Syy_int = 0; double Szz_int = 0;
	  double Sxy_int = 0; double Sxz_int = 0; double Syz_int = 0;
	  double ux_int  = 0; double uy_int  = 0; double uz_int  = 0;

	  double pr_int_t1 = 0.0, xc = 0.0, yc = 0.0, zc = 0.0;
	  xc = x + Delta*pp.x, yc = yp + Delta*pp.y, zc = z + Delta*pp.z;
          pr_int_t1 = my_interpolation(p, xc, yc, zc);

	  if(pr_int_t1 != nodata) { // best case: we can locate the centroid location
            //printf("Case 0\n");
            pr_int  = pr_int_t1;
	    Sxx_int = my_interpolation(Sxx, xc, yc, zc);
	    Syy_int = my_interpolation(Syy, xc, yc, zc);
	    Szz_int = my_interpolation(Szz, xc, yc, zc);
	    Sxy_int = my_interpolation(Sxy, xc, yc, zc);
	    Sxz_int = my_interpolation(Sxz, xc, yc, zc);
	    Syz_int = my_interpolation(Syz, xc, yc, zc);
	    ux_int  = my_interpolation(u.x, xc, yc, zc);
	    uy_int  = my_interpolation(u.y, xc, yc, zc);
#if dimension > 2
            uz_int  = my_interpolation(u.z, xc, yc, zc);
#else
	    uz_int  = 0.;
#endif
	  }
          else { // we avoid interpolation
            //printf("Case 2\n");
	    pr_int  = p[];
	    Sxx_int = Sxx[];
	    Syy_int = Syy[];
	    Szz_int = Szz[];
	    Sxy_int = Sxy[];
	    Sxz_int = Sxz[];
	    Syz_int = Syz[];
	    ux_int  = u.x[];
	    uy_int  = u.y[];
#if dimension > 2          
            uz_int  = u.z[];
#else
	    uz_int  = 0.;
#endif
	  }

	  mf_px += ar*( -(pr_int-pr_m)*n.x ); // x-Mom. flux due to pressure work
	  mf_py += ar*( -(pr_int-pr_m)*n.y ); // y-Mom. flux due to pressure work
	  mf_pz += ar*( -(pr_int-pr_m)*n.z ); // z-Mom. flux due to pressure work
          mf_vx += ar*mu2*(n.x*Sxx_int + n.y*Sxy_int + n.z*Sxz_int); // x-Mom. flux due to viscous dissipation
          mf_vy += ar*mu2*(n.x*Sxy_int + n.y*Syy_int + n.z*Syz_int); // y-Mom. flux due to viscous dissipation
          mf_vz += ar*mu2*(n.x*Sxz_int + n.y*Syz_int + n.z*Szz_int); // z-Mom. flux due to viscous dissipation
          
	  /*
          ef_p  += ar*( -pr_int*(ux_int*n.x+uy_int*n.y+uz_int*n.z) ); // En. flux due to pressure
          ef_v  += ar*mu2*( ( Sxx_int*n.x + Sxy_int*n.y + Sxz_int*n.z )*ux_int +
	                    ( Sxy_int*n.x + Syy_int*n.y + Syz_int*n.z )*uy_int +
	                    ( Sxz_int*n.x + Syz_int*n.y + Szz_int*n.z )*uz_int ); // En. flux due to viscous dissipation
          */

          ef_px += ar*( -(pr_int-pr_m)*(ux_int-u_xi) )*n.x; // x contr. to En. flux pressure
          ef_py += ar*( -(pr_int-pr_m)*(uy_int-u_yi) )*n.y; // y contr. to En. flux pressure
          ef_pz += ar*( -(pr_int-pr_m)*(uz_int-u_zi) )*n.z; // z contr. to En. flux pressure

          ef_vx += ar*mu2*(n.x*Sxx_int + n.y*Sxy_int + n.z*Sxz_int)*(ux_int-u_xi); // x contr. to En. flux diss.
          ef_vy += ar*mu2*(n.x*Sxy_int + n.y*Syy_int + n.z*Syz_int)*(uy_int-u_yi); // y contr. to En. flux diss.
          ef_vz += ar*mu2*(n.x*Sxz_int + n.y*Syz_int + n.z*Szz_int)*(uz_int-u_zi); // z contr. to En. flux diss.

	  ef_pp += ar*( -(pr_int-pr_m)*( (ux_int-u_xi)*n.x+(uy_int-u_yi)*n.y+(uz_int-u_zi)*n.z) ); // En. flux due to pressure
          ef_vp += ar*mu2*( ( Sxx_int*n.x + Sxy_int*n.y + Sxz_int*n.z )*(ux_int-u_xi) +
	                    ( Sxy_int*n.x + Syy_int*n.y + Syz_int*n.z )*(uy_int-u_yi) +
	                    ( Sxz_int*n.x + Syz_int*n.y + Szz_int*n.z )*(uz_int-u_zi) ); // En. flux due to viscous dissipation
 
	}
      }
    }
  }
  amp2  = amp2/area; 
  mf_px = mf_px/area; mf_py = mf_py/area; mf_pz = mf_pz/area;
  mf_vx = mf_vx/area; mf_vy = mf_vy/area; mf_vz = mf_vz/area;
  ef_px = ef_px/area; ef_py = ef_py/area; ef_pz = ef_pz/area;
  ef_vx = ef_vx/area; ef_vy = ef_vy/area; ef_vz = ef_vz/area;
  ef_pp = ef_pp/area; ef_vp = ef_vp/area;

  if (pid() == 0) {
  
    fflush(stderr);
    FILE * global_obs = fopen (fname, "a");
    fprintf (global_obs, "%8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E\n", 
                         t, 1.0*istep, area, eta_m, k_*sqrt(amp2), area_my, eta_my, k_*sqrt(amp2_my), pr_m,  
			 mf_px, mf_py, mf_pz, mf_vx, mf_vy, mf_vz, 
			 u_xi, u_yi, u_zi, 
			 ef_px, ef_py, ef_pz, ef_vx, ef_vy, ef_vz, 
			 ef_pp, ef_vp);
    fclose(global_obs);

  }

}

/**
   We want to compute some global observables - second way */

void output_global_obs_2 (char * fname, int istep, scalar c, double stp_eta) {

  // Preliminary calculations
 
  // --> position
  scalar pos[];
#if dimension > 2
  coord G = {0.,1.,0.}, Z = {0.,0.,0.};
#else
  coord G = {0.,1.}, Z = {0.,0.};
#endif
  position (c, pos, G, Z);
  boundary ({pos});

  // --> my diagnosis of eta_my, area_my, amp_my
  double area_my = 0; double eta_my = 0;
  foreach(reduction(+:area_my), reduction(+:eta_my)) {
    if (interfacial (point, c)) {
      //if (point.level == MAXLEVEL) {
      if (point.level == MAXLEVEL && abs(pos[]-eta_m0) < cirp_th) {
        coord n;
	n = interface_normal(point, c);
	double deltaD = sqrt(sq(n.x) + sq(n.y) + sq(n.z));
	double ar     = deltaD*dv();
	area_my += ar;
	eta_my  += ar*pos[]; // already defined at the interface
      }
    }
  }
  area_my += 1.0e-12; // to prevent arithmetic failure if area = 0
  eta_my   = eta_my/area_my;
  fprintf(stderr, "area %8E, eta %8E\n", area_my, eta_my);

  // --> my diagnosis of eta_my, area_my, amp_my (we want to just do it at the interface)
  double amp2_my = 0;
  foreach(reduction(+:amp2_my)) {
    if (interfacial (point, c)) {
      //if (point.level == MAXLEVEL) {
      if (point.level == MAXLEVEL && abs(pos[]-eta_m0) < cirp_th) {
        coord n;
	n = interface_normal(point, c);
	double deltaD = sqrt(sq(n.x) + sq(n.y) + sq(n.z));
	double ar     = deltaD*dv();
        amp2_my += ar*2.0*sq(pos[]-eta_my);
      }
    }
  }
  amp2_my = amp2_my/area_my;

  // --> compute some mean quantities to be used later
  double area = 0; double eta_m = 0; double pr_m = 0;
  double u_xi = 0; double u_yi  = 0; double u_zi = 0;
  foreach(reduction(+:area), reduction(+:eta_m), reduction(+:pr_m),
	  reduction(+:u_xi), reduction(+:u_yi) , reduction(+:u_zi)) {
    if (interfacial (point, c)) {
      //if (point.level == MAXLEVEL) {
      if (point.level == MAXLEVEL && abs(pos[]-eta_m0) < cirp_th) {
        coord n, pp;
	n = interface_normal(point, c);
        double alpha1 = plane_alpha (c[], n);
        plane_area_center(n, alpha1, &pp);
	double deltaD = sqrt(sq(n.x) + sq(n.y) + sq(n.z));
	double ar     = deltaD*dv();
	double eta = pos[]; // must be here since here it is defined

	double yp = y + stp_eta;
	/*
	//double yp = y + off_set*Delta;
#if dimension > 2
	Point point = locate (x, yp, z);
#else
	Point point = locate (x, yp);
#endif
	if (point.level > 0) {
        
	  POINT_VARIABLES;
	*/

	  area += ar;
	  eta_m += ar*eta; // already defined at the interface

	  double pr_m_t1 = 0.0, xc = 0.0, yc = 0.0, zc = 0.0;
	  xc = x + Delta*pp.x, yc = yp + Delta*pp.y, zc = z + Delta*pp.z;
          pr_m_t1 = my_interpolation(p, xc, yc, zc);

	  if(pr_m_t1 != nodata) { // best case: we can locate the centroid location
            pr_m += ar*pr_m_t1;
            u_xi += ar*my_interpolation(u.x, xc, yc, zc);
            u_yi += ar*my_interpolation(u.y, xc, yc, zc);
#if dimension > 2
            u_zi += ar*my_interpolation(u.z, xc, yc, zc);
#else
	    u_zi = 0.;
#endif
	  }
	  else { // we avoid interpolation
	    pr_m += ar*p[];
	    u_xi += ar*u.x[];
	    u_yi += ar*u.y[];
	    u_zi += ar*u.z[];
#if dimension > 2
            u_zi += ar*u.z[];
#else
	    u_zi += 0.;
#endif
	  }

	//}
      }
    }
  }
  area  += 1.0e-12; // to prevent arithmetic failure if area = 0
  eta_m  = eta_m/area;
  pr_m   = pr_m/area;
  u_xi   = u_xi/area; u_yi = u_yi/area; u_zi = u_zi/area;
  //fprintf(stderr, "area %8E, eta %8E, pr_m %8E\n", area, eta_m, pr_m);
  //fprintf(stderr, "u.x %8E, u.y %8E, u.z %8E\n", u_xi, u_yi, u_zi);

  // --> quantities to be probed at the interface
  // n.b.: for the stress we can use u since it contains derivatives
  scalar Sxx[]; scalar Syy[]; scalar Szz[];
  scalar Sxy[]; scalar Sxz[]; scalar Syz[];
  foreach() {
    double dudx = (u.x[1]     - u.x[-1]    )/(2.*Delta);
    double dudy = (u.x[0,1]   - u.x[0,-1]  )/(2.*Delta);
    double dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
    double dvdx = (u.y[1]     - u.y[-1]    )/(2.*Delta);
    double dvdy = (u.y[0,1]   - u.y[0,-1]  )/(2.*Delta);
    double dvdz = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
    double dwdx = (u.z[1]     - u.z[-1]    )/(2.*Delta);
    double dwdy = (u.z[0,1]   - u.z[0,-1]  )/(2.*Delta);
    double dwdz = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);
    Sxx[] = dudx+dudx; 
    Syy[] = dvdy+dvdy;
    Szz[] = dwdz+dwdz;
    Sxy[] = dudy+dvdx;
    Sxz[] = dudz+dwdx;
    Syz[] = dvdz+dwdy;
  }
  boundary({Sxx,Syy,Szz,Sxy,Sxz,Syz});

  // --> compute the amplitude, stress, velocity at the interface, energy fluxes
  double amp2  = 0.0;
  double mf_px = 0.0; double mf_py = 0.0; double mf_pz = 0.0;
  double mf_vx = 0.0; double mf_vy = 0.0; double mf_vz = 0.0;
  double ef_px = 0.0; double ef_py = 0.0; double ef_pz = 0.0;
  double ef_vx = 0.0; double ef_vy = 0.0; double ef_vz = 0.0;
  double ef_pp = 0.0; double ef_vp = 0.0;
  foreach(reduction(+:amp2), // amplitude 
	  reduction(+:mf_px), reduction(+:mf_py), reduction(+:mf_pz), // pressure - momentum flux (x,y,z)
	  reduction(+:mf_vx), reduction(+:mf_vy), reduction(+:mf_vz), // viscous  - momentum flux (x,y,z)
	  reduction(+:ef_px), reduction(+:ef_py), reduction(+:ef_pz), // pressure - en flux (x,y,z)
	  reduction(+:ef_vx), reduction(+:ef_vy), reduction(+:ef_vz), // viscous  - en flux (x,y,z)
	  reduction(+:ef_pp), reduction(+:ef_vp)) { // total en flux
    if (interfacial (point, c)) {
      //if (point.level == MAXLEVEL) {
      if (point.level == MAXLEVEL && abs(pos[]-eta_m0) < cirp_th) {
        coord n, pp;
	n = interface_normal(point, c);
        double alpha1 = plane_alpha (c[], n);
        plane_area_center(n, alpha1, &pp);
	double deltaD = sqrt(sq(n.x) + sq(n.y) + sq(n.z));
	double ar     = deltaD*dv();
	double eta    = pos[]; // must be here since here it is defined
	normalize (&n); // |n| = 1, we should normalize since we use the normals for online calculations

	double yp = y + stp_eta;
	/*
	//double yp = y + off_set*Delta;
#if dimension > 2
	Point point = locate (x, yp, z);
#else
	Point point = locate (x, yp);
#endif
	if (point.level > 0) {
        
	  POINT_VARIABLES;
	*/

	  amp2 += ar*2.0*sq(eta-eta_m); // amplitude
 
          double pr_int = 0;   
	  double Sxx_int = 0; double Syy_int = 0; double Szz_int = 0;
	  double Sxy_int = 0; double Sxz_int = 0; double Syz_int = 0;
	  double ux_int  = 0; double uy_int  = 0; double uz_int  = 0;

	  double pr_int_t1 = 0.0, xc = 0.0, yc = 0.0, zc = 0.0;
	  xc = x + Delta*pp.x, yc = yp + Delta*pp.y, zc = z + Delta*pp.z;
          pr_int_t1 = my_interpolation(p, xc, yc, zc);

	  if(pr_int_t1 != nodata) { // best case: we can locate the centroid location
            //printf("Case 0\n");
            pr_int  = pr_int_t1;
	    Sxx_int = my_interpolation(Sxx, xc, yc, zc);
	    Syy_int = my_interpolation(Syy, xc, yc, zc);
	    Szz_int = my_interpolation(Szz, xc, yc, zc);
	    Sxy_int = my_interpolation(Sxy, xc, yc, zc);
	    Sxz_int = my_interpolation(Sxz, xc, yc, zc);
	    Syz_int = my_interpolation(Syz, xc, yc, zc);
	    ux_int  = my_interpolation(u.x, xc, yc, zc);
	    uy_int  = my_interpolation(u.y, xc, yc, zc);
#if dimension > 2
            uz_int  = my_interpolation(u.z, xc, yc, zc);
#else
	    uz_int  = 0.;
#endif
	  }
          else { // we avoid interpolation
            //printf("Case 2\n");
	    pr_int  = p[];
	    Sxx_int = Sxx[];
	    Syy_int = Syy[];
	    Szz_int = Szz[];
	    Sxy_int = Sxy[];
	    Sxz_int = Sxz[];
	    Syz_int = Syz[];
	    ux_int  = u.x[];
	    uy_int  = u.y[];
#if dimension > 2          
            uz_int  = u.z[];
#else
	    uz_int  = 0.;
#endif
	  }

	  mf_px += ar*( -(pr_int-pr_m)*n.x ); // x-Mom. flux due to pressure work
	  mf_py += ar*( -(pr_int-pr_m)*n.y ); // y-Mom. flux due to pressure work
	  mf_pz += ar*( -(pr_int-pr_m)*n.z ); // z-Mom. flux due to pressure work
          mf_vx += ar*mu2*(n.x*Sxx_int + n.y*Sxy_int + n.z*Sxz_int); // x-Mom. flux due to viscous dissipation
          mf_vy += ar*mu2*(n.x*Sxy_int + n.y*Syy_int + n.z*Syz_int); // y-Mom. flux due to viscous dissipation
          mf_vz += ar*mu2*(n.x*Sxz_int + n.y*Syz_int + n.z*Szz_int); // z-Mom. flux due to viscous dissipation
          
	  /*
          ef_p  += ar*( -pr_int*(ux_int*n.x+uy_int*n.y+uz_int*n.z) ); // En. flux due to pressure
          ef_v  += ar*mu2*( ( Sxx_int*n.x + Sxy_int*n.y + Sxz_int*n.z )*ux_int +
	                    ( Sxy_int*n.x + Syy_int*n.y + Syz_int*n.z )*uy_int +
	                    ( Sxz_int*n.x + Syz_int*n.y + Szz_int*n.z )*uz_int ); // En. flux due to viscous dissipation
          */

          ef_px += ar*( -(pr_int-pr_m)*(ux_int-u_xi) )*n.x; // x contr. to En. flux pressure
          ef_py += ar*( -(pr_int-pr_m)*(uy_int-u_yi) )*n.y; // y contr. to En. flux pressure
          ef_pz += ar*( -(pr_int-pr_m)*(uz_int-u_zi) )*n.z; // z contr. to En. flux pressure

          ef_vx += ar*mu2*(n.x*Sxx_int + n.y*Sxy_int + n.z*Sxz_int)*(ux_int-u_xi); // x contr. to En. flux diss.
          ef_vy += ar*mu2*(n.x*Sxy_int + n.y*Syy_int + n.z*Syz_int)*(uy_int-u_yi); // y contr. to En. flux diss.
          ef_vz += ar*mu2*(n.x*Sxz_int + n.y*Syz_int + n.z*Szz_int)*(uz_int-u_zi); // z contr. to En. flux diss.

	  ef_pp += ar*( -(pr_int-pr_m)*( (ux_int-u_xi)*n.x+(uy_int-u_yi)*n.y+(uz_int-u_zi)*n.z) ); // En. flux due to pressure
          ef_vp += ar*mu2*( ( Sxx_int*n.x + Sxy_int*n.y + Sxz_int*n.z )*(ux_int-u_xi) +
	                    ( Sxy_int*n.x + Syy_int*n.y + Syz_int*n.z )*(uy_int-u_yi) +
	                    ( Sxz_int*n.x + Syz_int*n.y + Szz_int*n.z )*(uz_int-u_zi) ); // En. flux due to viscous dissipation
 
	//}
      }
    }
  }
  amp2  = amp2/area; 
  mf_px = mf_px/area; mf_py = mf_py/area; mf_pz = mf_pz/area;
  mf_vx = mf_vx/area; mf_vy = mf_vy/area; mf_vz = mf_vz/area;
  ef_px = ef_px/area; ef_py = ef_py/area; ef_pz = ef_pz/area;
  ef_vx = ef_vx/area; ef_vy = ef_vy/area; ef_vz = ef_vz/area;
  ef_pp = ef_pp/area; ef_vp = ef_vp/area;

  if (pid() == 0) {
  
    fflush(stderr);
    FILE * global_obs = fopen (fname, "a");
    fprintf (global_obs, "%8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E\n", 
                         t, 1.0*istep, area, eta_m, k_*sqrt(amp2), area_my, eta_my, k_*sqrt(amp2_my), pr_m,  
			 mf_px, mf_py, mf_pz, mf_vx, mf_vy, mf_vz, 
			 u_xi, u_yi, u_zi, 
			 ef_px, ef_py, ef_pz, ef_vx, ef_vy, ef_vz, 
			 ef_pp, ef_vp);
    fclose(global_obs);

  }

}


//event bulk_energy (t += T0_/tout_glo) {
event bulk_energy (t += T0_/tout_glo_my) {
//event bulk_energy (i += 2) {

  fprintf(stderr, "I output the bulk energy and dissipation every t=%8E\n", T0_/tout_glo_my);

  // --> Kinetic, Potential energy, dissipation of water and air (bulk)
  double ke = 0., gpe = 0.;
  double keAir = 0., gpeAir = 0.;
  foreach(reduction(+:ke) reduction(+:gpe) 
	  reduction(+:keAir) reduction(+:gpeAir)) {
    double norm2 = 0.;
    foreach_dimension() {
      norm2 += sq(u.x[]);
    }
    ke     += rho[]*norm2*f[]*dv();
    keAir  += rho[]*norm2*(1.0-f[])*dv();
    gpe    += rho1*g_*y*f[]*dv();
    gpeAir += rho2*g_*y*(1.0-f[])*dv();
  }
  double rates[2];
  dissipation_rate(rates);
  double dissWater = rates[0];
  double dissAir   = rates[1];
 
  if(pid() == 0) {

    fflush(stderr);
    char name_2[80], name_3[80];
    sprintf (name_2,"budgetWaterwind.out");
    sprintf (name_3,"budgetAirwind.out");
    FILE * fpwater = fopen(name_2,"a");
    FILE * fpair   = fopen(name_3,"a");
    fprintf (fpwater, "%8E %8E %8E %8E %8E\n",
             t, 1.0*i, ke/2., gpe + 0.125*L0*L0*L0, dissWater);
    fprintf (fpair, "%8E %8E %8E %8E %8E\n",
             t, 1.0*i, keAir/2., gpeAir + 0.125*L0*L0*L0, dissAir);
    fclose(fpwater);
    fclose(fpair);

  }

}

int compare (const void *a, const void *b) {
  
  double x = *(double *)a;
  double y = *(double *)b;

  if(     x < y) { 
    return -1;
  }
  else if(x > y) {
    return +1;
  }
  else {
    return +0;
  }

}

/**
   We want to compute some quantities at the interface. */

void output_int_qtn (char * fname, int istep, scalar c, double stp_eta) {

  fprintf(stderr, "Preliminary calculation\n");

  // We first loop over all the interfacial points 
  // and we count them (per processor)
  
  int int_pt = 0;
  foreach(serial) { 
    if (interfacial (point, f)) {
      if (point.level == MAXLEVEL) {
	double yp = y + stp_eta;
#if dimension > 2
	Point point = locate (x, yp, z);
#else
	Point point = locate (x, yp);
#endif
	if (point.level > 0) {
	  POINT_VARIABLES;
	  int_pt++;
	}
      }
    }
  }

  //int tot_column = 21;
  int tot_column = 18;
  double t_mat[int_pt][tot_column];

  for (int j = 0; j < tot_column; j++) {
    for (int i = 0; i < int_pt; i++) {
      t_mat[i][j] = 0;
    }
  }
  fprintf(stderr, "First pass over interfacial points\n");

  scalar pos[];
#if dimension > 2
  coord G = {0.,1.,0.}, Z = {0.,0.,0.};
#else
  coord G = {0.,1.}, Z = {0.,0.};
#endif
  position (c, pos, G, Z);
  boundary ({pos});

  scalar curv[];
  curvature (c, curv);

  // --> quantities to be probed at the interface
  // n.b.: for the stress we can use u since it contains derivatives
  scalar Sxx[]; scalar Syy[]; scalar Szz[];
  scalar Sxy[]; scalar Sxz[]; scalar Syz[];
  foreach() {
    double dudx = (u.x[1]     - u.x[-1]    )/(2.*Delta);
    double dudy = (u.x[0,1]   - u.x[0,-1]  )/(2.*Delta);
    double dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
    double dvdx = (u.y[1]     - u.y[-1]    )/(2.*Delta);
    double dvdy = (u.y[0,1]   - u.y[0,-1]  )/(2.*Delta);
    double dvdz = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
    double dwdx = (u.z[1]     - u.z[-1]    )/(2.*Delta);
    double dwdy = (u.z[0,1]   - u.z[0,-1]  )/(2.*Delta);
    double dwdz = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);
    Sxx[] = dudx+dudx; 
    Syy[] = dvdy+dvdy;
    Szz[] = dwdz+dwdz;
    Sxy[] = dudy+dvdx;
    Sxz[] = dudz+dwdx;
    Syz[] = dvdz+dwdy;
  }
  boundary({Sxx,Syy,Szz,Sxy,Sxz,Syz});

  int count = 0;
  foreach(serial) {
    if (interfacial (point, c)) {
      if (point.level == MAXLEVEL) {
        coord n, pp;
	n = interface_normal(point, c);
        double alpha1 = plane_alpha (c[], n);
        plane_area_center(n, alpha1, &pp);
	double eta = pos[];
	double cur = curv[];
	double norm_2 = sqrt(sq(n.x) + sq(n.y) + sq(n.z));
	double yp = y + stp_eta;
#if dimension > 2
	Point point = locate (x, yp, z);
#else
	Point point = locate (x, yp);
#endif
	if (point.level > 0) {
	 
	  POINT_VARIABLES;

          double pr_int_t1 = 0, xc = 0, yc = 0, zc = 0;
	  xc = x + Delta*pp.x, yc = yp + Delta*pp.y, zc = z + Delta*pp.z;
          pr_int_t1 = my_interpolation(p, xc, yc, zc);

          if(pr_int_t1 != nodata) { // best case: we can locate the centroid location
            //printf("Case 0\n");
	    t_mat[count][0]  = xc;
#if dimension > 2
	    t_mat[count][1]  = zc;
#else
	    t_mat[count][1]  = 0;
#endif
	    t_mat[count][2]  = pr_int_t1;
	    t_mat[count][3]  = my_interpolation(Sxx, xc, yc, zc);
	    t_mat[count][4]  = my_interpolation(Syy, xc, yc, zc);
	    t_mat[count][5]  = my_interpolation(Szz, xc, yc, zc);
	    t_mat[count][6]  = my_interpolation(Sxy, xc, yc, zc);
	    t_mat[count][7]  = my_interpolation(Sxz, xc, yc, zc);
	    t_mat[count][8]  = my_interpolation(Syz, xc, yc, zc);
	    t_mat[count][9]  = my_interpolation(u.x, xc, yc, zc);
	    t_mat[count][10] = my_interpolation(u.y, xc, yc, zc);
#if dimension > 2
	    t_mat[count][11] = my_interpolation(u.z, xc, yc, zc);
#else
	    t_mat[count][11] = 0;
#endif
	  }
	  else {
	    t_mat[count][0]  = x;
#if dimension > 2
	    t_mat[count][1]  = z;
#else
	    t_mat[count][1]  = 0;
#endif
	    t_mat[count][2]  = p[];
	    t_mat[count][3]  = Sxx[];
	    t_mat[count][4]  = Syy[];
	    t_mat[count][5]  = Szz[];
	    t_mat[count][6]  = Sxy[];
	    t_mat[count][7]  = Sxz[];
	    t_mat[count][8]  = Syz[];
	    t_mat[count][9]  = u.x[];
	    t_mat[count][10] = u.y[];
#if dimension > 2
	    t_mat[count][11] = u.z[];
#else
	    t_mat[count][11] = 0;
#endif
	  }

	  t_mat[count][12] = eta;           // already defined at the interface
	  t_mat[count][13] = cur;           // already defined at the interface
	  t_mat[count][14] = n.x/norm_2;    // already defined at the interface
	  t_mat[count][15] = n.y/norm_2;    // already defined at the interface
	  t_mat[count][16] = n.z/norm_2;    // already defined at the interface
	  t_mat[count][17] = Delta;         // no interpolation is meaningful
	  count++;

	}
      }
    }
  }
  fprintf(stderr, "Second pass over interfacial points\n");

  // We sort locally t_mat by the x coordinate (the first, i.e. 0-th, column) 
  qsort(t_mat, int_pt, tot_column*sizeof(double), compare);
  fprintf(stderr, "First sort\n");

#if _MPI

  // On multiple cores, we gather all the int_pt to the root pid 
  // and, then, we broadcast this information to all the processes
  
  int nproc;
  MPI_Comm_size (MPI_COMM_WORLD, &nproc);
  int counts_it[nproc];
  if( pid() == 0) {
    MPI_Gather(&int_pt,1,MPI_INT,counts_it,1,MPI_INT,0,MPI_COMM_WORLD); // MPI_gather gathers by rank order
  }
  else {
    MPI_Gather(&int_pt,1,MPI_INT,NULL     ,1,MPI_INT,0,MPI_COMM_WORLD); // MPI_gather gathers by rank order
  }
  MPI_Bcast(counts_it,nproc,MPI_INT,0,MPI_COMM_WORLD);

  // Each processor knows the int_pt of the others. So we can compute 
  // the displacement, the total number of interfacial points and 
  // the number of elements owned by each processor

  int tot_int_pt = 0;
  int tot_el_p[nproc], disp_r[nproc], disp[nproc];
  for (int i = 0; i < nproc; i++) {
    disp_r[i]    = tot_int_pt;
    disp[i]      = disp_r[i]*tot_column;
    tot_int_pt  += counts_it[i];
    tot_el_p[i]  = counts_it[i]*tot_column;
  }

  // --> Gather to the root pid 
  // --> Sort by first column

  double t_mat_tot[tot_int_pt][tot_column];

  if( pid() == 0 ) {
    //int tot_el_p0[nproc], disp_r0[nproc], disp0[nproc];
    int tot_el_p0[nproc], disp0[nproc];
    for (int i = 0; i < nproc; i++) {
      tot_el_p0[i] = tot_el_p[i];
      //disp_r0[i]   = disp_r[i];
      disp0[i]     = disp[i];
    }
    MPI_Gatherv(&t_mat,int_pt*tot_column,MPI_DOUBLE,t_mat_tot,tot_el_p0,disp0,MPI_DOUBLE,0,MPI_COMM_WORLD);
    fprintf(stderr, "GatherV 0\n");
    qsort(t_mat_tot, tot_int_pt, tot_column*sizeof(double), compare);
    fprintf(stderr, "Second sort 0\n");
  }
  else {
    int tot_el_p1[nproc], disp1[nproc];
    for (int i = 0; i < nproc; i++) {
      tot_el_p1[i] = tot_el_p[i];
      disp1[i]     = disp[i];
    }
    MPI_Gatherv(&t_mat,int_pt*tot_column,MPI_DOUBLE,NULL     ,tot_el_p1,disp1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    //printf("GatherV Non 0\n");
  }

#else

  for (int j = 0; j < tot_column; j++) {
    for (int i = 0; i < int_pt; i++) {
      t_mat_tot[i][j] = t_mat[i][j];
    }
  }

#endif

  if (pid() == 0) {
    
    fflush(stderr);

    // Print in binary format
    FILE * fpver = fopen (fname,"w");
    for (int i = 0; i < tot_int_pt; i++) {
      for (int j = 0; j < tot_column; j++) {
	fwrite ( &t_mat_tot[i][j], sizeof(double), 1, fpver );
      }
    }
    fclose(fpver);
    fflush(fpver);

    // Log some info about the binary
    char name_1[100];
    sprintf (name_1, "./eta/global_int.out");
    FILE * global_int = fopen (name_1, "a");
    fprintf (global_int, "%8E %09d %09d %09d\n", 
		          t, istep, tot_int_pt, tot_column);
    fclose(global_int);

  }
  fprintf(stderr, "Print\n");

}

/**
We also want to count the drops and bubbles in the flow. */

void countDropsBubble (char * name_1, char * name_2, char * name_3, int istep, scalar c) {

  scalar m1[]; // droplets
  scalar m2[]; // bubbles
  foreach() {
    m1[] = c[] > 0.5; // i.e. set m true if f[] is close to unity (droplets)
    m2[] = c[] < 0.5; // m true if f[] close to zero (bubbles)
  }
  int n1 = tag(m1);
  int n2 = tag(m2);
  
  /**
  Having counted the bubbles, now we find their size. This example
  is similar to the jet atomization problem. We are interested in
  the volumes and positions of each droplet/bubble.*/

  double v1[n1]; // droplet
  //coord b1[n1];  // droplet
  double v2[n2]; // bubble
  //coord b2[n2];  // bubble
  
  /**
  We initialize: */
  
  for (int j=0; j<n1; j++)
      //v1[j] = b1[j].x = b1[j].y = b1[j].z = 0.0;
      v1[j] = 0.0;
  for (int j=0; j<n2; j++)
      //v2[j] = b2[j].x = b2[j].y = b2[j].z = 0.0;
      v2[j] = 0.0;
  
  /**
  We proceed with calculation. */

  foreach_leaf() {

    // droplets
    if (m1[] > 0) {
      int j = m1[] - 1;
      v1[j] += dv()*c[]; //increment the volume of the droplet
      /*
      coord p = {x,y,z};
      foreach_dimension()
	b1[j].x += dv()*c[]*p.x; */
    }
    // bubbles
    if (m2[] > 0) {
      int j = m2[] - 1;
      v2[j] += dv()*(1.0-c[]);
      /* coord p = {x,y,z};
      foreach_dimension()
	b2[j].x += dv()*(1.0-c[])*p.x; */
    }

  }
  
  /**
  Reduce for MPI. */ 
  
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, v1, n1  , MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  fprintf(stderr, "I reduced v1\n");
  //MPI_Allreduce (MPI_IN_PLACE, b1, 3*n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  //fprintf(stderr, "I reduced b1\n");
  MPI_Allreduce (MPI_IN_PLACE, v2, n2  , MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  fprintf(stderr, "I reduced v2\n");
  //MPI_Allreduce (MPI_IN_PLACE, b2, 3*n2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  //fprintf(stderr, "I reduced b2\n");
#endif

  if (pid() == 0) {

    fflush(stderr);

    /**
    We first output the number of droplets and bubbles. */
    
    FILE * ftot = fopen(name_1,"a");
    fprintf (ftot, "%8E %8E %8E %8E\n", t, 1.0*istep, 1.0*(n1-1), 1.0*(n2-1)); // we remove the main region
    fclose(ftot);
    
    /**
    We output separately the volume and position of each droplet and bubble to file. */

    FILE * fdrop = fopen(name_2,"a");
    for (int j=0; j<n1; j++) {
      /* fprintf (fdrop, "%d %8E %8E %8E %8E\n", 
	                 j, v1[j], b1[j].x/v1[j], b1[j].y/v1[j], b1[j].z/v1[j]); */
      fprintf (fdrop, "%d %8E\n", j, v1[j]);
    }
    fclose(fdrop);
    fflush(fdrop);

    FILE * fbubb = fopen(name_3,"a");
    for (int j=0; j<n2; j++) {
      /* fprintf (fbubb, "%d %8E %8E %8E %8E\n", 
	                j, v2[j], b2[j].x/v2[j], b2[j].y/v2[j], b2[j].z/v2[j]); */
      fprintf (fdrop, "%d %8E\n", j, v2[j]);
    }
    fclose(fbubb);
    fflush(fbubb);
  
  }
	     
}

//event output_glo_obs (t += T0_/tout_glo) {
event output_glo_obs (t += T0_/tout_glo_my) {
//event output_glo_obs (i += 2) {

  scalar f2[];
  scalar_clone (f2, f);
  foreach() {
    f2[] = f[];
  }
  boundary ({f2});

  if(do_remove) {
    remove_droplets (f2,min_size,threshold,false); // first remove droplets
    remove_droplets (f2,min_size,threshold,true);  // then remove bubbles
  }

  fprintf(stderr, "I output the global observables file every t=%8E\n", T0_/tout_glo_my);

  char name_0[100];

  double stp_eta_1 = L0/256.0; // on L10 grid, it corresponds to 4*L0/2**10 = 4*Delta
  sprintf (name_0, "./global_obs_1.out");
  output_global_obs_1 (name_0, i, f2, stp_eta_1);
  
  double stp_eta_2 = 0.0; // we removed the point.level > 0
  sprintf (name_0, "./global_obs_2.out");
  output_global_obs_2 (name_0, i, f2, stp_eta_2);

  char name_1[100], name_2[100], name_3[100];
  
  fprintf(stderr, "I do first tagging every t=%8E\n", T0_/tout_glo_my);

  sprintf (name_1, "./num_bubbles_droplets_f1.out");
  sprintf (name_2, "./tagging_f1/droplets_f1_%09d.out", i);
  sprintf (name_3, "./tagging_f1/bubbles_f1_%09d.out", i);
  countDropsBubble (name_1, name_2, name_3, i, f);
  
  fprintf(stderr, "I do second tagging every t=%8E\n", T0_/tout_glo_my);
  
  sprintf (name_1, "./num_bubbles_droplets_f2.out");
  sprintf (name_2, "./tagging_f2/droplets_f2_%09d.out", i);
  sprintf (name_3, "./tagging_f2/bubbles_f2_%09d.out", i);
  countDropsBubble (name_1, name_2, name_3, i, f2);
}

//event eta_loc (t += T0_/tout_eta) {
event eta_loc (t += T0_/tout_eta_my) {
//event eta_loc (i += 2) {

  if (do_eta_loc == 1) {

    //fprintf(stderr, "I output the eta_loc file every t=%8E\n", T0_/tout_eta);
    fprintf(stderr, "I output the eta_loc file every t=%8E\n", T0_/tout_eta_my);
    //output_twophase_locate (i);
   
    /*
    scalar f2[];
    scalar_clone (f2, f);
    foreach() {
      f2[] = f[];
    }
    boundary ({f2});
    
    if(do_remove) {
      remove_droplets (f2,min_size,threshold,false); // first remove droplets
      remove_droplets (f2,min_size,threshold,true);  // then remove bubbles
    }
    output_int_qtn (i, f2);
    */
    
    fflush(stderr);

    char eta_out[100];
    sprintf (eta_out, "./eta/eta_loc/eta_loc_t%09d.bin", i);
    double stp_eta = L0/256.0; // on L10 grid, it corresponds to 4*L0/2**10 = 4*Delta 
    output_int_qtn (eta_out, i, f, stp_eta);
 
    if (pid() == 0) {
    
      char name_1[80];
      sprintf (name_1,"./eta/log_eta.out");
      FILE * log_sim = fopen(name_1,"a");
      fprintf (log_sim, "%8E %8E\n", t, 1.0*i);
      fclose(log_sim);
    
    }

  }

}

/**
## End 

   We want to run up end_sim physical time. */

event end (t = end_sim*T0_) {

  fprintf (fout, "i = %d t = %8E\n", i, t); fflush(fout);
  dump ("end.bin");

  if ( pid() == 0 ) {

    char comm[80];
    sprintf (comm, "ln -sf end.bin restart.bin");
    system(comm);

  }

  return 1; // exit

}

/** 
   Dump every tout_res period */

//event dumpstep (t += T0_/tout_res) {
event dumpstep (t += T0_/tout_res_my) {

  //fprintf(stderr, "I output the restarting files every t=%8E\n", T0_/tout_res), fflush (stderr);
  fprintf(stderr, "I output the restarting files every t=%8E\n", T0_/tout_res_my), fflush (stderr);
  
  if(counter < counter_max) {
    counter++;
  }
  else {
    counter = 1;
  }

  /*
  if (pid () == 0) { // before writing we remove the file that already exists (if any)

    char comm_1[80];
    sprintf (comm_1, "rm -f dump_%d.bin", counter);
    system(comm_1);

  }
  MPI_Barrier(MPI_COMM_WORLD);
  */

  char dname[100];
  u_water.x.nodump = true;
  u_water.y.nodump = true;
  u_water.z.nodump = true;
  sprintf (dname, "dump_%d.bin", counter);
  dump (dname);

  /** 
     Add a symbolic link, log restarting info and size of the bin */

  if (pid () == 0) {

    char comm_2[80];
    sprintf (comm_2, "ln -sf dump_%d.bin restart.bin", counter);
    system(comm_2);

    fflush(stderr);
    char name_0[80];
    sprintf (name_0, "dump_%d.bin", counter);
    FILE * fp = fopen(name_0, "r");
    fseek(fp, 0L, SEEK_END);
    int res = ftell(fp);
    fclose(fp);

    fflush(stderr);
    char name_1[80];
    sprintf (name_1,"log_restart.out");
    FILE * log_sim = fopen(name_1,"a");
    fprintf (log_sim, "%8E %8E %8E %8E\n", t, 1.0*i, 1.0*counter, 1.0*res);
    fclose(log_sim);

  }

}

//event dump_backup (t += T0_/tout_rbk) {
event dump_backup (t += T0_/tout_rbk_my) {
//event dump_backup (i += 2) {
 
  //fprintf(stderr, "I output the restarting for backup files every t=%8E\n", T0_/tout_rbk), fflush (stderr);
  fprintf(stderr, "I output the restarting files for backup every t=%8E\n", T0_/tout_rbk_my), fflush (stderr);
  
  char dname[100];
  u_water.x.nodump = true;
  u_water.y.nodump = true;
  u_water.z.nodump = true;
  sprintf (dname, "./restart_bk/dump_%09d.bin", i);
  dump (dname);

}

/*
event dumpforrestart (t += 1) { // There seems to be problem with tiger with dumping too frequently. Otherwise could do t+=0.1
  char dname[100];
  u_water.x.nodump = true;
  u_water.y.nodump = true;
  u_water.z.nodump = true;
  // p.nodump = false;
  sprintf (dname, "restart");
  dump (dname);
}

event dumpstep (t += 0.5) {
  char dname[100];
  u_water.x.nodump = true;
  u_water.y.nodump = true;
  u_water.z.nodump = true;
  pair.nodump = true;
  // p.nodump = false;
  sprintf (dname, "dump%8E", t);
  dump (dname);
}
*/

/** 
   Adaptive function. uemax is tuned to be 0.3Ustar in most cases. We need a more strict criteria for water speed once the waves starts moving, i.e. t >= RELEASETIME. */ 
   
#if TREE
event adapt (i++) {

  if (i == 5) {
    fprintf(stderr, "uemaxRATIO = %8E\n", uemaxRATIO);
  }

  if (t < RELEASETIME) {
    adapt_wavelet ({f,u}, (double[]){femax,uemax,uemax,uemax}, MAXLEVEL);
  }

  if (t >= RELEASETIME) {
    foreach () {
      foreach_dimension ()
	u_water.x[] = u.x[]*f[];
    }
    adapt_wavelet ({f,u,u_water}, (double[]){femax,uemax,uemax,uemax,uwemax,uwemax,uwemax}, MAXLEVEL);
  }

}
#endif
