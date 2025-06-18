/** 
   This is the wave wind interaction simulation with logarithmic wind profile and adaptive grid. 
   This is a turbulent configuration. */

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"  // reduced gravity
#include "view.h"
#include "tag.h"
#include "sandbox/redistance2.h"
#include "sandbox/basic_funcs.h"   
#include "sandbox/profile6.h"   // From Antoon
#include "sandbox/my_funcs.h"   
#include "sandbox/maxruntime.h"  
#include "sandbox/colormap.h"
#include "sandbox/iso3D.h"
//#include "sandbox/perfs.h"   

/** 
   Input parameters of the simulations. 
   (we preinitialize them to some values, which will be overwritten by
    the ones passed on the command line) */

double Re_ast = 720;                   // Friction Reynolds number
double BO = 200;                       // Bond number
double Re_wave = 2.55e4;               // Wave Reynolds number
double UstarRATIO = 0.5;               // Friction velocity over wave speed
double ak = 0.3;                       // initial steepness
double r_L0lam = 4.0;                  // ratio between the box size and one wavelength
double rho_r = 1.225/1000.0;           // reference density (air)
double mu_r  = 2.2471881948940954e-06; // reference dynamic viscosity (air)
int MAXLEVEL = 10;                     // max level of the simulation
int MINLEVEL = 6;                      // min level of the simulation
double RELEASETIME = 1000;             // physical time before which the precursor simulation should be run
double uemaxRATIO = 0.3;               // Ratio to control the maximum error on the air velocity
double end_sim = 1000;                 // end of the simulation (physical time unit)
int do_eta_loc = 1;                    // output or not eta_loc
int do_profile = 1;                    // output or not profiles
int do_fields  = 1;                    // output or not fields
int do_tagging = 1;                    // output or not tagging
int prt_res = 9;                       // printing resolution
int from_pr = 1;                       // from precursor simulation

/**
   We define these values: the wave number, fluid depth, wave period, gravity acceleration,
   the wave speed, the friction velocity, the wavelength, alter_MU. */

double k_  = 1.0; // we change later 
double h_  = 1.0; // we change later 
double T0_ = 1.0; // we change later
double g_  = 1.0; // we change later
double c_  = 1.0; // we change later
double Ustar = 1.0; // we change later
double lam_ = 1.0; // we change later
double alter_MU = 1.0; // we change later

/**
   For some statistics. */

bool is_first = true;
int i_sp = 0;
double t_sp = 0.0;
#include "sandbox/perfs.h"  

/**
   Clean the VoF from bubbles/droplets. */

bool do_remove   = true;   // do it or not
int min_size     = 3;      // size of the bubble/droplet in grid cells
double threshold = 1.0e-4; // VoF and 1.0-VoF threshold
double cirp_th   = 0.20;   // minimum threshold for position
double eta_m0    = 1.0;    // initial mean eta0

/**
   Define some output frequencies (fraction of the wave period, T0) */

double tout_glo_my = 64.0; // output frequency of global observables
double tout_tag_my = 64.0; // output frequency of tagging
double tout_eta_my = 32.0; // output frequency of interfacial quantity
double tout_fld_my = 32.0; // output frequency of fields
double tout_pro_my = 32.0; // output frequency of 1d profile
double tout_mov_my = 4.00; // output frequency of the movie
double tout_res_my = 4.00; // output frequency of restart
double tout_cut_my = 2.00; // output frequency of some 2D cuts (keep <<tout_cut_my>> equal to <<tout_rbk_my>> for future reuse)
double tout_rbk_my = 2.00; // output frequency of the backup restarting files (keep <<tout_rbk_my>> equal to <<tout_cut_my>> for future reuse)

/**
   For the restarting step. */

int counter_max = 2;
int counter = 0;

/**
   The error on the components of the velocity field used for adaptive
   refinement. */

double uemax  = 1.0e-2; // we are going to change later
double uwemax = 1.0e-3; // we are going to change later
double femax  = 1.0e-4; // we keep this value

/**
   The density and viscosity ratios are those of air and water. 
   note: the ratio is air properties over water properties */

double RHO_RATIO = 1.225/1000.;
double MU_RATIO  = 18.31e-6/10.0e-4;

/**
   We define some auxiliary fields:
   --> the shifted VoFs and the associated variables;
   --> the hydrodynamic pressure with the hydrostatic load; 
   --> the water velocity field;
   --> the forcing term. */

//scalar f2s[], f1s[];
scalar f2s[];
double my_stp_eta_f2s = 0.1; // we will change later
//double my_stp_eta_f1s = 0.1; // we will change later
scalar p_hd[];
vector u_water[]; 
face vector av[];

/**
   Input paramters are passed in from the command line. */

int main(int argc, char *argv[]) {

  /**
     Provide the inputs */

  maxruntime (&argc, argv);

  if (argc > 1)
    Re_ast  = atof(argv[1]);
  if (argc > 2)
    BO = atof(argv[2]);
  if (argc > 3)
    Re_wave = atof(argv[3]);
  if (argc > 4)
    UstarRATIO = atof(argv[4]);
  if (argc > 5)
    ak = atof(argv[5]);
  if (argc > 6)
    r_L0lam = atof(argv[6]);
  if (argc > 7)
    rho_r = atof(argv[7]);
  if (argc > 8)
    mu_r = atof(argv[8]);
  if (argc > 9)
    MAXLEVEL = atoi(argv[9]);
  if (argc > 10)
    MINLEVEL = atoi(argv[10]);
  if (argc > 11)
    RELEASETIME = atof(argv[11]);
  if (argc > 12)
    uemaxRATIO = atof(argv[12]);
  if (argc > 13)
    end_sim  = atof(argv[13]);
  if (argc > 14)
    do_eta_loc = atoi(argv[14]); // do_eta_loc = 1 --> true, else --> false
  if (argc > 15)
    do_profile = atoi(argv[15]); // do_profile = 1 --> true, else --> false
  if (argc > 16)
    do_fields  = atoi(argv[16]); // do_fields  = 1 --> true, else --> false
  if (argc > 17)
    do_tagging = atoi(argv[17]); // do_tagging = 1 --> true, else --> false
  if (argc > 18)
    prt_res    = atoi(argv[18]); 
  if (argc > 19)
    from_pr    = atoi(argv[19]);
  
  if (argc < 20) {

    fprintf(ferr, "Lack of command line arguments. Check! Need %d more arguments\n", 20-argc);
    return 1;

  }

  fprintf(stderr, "************************\n"), fflush (stderr);
  fprintf(stderr, "maximum runtime = %.10e seconds\n", _maxruntime), fflush (stderr);
  fprintf(stderr, "Check of input parameters only\n"), fflush (stderr);
  fprintf(stderr, " Re_ast = %.10e\n ", Re_ast), fflush (stderr);
  fprintf(stderr, " Bo = %.10e\n ", BO), fflush (stderr);
  fprintf(stderr, " Re_wave = %.10e\n ", Re_wave), fflush (stderr);
  fprintf(stderr, " UstarRATIO = %.10e\n ", UstarRATIO), fflush (stderr);
  fprintf(stderr, " a_0k = %.10e\n ", ak), fflush (stderr);
  fprintf(stderr, " r_L0lam = %.10e\n ", r_L0lam), fflush (stderr);
  fprintf(stderr, " Reference density = %8E\n ", rho_r), fflush (stderr);
  fprintf(stderr, " Reference dynamic viscosity = %8E\n ", mu_r), fflush (stderr);
  fprintf(stderr, " MAX LEVEL = %d\n ", MAXLEVEL), fflush (stderr);
  fprintf(stderr, " MIN LEVEL = %d\n ", MINLEVEL), fflush (stderr);
  fprintf(stderr, " RELEASETIME = %.10e\n ", RELEASETIME), fflush (stderr);
  fprintf(stderr, " uemaxRatio = %.10e\n ", uemaxRATIO), fflush (stderr);
  fprintf(stderr, " tend = %.10e\n ", end_sim), fflush (stderr);
  fprintf(stderr, " do_eta_loc = %d\n ", do_eta_loc), fflush (stderr);
  fprintf(stderr, " do_profile = %d\n ", do_profile), fflush (stderr);
  fprintf(stderr, " do_fields  = %d\n ", do_fields), fflush (stderr);
  fprintf(stderr, " do_tagging = %d\n ", do_tagging), fflush (stderr);
  fprintf(stderr, " Printing resolution = %d\n ", prt_res), fflush (stderr);
  fprintf(stderr, " Restart from precursor = %d\n ", from_pr), fflush (stderr);
  fprintf(stderr, "************************\n"), fflush (stderr);

  /**
     The domain is a cubic box centered on the origin and of length
     $L_0=2*\pi$, periodic in the x- and z-directions. */

  L0   = 2.0*pi;
  h_   = L0/(2.0*pi); // Water depth set L0/(2*pi) following Wu, Popinet & Deike JFM2022
  lam_ = L0/r_L0lam;
  k_   = 2.0*pi/lam_;
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
     To set the dimensional parameters we can use three strategies:
       1. $u_\ast$ to a fixed value and we modify $c$. This strategy allows to keep the same $Re_\ast$ and (1a) modify either the wave Reynolds number and keep fixed the viscosity ratio or (2a) modify the viscosity ratio and keep fixed the wave Reynolds number (strategy of Wu, Popinet and Deike, JFM2022);
       2. $c$ to a fixed value and we modify $u_\ast$. This strategy allows to keep the same viscosity ratio and wave Reynolds number, but we change $Re_\ast$.
     */

  rho2     = rho_r;
  mu2      = mu_r;
  rho1     = rho2/RHO_RATIO;
  Ustar    = (mu2*Re_ast)/(rho2*(L0-h_));
  c_       = Ustar/UstarRATIO; 
  g_       = k_*sq(c_)/( 1.0+(rho1-rho2)/(rho1*BO) ); // tune g_ based on c_ and the Bond number 
  f.sigma  = g_*(rho1-rho2)/(BO*sq(k_));
  G.y      = -g_;
  mu1      = rho1*c_*(2.*pi/k_)/Re_wave;
  alter_MU = (mu2/mu1)/MU_RATIO;
  T0_      = 2.*pi/sqrt(g_*k_ + f.sigma*k_*k_*k_/rho1);
    
  /**
     Give the address of av[] to a[] so that acceleration can be changed */

  a = av;
  
  /**
     We set the refinement level */

  // Refine according to 
  uemax  = uemaxRATIO*Ustar;
  uwemax = 0.001*c_;
  fprintf (stderr, "RELEASETIME = %.10e, uemax = %.10e, uwemax = %.10e, femax = %.10e\n", 
		   RELEASETIME, uemax, uwemax, femax);
  
  /**
     We run! */

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
  return eta1 + ak*eta2 + sq(ak)*eta3 + h_;

}

/** 
   Specify the initial velocity field in the water
   Note:
   --> the analytical expression of u_x and u_y 
       are consistent with a third-order Stokes wave. 
   --> we added a correction for g_ to account for surface tension effects
       (this correction is negligible for large Bo) */

double u_x (double x, double y) {

  double gt_  = g_ + f.sigma*k_*k_/rho1;
  double alpa = 1./tanh(k_*h_);
  double a_   = ak/k_;
  double sgma = sqrt(gt_*k_*tanh(k_*h_)*
		     (1. + k_*k_*a_*a_*(9./8.*(sq(alpa) - 1.)*
					(sq(alpa) - 1.) + sq(alpa))));
  double A_ = a_*gt_/sgma;
  return A_*cosh(k_*(y + h_))/cosh(k_*h_)*k_*cos(k_*x) +
    ak*3.*ak*A_/(8.*alpa)*(sq(alpa) - 1.)*(sq(alpa) - 1.)*
    cosh(2.0*k_*(y + h_))*2.*k_*cos(2.0*k_*x)/cosh(2.0*k_*h_) +
    ak*ak*1./64.*(sq(alpa) - 1.)*(sq(alpa) + 3.)*
    (9.*sq(alpa) - 13.)*
    cosh(3.*k_*(y + h_))/cosh(3.*k_*h_)*a_*a_*k_*k_*A_*3.*k_*cos(3.*k_*x);

}

double u_y (double x, double y) {

  double gt_  = g_ + f.sigma*k_*k_/rho1;
  double alpa = 1./tanh(k_*h_);
  double a_   = ak/k_;
  double sgma = sqrt(gt_*k_*tanh(k_*h_)*
		     (1. + k_*k_*a_*a_*(9./8.*(sq(alpa) - 1.)*
					(sq(alpa) - 1.) + sq(alpa))));
  double A_ = a_*gt_/sgma;
  return A_*k_*sinh(k_*(y + h_))/cosh(k_*h_)*sin(k_*x) +
    ak*3.*ak*A_/(8.*alpa)*(sq(alpa) - 1.)*(sq(alpa) - 1.)*
    2.*k_*sinh(2.0*k_*(y + h_))*sin(2.0*k_*x)/cosh(2.0*k_*h_) +
    ak*ak*1./64.*(sq(alpa) - 1.)*(sq(alpa) + 3.)*
    (9.*sq(alpa) - 13.)*
    3.*k_*sinh(3.*k_*(y + h_))/cosh(3.*k_*h_)*a_*a_*k_*k_*A_*sin(3.*k_*x);

}

/** 
   Initial event */

# define POPEN(name, mode) fopen (name ".ppm", mode)

event init (i = 0) {

  /**
     Print relevant simulation parameters */

  fprintf(stderr, "************************\n"), fflush (stderr);
  fprintf(stderr, "A-posteriori check of simulation parameters\n"), fflush (stderr);
  fprintf(stderr, " Re_ast = %.10e\n ", rho2*Ustar*(L0-h_)/mu2), fflush (stderr);
  fprintf(stderr, " Bo = %.10e\n ", g_*(rho1-rho2)/(f.sigma*sq(k_))), fflush (stderr);
  fprintf(stderr, " Re_wave = %.10e\n ", rho1*c_*(2*pi/k_)/mu1), fflush (stderr);
  fprintf(stderr, " UstarRATIO = %.10e\n ", Ustar/c_), fflush (stderr);
  fprintf(stderr, " a_0k = %.10e\n ", ak), fflush (stderr);
  fprintf(stderr, " r_L0lam = %.10e\n ", r_L0lam), fflush (stderr);
  fprintf(stderr, " Reference density = %8E\n ", rho_r), fflush (stderr);
  fprintf(stderr, " Reference dynamic viscosity = %8E\n ", mu_r), fflush (stderr);
  fprintf(stderr, " MAX LEVEL = %d\n ", MAXLEVEL), fflush (stderr);
  fprintf(stderr, " MIN LEVEL = %d\n ", MINLEVEL), fflush (stderr);
  fprintf(stderr, " RELEASETIME = %.10e\n ", RELEASETIME), fflush (stderr);
  fprintf(stderr, " uemaxRATIO = %.10e\n ", uemaxRATIO), fflush (stderr);
  fprintf(stderr, " alter_MU = %.10e\n ", alter_MU), fflush (stderr);
  fprintf(stderr, " end_sim  = %.10e\n ", end_sim), fflush (stderr);
  fprintf(stderr, " do_eta_loc = %d\n ", do_eta_loc), fflush (stderr);
  fprintf(stderr, " do_profile = %d\n ", do_profile), fflush (stderr);
  fprintf(stderr, " do_fields  = %d\n ", do_fields), fflush (stderr);
  fprintf(stderr, " do_tagging = %d\n ", do_tagging), fflush (stderr);
  fprintf(stderr, " gravity = %.10e\n ", g_), fflush (stderr);
  fprintf(stderr, " RHO_RATIO = %.10e\n ", (rho1/rho2)), fflush (stderr);
  fprintf(stderr, " MU_RATIO  = %.10e\n ", (mu1/mu2)), fflush (stderr);
  fprintf(stderr, " Surface tension  = %.10e\n ", f.sigma), fflush (stderr);
  fprintf(stderr, " Friction velocity  = %.10e\n ", Ustar), fflush (stderr);
  fprintf(stderr, " T0 = %.10e\n ", T0_), fflush (stderr);

  fprintf(stderr, " We_w = %.10e\n ", rho1*sq(c_)*(2.0*pi/k_)/f.sigma ), fflush (stderr);
  fprintf(stderr, " Fr_w = %.10e\n ", sq(c_)/(g_*2.0*pi/k_) ), fflush (stderr);
  fprintf(stderr, " We_a = %.10e\n ", rho2*sq(Ustar)*(L0-h_)/f.sigma ), fflush (stderr);
  fprintf(stderr, " Fr_a = %.10e\n ", sq(Ustar)/(g_*(L0-h_)) ), fflush (stderr);

  fprintf(stderr, " Printing resolution = %d\n ", prt_res), fflush (stderr);
  fprintf(stderr, " Restart from precursor = %d\n ", from_pr), fflush (stderr);
  fprintf(stderr, "************************\n"), fflush (stderr);

  /**
     Ensure that u_water is not dumped */

  foreach_dimension() {
    u_water.x.nodump = true;
  }
   
  /**
     Specify the surface elevation shift */
  
  my_stp_eta_f2s = 0.1/k_; // following Wu et al. (JFM2022)
  //my_stp_eta_f1s = 0.1/k_; // following Wu et al. (JFM2022)

  /**
     Create directories for better organization of the output */
  
  fprintf(stderr, "Create directories if needed\n"), fflush (stderr);
  char * dir[] = {"profiles","eta","eta/eta_loc","field","restart_bk",
	          "tagging","tagging/tagging_f1","tagging/tagging_f2",
		  "budgets","slices"};
  for (int i = 0; i < (int)sizeof(dir)/sizeof(dir[0]); i++) {
    char command[256];
    sprintf(command, "mkdir -p %s", dir[i]);
    system(command);
  }

  /**
     Restart from the latest dump files or initialize a new two-phase simulation */
    
  if (from_pr == 1) {

    // Check that the restart_sp.bin exists
    FILE * fp = fopen("restart_sp.bin", "r");
    if( fp == NULL ) {
      fprintf(stderr, "The precursor field to restart the simulation is absent, please provide it!\n"), fflush (stderr);
      return 1; 
    }

    // Restore the simulation
    scalar cs = new scalar; // temporary file just for restore 
    fprintf(stderr, "Scalar created\n"), fflush (stderr);
    restore ("restart_sp.bin");
    delete ({cs}); // deleted since we use f[] from now on
    fprintf(stderr, "Simulation restarts from a dumped file (SINGLE-PHASE)\n"), fflush (stderr);

    // Initalize the VoF profile with density
    fraction (f, WaveProfile(x,z)-y); // we initialize the profile and update the density
    foreach() {
      rhov[] = cm[]*rho(sf[]); 
    }
    
    // Check air velocity
    vector uA_tmp = new vector;
    foreach() {
      foreach_dimension() {
	uA_tmp.x[] = u.x[]*(1.0-f[]);
      }
    }
    stats uxA_tmp = statsf (uA_tmp.x);
    stats uyA_tmp = statsf (uA_tmp.y);
    stats uzA_tmp = statsf (uA_tmp.z);
    delete({uA_tmp.x,uA_tmp.y,uA_tmp.z});
    fprintf(stderr, "Air Max ux = %.10e\n", uxA_tmp.max), fflush (stderr);
    fprintf(stderr, "Air Min ux = %.10e\n", uxA_tmp.min), fflush (stderr);
    fprintf(stderr, "Air Max uy = %.10e\n", uyA_tmp.max), fflush (stderr);
    fprintf(stderr, "Air Min uy = %.10e\n", uyA_tmp.min), fflush (stderr);
    fprintf(stderr, "Air Max uz = %.10e\n", uzA_tmp.max), fflush (stderr);
    fprintf(stderr, "Air Min uz = %.10e\n", uzA_tmp.min), fflush (stderr);

    // Compute some quantities
    vector omg[];
    vorticity3D (u, omg); 

    scalar u_x[], omgz[];
    foreach() {
      u_x[]  = u.x[]/Ustar;
      omgz[] = omg.z[]/(2.0*pi/T0_);
    }
    boundary({u_x,omgz});
    
    // Output an image for checking 
    char stg[80];
    char file[99];
    clear();
    view (fov = 27.5,  
          theta = pi/4.0, phi = pi/50.0, psi = 0.0, ty = -0.5,
          width = 1300, height = 1300, bg = {1,1,1}, samples = 4);
    box(false, lc = {1,1,1}, lw = 0.1);
    draw_vof ("f", color = "u.x");
    squares ("u_x" , linear = true, min = -2.0 , max = 18.0, map = rain_cm   , n = {0,0,1}, alpha = -L0/2.0);
    squares ("omgz", linear = true, min = -60.0, max = 30.0, map = balance_cm, n = {1,0,0}, alpha = -L0/2.0);
    sprintf (stg, "t = %0.3f", 2.0*pi*t/T0_);
    draw_string (stg, size = 30); 
    sprintf (file, "./3D_sp_%09d.ppm", i);
    save(file);

    /* 
    // Test 
    dump ("test_tw.bin");
    return 1;
    */

  }
  else {
    if (restore ("restart.bin")) {
      fprintf(stderr, "Simulation restarts from a dumped file (TWO-PHASE)\n"), fflush (stderr);
    }
    else {
      fprintf(stderr, "Restarting file for two-phase simulation is absent! Please provide it!\n"), fflush (stderr);
      return 1;
    }
  }
   
}

/**
   We log some useful information. */

event log_simulation (i += 10) {
 
  int min_level = +100, max_level = -100;
  foreach(reduction(min:min_level) reduction(max:max_level)) {
    if (level > max_level) max_level = level;
    if (level < min_level) min_level = level;
  }

  if (pid() == 0) { 
    fflush(stderr);
    char name_1[80];
    sprintf (name_1,"log_simulation.out");
    FILE * log_sim = fopen(name_1,"a");
    fprintf (log_sim, "%.10e %.10e %.10e %.10e %d %d\n", t, 1.0*i, t-RELEASETIME, dt, max_level, min_level);
    fclose(log_sim);
  }

  /**
     We continue to log the quantities as in the precusor case
     Remember to copy log_forcing.out of precursor! */

  double u_air_mean = 0.;
  double u_air_diff = 0.;
  double u_wat_mean = 0.;
  double voll_mean  = 0.;
  double volg_mean  = 0.;
  double volg_diff  = 0.;
  foreach(reduction(+:u_air_mean) reduction(+:u_wat_mean)
          reduction(+:volg_mean ) reduction(+:voll_mean )
	  reduction(+:u_air_diff) reduction(+:volg_diff )) {
    u_air_mean += u.x[]*(0.0+f[])*dv();
    u_wat_mean += u.x[]*(1.0-f[])*dv();
    u_air_diff += u.x[]*dv();
    volg_mean  += (0.0+f[])*dv();
    voll_mean  += (1.0-f[])*dv();
    volg_diff  += dv();
  }

  if (pid() == 0) {

    fflush(stderr);
    char name_1[80];
    sprintf (name_1,"log_forcing.out");
    FILE * log_sim = fopen(name_1,"a");
    fprintf (log_sim, "%8E %8E %8E %8E %8E %8E %8E %8E\n", 
		      t, 1.0*i, volg_mean, u_air_mean/volg_mean,
		      voll_mean, u_wat_mean/voll_mean,
		      u_air_diff/volg_diff,volg_diff);
    fclose(log_sim);

  }

}

/**
  Keep the wave profile and set the wave velocity to 0, until the RELEASETIME. */
  
event set_wave(i=0; i++; t<RELEASETIME) {

  fraction (f, WaveProfile(x,z)-y);
  foreach() {
    foreach_dimension() {
      u.x[] = (1.0 - f[])*u.x[];
    }
  }

}

/**
   Release the wave at RELEASETIME. */

event start(t = RELEASETIME) {
  
  fprintf(stderr, "t = RELEASETIME, we initialize the liquid velocity\n"), fflush (stderr);
  fraction (f, WaveProfile(x,z)-y);
  foreach () {
    u.x[] += u_x(x, y-h_)*f[];
    u.y[] += u_y(x, y-h_)*f[];
  }

}

/**
   Forcing term to sustain the flow equivalent to the pressure gradient in x direction. */

event acceleration (i++) {

  double ampl = sq(Ustar)/(L0-h_);
  foreach_face(x) {
    double ff = (f[] + f[-1])/2.;
    av.x[] += ampl*(1.0-ff);
  }

}

event cmpt_f_shifted (t = RELEASETIME; t <= T0_*end_sim; t += T0_/tout_glo_my) {

  // We compute the shifted VoF
  //coord ufv_f1s = {0.,-1.,0.};
  //shift_vof (f, ufv_f1s, i, my_stp_eta_f1s, L0, f1s);
  coord ufv_f2s = {0.,+1.,0.};
  shift_vof (f, ufv_f2s, i, my_stp_eta_f2s, L0, f2s);

  /*
  // We remove debris from f1s
  if(do_remove) {
    remove_droplets (f1s,min_size,threshold,false); // first remove droplets
    remove_droplets (f1s,min_size,threshold,true ); // then remove bubbles
  }
  */
  
  // We remove debris from f2s
  if(do_remove) {
    remove_droplets (f2s,min_size,threshold,false); // first remove droplets
    remove_droplets (f2s,min_size,threshold,true ); // then remove bubbles
  }

  // We estimate the total pressure p_hd with the hydrostatic load
  // Note we need to use f[] to remove the load correctly
  // Note that p_hd is a continous function at the two-phase interface,
  // and, therefore, taking its gradient is less noisy
  scalar phi_1[];
#if dimension > 2
  coord Z_hd = {0.,0.,0.};
#else
  coord Z_hd = {0.,0.};
#endif
  position (f, phi_1, G, Z_hd, add = false);
  scalar_clone(p_hd,p);
  foreach() {
    p_hd[] = p[];
  }
  p_hd.nodump = true; // we need to put here
  foreach() { 
    ///*
    if(interfacial(point,f)) {
      p_hd[] += rhov[]*phi_1[];
    }
    else {
      coord o = {x,y,z};
      foreach_dimension() {
        p_hd[] += rhov[]*G.x*o.x;
      }
    }
    //*/
  }

  /*
  int int_pt1 = 0;
  foreach(reduction(+:int_pt1)) { 
    if (interfacial (point, f)) {
      //if (point.level == MAXLEVEL) {
        coord n = interface_normal(point, f), pp;
        double alpha1 = plane_alpha (f[], n);
        plane_area_center(n, alpha1, &pp);
        double xc = x + Delta*pp.x;
        double yc = y + Delta*pp.y;
        double zc = z + Delta*pp.z;
	zc *= (dimension - 2);
	Point point1 = locate (xc, yc, zc);
	if (point1.level > 0) { // best case
	  POINT_VARIABLES;
	  int_pt1++;
	}
      //}
    }
  }

  int int_pt2 = 0;
  foreach(reduction(+:int_pt2)) { 
    if (interfacial (point, f2s)) {
      //if (point.level == MAXLEVEL) {
        coord n = interface_normal(point, f2s), pp;
        double alpha1 = plane_alpha (f2s[], n);
        plane_area_center(n, alpha1, &pp);
        double xc = x + Delta*pp.x;
        double yc = y + Delta*pp.y;
        double zc = z + Delta*pp.z;
	zc *= (dimension - 2);
	Point point1 = locate (xc, yc, zc);
	if (point1.level > 0) { // best case
	  POINT_VARIABLES;
	  int_pt2++;
	}
      //}
    }
  }

  fprintf(stderr, "# int pt in f[]=   %09d\n", int_pt1); fflush (stderr);
  fprintf(stderr, "# int pt in f2s[]= %09d\n", int_pt2); fflush (stderr);

  // Output an image for checking 
  char stg[80];
  char file[99];
  clear();
  view (fov = 27.5,  
        theta = pi/4.0, phi = pi/50.0, psi = 0.0, ty = -0.5,
        width = 1300, height = 1300, bg = {1,1,1}, samples = 4);
  box(false, lc = {1,1,1}, lw = 0.1);
  draw_vof ("f"  , color = "u.x");
  draw_vof ("f2s", color = "u.x");
  draw_vof ("f1s", color = "u.x");
  sprintf (stg, "t = %0.3f", 2.0*pi*t/T0_);
  draw_string (stg, size = 30); 
  sprintf (file, "./3D_shift_%09d.ppm", i);
  save(file);
  */

}

/** 
   Output video and field in uncompressed format 
   (we can compress later using convert <FILE>.ppm to <FILE>.mpg and open with mplayer) */

event movies (t = RELEASETIME; t <= T0_*end_sim; t += T0_/tout_mov_my) {

  fprintf(stderr, "I output the movies every t=%.10e\n", T0_/tout_mov_my); fflush (stderr);
  
  // Compute some quantities
  scalar u_x[], omgz[];
  vector omg[];
  vorticity3D (u, omg); 

  foreach() {
    u_x[]  = u.x[]/Ustar;
    omgz[] = omg.z[]/(2.0*pi/T0_);
  }
  boundary({u_x,omgz});

  char stg[80];
 
  clear();
  view (fov = 27.5,  
        theta = pi/4.0, phi = pi/50.0, psi = 0.0, ty = -0.5,
	width = 1300, height = 1300, bg = {1,1,1}, samples = 4);
  box(false, lc = {1,1,1}, lw = 0.1);
  draw_vof ("f", color = "u.x");
  squares ("u_x" , linear = true, min = -2.0 , max = 18.0, map = rain_cm   , n = {0,0,1}, alpha = -L0/2.0);
  squares ("omgz", linear = true, min = -60.0, max = 30.0, map = balance_cm, n = {1,0,0}, alpha = -L0/2.0);
  sprintf (stg, "t = %0.3f", 2.0*pi*(t-RELEASETIME)/T0_);
  draw_string (stg, size = 30); 
  {
    static FILE * fp = POPEN ("3D", "a");
    save (fp = fp);
    fflush (fp);
  }

  clear();
  view (fov = 10.5, camera = "front", ty = -0.2,
        width = 1300, height = 1300, bg = {1,1,1}, samples = 4);
  box(false, lc = {1,1,1}, lw = 0.1);
  isoline2 ("f", val = 0.5, alpha = 1e-4, lc = {1,0,0}, lw = 3);
  squares ("u_x" , linear = true, min = -2.0 , max = 20.0, n = {0,0,1}, alpha = 1e-4, map = rain_cm);
  cells(n={0,0,1}, alpha = 1e-4);
  sprintf (stg, "t = %0.3f", 2.0*pi*(t-RELEASETIME)/T0_);
  draw_string (stg, size = 30); 
  {
    static FILE * fp = POPEN ("zoom_int", "a");
    save (fp = fp);
    fflush (fp);
  }

}

event out_pro (t = RELEASETIME; t <= T0_*end_sim; t += T0_/tout_pro_my) {

  if (do_profile == 1) {

    fprintf(stderr, "I output the vertical profiles file every t=%.10e\n", T0_/tout_pro_my);

    // We compute the vertical coordinate
    vertex scalar phi_v1[];
    foreach_vertex() {
      phi_v1[] = y;
    }

    // We compute with the level-set
    /*
    //vertex scalar phi_v2[];
    scalar phi_v2[];
    int imax = 512;
    vof_to_ls (f, phi_v2, imax);
    foreach_vertex() {
      phi_v2[] = fabs(y-1.0) > 0.2 ? y : -phi_v2[]; 
      //phi_v2[] = phi_v2[]; 
    }
    */
    
    char file[99];

    sprintf (file, "./profiles/prof_%09d_v1.out", i); // only vertical coordinate
    profile_output (u, p, p_hd, phi_v1, i, MAXLEVEL, file);
    fflush (stderr);
   
    /* 
    sprintf (file, "./profiles/prof_%09d_v2.out", i); // level-set
    profile_output (u, p, p_hd, phi_v2, i, MAXLEVEL, file);
    fflush (stderr);
    */
  
    if (pid() == 0) {
    
      char name_1[80];
      sprintf (name_1,"./profiles/log_pro.out");
      FILE * log_sim = fopen(name_1,"a");
      fprintf (log_sim, "%.10e %.10e %.10e\n", t, 1.0*i, t-RELEASETIME);
      fclose(log_sim);
    
    }
    //return 1;

  }

}

event fields (t = RELEASETIME; t <= T0_*end_sim; t += T0_/tout_fld_my) {

  if (do_fields == 1) {

    fprintf(stderr, "I output the fields every t=%.10e\n", T0_/tout_fld_my), fflush (stderr);

    int res        = prt_res;
    bool do_linear = true;
    bool print_bin = true;
    
    // Compute some quantities
    scalar uv[], duy[], diss[], ke[], vke[];
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
      double sqterm = 2.*( sq(SDeformxx) + sq(SDeformxy) + sq(SDeformxz) +
          		   sq(SDeformyx) + sq(SDeformyy) + sq(SDeformyz) +
        		   sq(SDeformzx) + sq(SDeformzy) + sq(SDeformzz) );
      uv[]   = u.x[]*u.y[];
      duy[]  = dudy; 
      diss[] = sqterm;
      ke[] = 0.0;
      foreach_dimension() {
        ke[] += sq(u.x[]);
      }
      vke[] = u.y[]*ke[];
    }
    scalar dkdy[];
    foreach() {
      dkdy[] = (ke[0,1]-ke[0,-1])/(2.*Delta);
    }

    /*
    char * nfile[] = {"ux","uy","uz","fv","pr",
                      "uv","du","di","ke","vk","dk"};
    scalar * list_s = {u.x,u.y,u.z,f,p,uv,duy,diss,ke,vke,dkdy};
    if ((int)list_len(list_s) != (int)sizeof(nfile)/sizeof(nfile[0])) {
      fprintf(stderr, "Size of nfile is different from size of list_s. Please check!\n"), fflush (stderr);
      return 1;
    } 
    output_2d_span_avg_ls (nfile, list_s, res, 
		           do_linear, print_bin, i);
    */

    ///*
    char filename[100];
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
    sprintf (filename, "./field/uv_2d_avg_%09d.bin", i); // uv
    output_2d_span_avg (filename,uv   ,res, do_linear, print_bin);
    sprintf (filename, "./field/du_2d_avg_%09d.bin", i); // dudy
    output_2d_span_avg (filename,duy  ,res, do_linear, print_bin);
    sprintf (filename, "./field/ke_2d_avg_%09d.bin", i); // ke
    output_2d_span_avg (filename,ke   ,res, do_linear, print_bin);
    sprintf (filename, "./field/vk_2d_avg_%09d.bin", i); // vke
    output_2d_span_avg (filename,vke  ,res, do_linear, print_bin);
    sprintf (filename, "./field/di_2d_avg_%09d.bin", i); // dissipation
    output_2d_span_avg (filename,diss ,res, do_linear, print_bin);
    sprintf (filename, "./field/dk_2d_avg_%09d.bin", i); // dkdy
    output_2d_span_avg (filename,dkdy ,res, do_linear, print_bin);
    //*/
    
    if (pid() == 0) {
    
      char name_1[80];
      sprintf (name_1,"./field/log_fields.out");
      FILE * log_sim = fopen(name_1,"a");
      fprintf (log_sim, "%.10e %.10e %.10e\n", t, 1.0*i, t-RELEASETIME);
      fclose(log_sim);
    
    }

  }
  //return 1;

}

event slices (t = RELEASETIME; t <= T0_*end_sim; t += T0_/tout_cut_my) {

  if (do_fields == 1) {

    fprintf(stderr, "I output the slices every t=%.10e\n", T0_/tout_cut_my), fflush (stderr);

    int res        = prt_res;
    bool do_linear = true;
    bool print_bin = true;
    
    scalar uv[];
    foreach () {
      uv[] = u.x[]*u.y[];
    }
    
    vector omg[];
    vorticity3D (u, omg); 

    scalar diss[];
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
      double sqterm = 2.*( sq(SDeformxx) + sq(SDeformxy) + sq(SDeformxz) +
          		   sq(SDeformyx) + sq(SDeformyy) + sq(SDeformyz) +
        		   sq(SDeformzx) + sq(SDeformzy) + sq(SDeformzz) );
      diss[] = sqterm;
    }

    /*
    int len = 3;
    for (int s = 0; s < len; s++) {
      char filename[100];
      double stp = -L0/4.0+s*L0/4.0;
      sprintf (filename, "./slices/ux_2d_%03d_%09d.bin", s, i); // x-vel
      sliceXY (filename, u.x, stp, res, do_linear, print_bin);
      sprintf (filename, "./slices/uy_2d_%03d_%09d.bin", s, i); // y-vel
      sliceXY (filename, u.y, stp, res, do_linear, print_bin);
      sprintf (filename, "./slices/uz_2d_%03d_%09d.bin", s, i); // z-vel
      sliceXY (filename, u.z, stp, res, do_linear, print_bin);
      sprintf (filename, "./slices/fv_2d_%03d_%09d.bin", s, i); // VoF
      sliceXY (filename, f  , stp, res, do_linear, print_bin);
      sprintf (filename, "./slices/uv_2d_%03d_%09d.bin", s, i); // uv
      sliceXY (filename, uv , stp, res, do_linear, print_bin);
      sprintf (filename, "./slices/om_2d_%03d_%09d.bin", s, i); // vorticity
      sliceXY (filename, omega, stp, res, do_linear, print_bin);
    }
    */

    ///*
    char * nfile[] = {"ux","uy","uz","fv","uv","ox","oy","oz","di"}; // keep two letters
    scalar * list_s = {u.x,u.y,u.z,f,uv,omg.x,omg.y,omg.z,diss}; 
    if ((int)list_len(list_s) != (int)sizeof(nfile)/sizeof(nfile[0])) {
      fprintf(stderr, "Size of nfile is different from size of list_s. Please check!\n"), fflush (stderr);
      return 1;
    } 
    int len = 3;
    for (int s = 0; s < len; s++) {
      double stp = -L0/4.0+s*L0/4.0;
      sliceXY_ls(nfile, list_s, stp, res, do_linear, print_bin, s, i);
    }
    //*/

    if (pid() == 0) {
    
      char name_1[80];
      sprintf (name_1,"./slices/log_slices.out");
      FILE * log_sim = fopen(name_1,"a");
      fprintf (log_sim, "%.10e %.10e %.10e\n", t, 1.0*i, t-RELEASETIME);
      fclose(log_sim);
    
    }

  }

}

event bulk_budgets (t = RELEASETIME; t <= T0_*end_sim; t += T0_/tout_glo_my) {

  fprintf(stderr, "I output the bulk budgets every t=%.10e\n", T0_/tout_glo_my);

  // --> Kinetic, Potential energy, dissipation of water and air (bulk)
  double keWat = 0., gpeWat = 0.;
  double keAir = 0., gpeAir = 0.;
  foreach(reduction(+:keWat) reduction(+:gpeWat) 
	  reduction(+:keAir) reduction(+:gpeAir)) {
    double norm2_vel = 0.;
    foreach_dimension() {
      norm2_vel += sq(u.x[]);
    }
    keAir  += rho2*norm2_vel*(1.0-f2s[])*dv();
    //keWat  += rho1*norm2_vel*(0.0+f1s[])*dv();
    keWat  += rho1*norm2_vel*(0.0+f[])*dv();
    coord o = {x,y,z};
    double norm2_gpe = 0.;
    foreach_dimension() {
      norm2_gpe += fabs(o.x*G.x);
    }
    gpeAir += rho2*norm2_gpe*(1.0-f2s[])*dv();
    //gpeWat += rho1*norm2_gpe*(0.0+f1s[])*dv();
    //gpeAir += rho2*norm2_gpe*(1.0-f[])*dv();
    gpeWat += rho1*norm2_gpe*(0.0+f[])*dv();
  }
  double rates[2];
  //dissipation_rate(mu1,mu2,u,f1s,f2s,rates);
  dissipation_rate(mu1,mu2,u,f,f2s,rates);
  double dissWater = rates[0];
  double dissAir   = rates[1];
  double gpe_baseWat = 0.0;
  foreach_dimension() {
    //gpe_baseWat += 0.5*sq(h_+0.-my_stp_eta_f1s)*pow(L0,dimension-1)*fabs(G.x);
    gpe_baseWat += 0.5*sq(h_)*pow(L0,dimension-1)*fabs(G.x);
  }
  gpe_baseWat *= rho1;
  
  // --> Mean velocity (all three components) for both phases
  double vol_a_m = 0., vol_w_m = 0.;
  coord vel_a_m = {0.,0.,0.}, vel_w_m = {0.,0.,0.};
  foreach(reduction(+:vol_a_m) reduction(+:vol_w_m) 
	  reduction(+:vel_a_m) reduction(+:vel_w_m)) {
    vol_a_m += (1.0-f2s[])*dv();
    //vol_w_m += (0.0+f1s[])*dv();
    vol_w_m += (0.0+f[])*dv();
    foreach_dimension() {
      vel_a_m.x += u.x[]*(1.0-f2s[])*dv(); 
      //vel_w_m.x += u.x[]*(0.0+f1s[])*dv(); 
      vel_w_m.x += u.x[]*(0.0+f[])*dv(); 
    }
  }

  // --> Mean momentum (all three components)
  scalar phi_tmp[];
  curvature (f, phi_tmp);
  coord mom_m = {0.,0.,0.}, vel_m = {0.,0.,0}, sur_t = {0.,0.,0.};
  foreach(reduction(+:mom_m) reduction(+:vel_m) reduction(+:sur_t)) {
    foreach_dimension() {
      coord int_ft = {0.,0.,0};
      if(interfacial(point,f) && phi_tmp[] != nodata) {
        int_ft.x = phi_tmp[]*(f[1,0,0]-f[-1,0,0])/(2.*Delta);
      }
      mom_m.x += u.x[]*rhov[]*dv();
      vel_m.x += u.x[]*dv();
      sur_t.x += f.sigma*int_ft.x*dv(); 
    }
  }

  // --> compute the viscous stress
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

  // --> Pressure gradient integral, viscous term integral
  coord divc_a_int1 = {0.,0.,0.}, divc_w_int1 = {0.,0.,0.};
  coord grap_a_int1 = {0.,0.,0.}, grap_w_int1 = {0.,0.,0.};
  coord divu_a_int1 = {0.,0.,0.}, divu_w_int1 = {0.,0.,0.};
  foreach(reduction(+:divc_a_int1) reduction(+:divc_w_int1)
          reduction(+:grap_a_int1) reduction(+:grap_w_int1)
	  reduction(+:divu_a_int1) reduction(+:divu_w_int1)) {
    //
    double div_con = 0;
    div_con = (u.x[1,0,0]*u.x[1,0,0]-u.x[-1,0,0]*u.x[-1,0,0])/(2.*Delta) +
	      (u.x[0,1,0]*u.y[0,1,0]-u.x[0,-1,0]*u.y[0,-1,0])/(2.*Delta) +
              (u.x[0,0,1]*u.z[0,0,1]-u.x[0,0,-1]*u.z[0,0,-1])/(2.*Delta);
    divc_a_int1.x += div_con*(1.0-f2s[])*dv();
    //divc_w_int1.x += div_con*(0.0+f1s[])*dv();
    divc_w_int1.x += div_con*(0.0+f[])*dv();
    div_con = (u.y[1,0,0]*u.x[1,0,0]-u.y[-1,0,0]*u.x[-1,0,0])/(2.*Delta) +
	      (u.y[0,1,0]*u.y[0,1,0]-u.y[0,-1,0]*u.y[0,-1,0])/(2.*Delta) +
              (u.y[0,0,1]*u.z[0,0,1]-u.y[0,0,-1]*u.z[0,0,-1])/(2.*Delta);
    divc_a_int1.y += div_con*(1.0-f2s[])*dv();
    //divc_w_int1.y += div_con*(0.0+f1s[])*dv();
    divc_w_int1.y += div_con*(0.0+f[])*dv();
    div_con = (u.z[1,0,0]*u.x[1,0,0]-u.z[-1,0,0]*u.x[-1,0,0])/(2.*Delta) +
	      (u.z[0,1,0]*u.y[0,1,0]-u.z[0,-1,0]*u.y[0,-1,0])/(2.*Delta) +
              (u.z[0,0,1]*u.z[0,0,1]-u.z[0,0,-1]*u.z[0,0,-1])/(2.*Delta);
    divc_a_int1.z += div_con*(1.0-f2s[])*dv();
    //divc_w_int1.z += div_con*(0.0+f1s[])*dv();
    divc_w_int1.z += div_con*(0.0+f[])*dv();
    //
    //coord o = {x,y,z};
    double div_pre = 0;
    div_pre = (p_hd[1,0,0]-p_hd[-1,0,0])/(2.*Delta); // - G.x*o.x*(rhov[1,0,0]-rhov[-1,0,0])/(2.0*Delta); // - rhov[]*G.x;
    grap_a_int1.x += -(div_pre-rho2*G.x)*(1.0-f2s[])*dv();
    //grap_w_int1.x += -(div_pre-rho1*G.x)*(0.0+f1s[])*dv();
    grap_w_int1.x += -(div_pre-rho1*G.x)*(0.0+f[])*dv();
    div_pre = (p_hd[0,1,0]-p_hd[0,-1,0])/(2.*Delta); // - G.y*o.y*(rhov[0,1,0]-rhov[0,-1,0])/(2.0*Delta); // - rhov[]*G.y;
    grap_a_int1.y += -(div_pre-rho2*G.y)*(1.0-f2s[])*dv();
    //grap_w_int1.y += -(div_pre-rho1*G.y)*(0.0+f1s[])*dv();
    grap_w_int1.y += -(div_pre-rho1*G.y)*(0.0+f[])*dv();
    div_pre = (p_hd[0,0,1]-p_hd[0,0,-1])/(2.*Delta); // - G.z*o.z*(rhov[0,0,1]-rhov[0,0,-1])/(2.0*Delta); // - rhov[]*G.z;
    grap_a_int1.z += -(div_pre-rho2*G.z)*(1.0-f2s[])*dv();
    //grap_w_int1.z += -(div_pre-rho1*G.z)*(0.0+f1s[])*dv();
    grap_w_int1.z += -(div_pre-rho1*G.z)*(0.0+f[])*dv();
    //
    double div_tau = 0;
    div_tau = (Sxx[1,0,0]-Sxx[-1,0,0]+Sxy[0,1,0]-Sxy[0,-1,0]+Sxz[0,0,1]-Sxz[0,0,-1])/(2.*Delta);
    divu_a_int1.x += mu2*div_tau*(1.0-f2s[])*dv();
    //divu_w_int1.x += mu1*div_tau*(0.0+f1s[])*dv();
    divu_w_int1.x += mu1*div_tau*(0.0+f[])*dv();
    div_tau = (Sxy[1,0,0]-Sxy[-1,0,0]+Syy[0,1,0]-Syy[0,-1,0]+Syz[0,0,1]-Syz[0,0,-1])/(2.*Delta);
    divu_a_int1.y += mu2*div_tau*(1.0-f2s[])*dv();
    //divu_w_int1.y += mu1*div_tau*(0.0+f1s[])*dv();
    divu_w_int1.y += mu1*div_tau*(0.0+f[])*dv();
    div_tau = (Sxz[1,0,0]-Sxz[-1,0,0]+Syz[0,1,0]-Syz[0,-1,0]+Szz[0,0,1]-Szz[0,0,-1])/(2.*Delta); 
    divu_a_int1.z += mu2*div_tau*(1.0-f2s[])*dv();
    //divu_w_int1.z += mu1*div_tau*(0.0+f1s[])*dv();
    divu_w_int1.z += mu1*div_tau*(0.0+f[])*dv();
    //
  }

  // --> Compute the power vector
  vector t_vel[];
  foreach() {
    //
    double dudx = (u.x[1]     - u.x[-1]    )/(2.*Delta);
    double dudy = (u.x[0,1]   - u.x[0,-1]  )/(2.*Delta);
    double dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
    double dvdx = (u.y[1]     - u.y[-1]    )/(2.*Delta);
    double dvdy = (u.y[0,1]   - u.y[0,-1]  )/(2.*Delta);
    double dvdz = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
    double dwdx = (u.z[1]     - u.z[-1]    )/(2.*Delta);
    double dwdy = (u.z[0,1]   - u.z[0,-1]  )/(2.*Delta);
    double dwdz = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);
    double SDeformxx = dudx + dudx;
    double SDeformxy = dudy + dvdx;
    double SDeformxz = dudz + dwdx;
    double SDeformyx = SDeformxy;
    double SDeformyy = dvdy + dvdy;
    double SDeformyz = dvdz + dwdy;
    double SDeformzx = SDeformxz;
    double SDeformzy = SDeformyz;
    double SDeformzz = dwdz + dwdz;
    //
    t_vel.x[] = SDeformxx*u.x[] + SDeformxy*u.x[] + SDeformxz*u.x[]; 
    t_vel.y[] = SDeformyx*u.y[] + SDeformyy*u.y[] + SDeformyz*u.y[]; 
    t_vel.z[] = SDeformzx*u.z[] + SDeformzy*u.z[] + SDeformzz*u.z[]; 
    //
  }

  // --> Pressure power integral, viscous power integral (per phase and for both kinetic and gravitational forces)
  double grap_pa_int = 0., grap_pw_int = 0.;
  double taup_pa_int = 0., taup_pw_int = 0.;
  foreach(reduction(+:grap_pa_int) reduction(+:grap_pw_int)
	  reduction(+:taup_pa_int) reduction(+:taup_pw_int)) {
    /*
    coord o = {x,y,z};
    double gr_tot = 0;
    foreach_dimension() {
      gr_tot += G.x*o.x;
    }
    */
    foreach_dimension() {
      double dpvel_dl = (u.x[1,0,0]*p_hd[1,0,0]-u.x[-1,0,0]*p_hd[-1,0,0])/(2.0*Delta);
	                //(rhov[1,0,0]-rhov[-1,0,0])/(2.0*Delta)*gr_tot*u.x[] - G.x*u.x[]*rhov[];
      grap_pa_int += -(dpvel_dl-G.x*u.x[]*rho2)*(1.0-f2s[])*dv();
      //grap_pw_int += -(dpvel_dl-G.x*u.x[]*rho1)*(0.0+f1s[])*dv();
      grap_pw_int += -(dpvel_dl-G.x*u.x[]*rho1)*(0.0+f[])*dv();
      double dtvel_dl = (t_vel.x[1,0,0]-t_vel.x[-1,0,0])/(2.0*Delta);
      taup_pa_int += mu2*dtvel_dl*(1.0-f2s[])*dv();
      //taup_pw_int += mu1*dtvel_dl*(0.0+f1s[])*dv();
      taup_pw_int += mu1*dtvel_dl*(0.0+f[])*dv();
    }
  }

  if (pid() == 0) { 

    char name[80];

    fflush(stderr);
    sprintf (name,"./budgets/mom_bud_x_tot.out");
    FILE * log_mtx = fopen(name,"a");
    fprintf (log_mtx, "%8E %8E %8E %8E %8E %8E\n", t, 1.0*i, t-RELEASETIME,
		      vel_m.x, mom_m.x, sur_t.x);
    fclose(log_mtx);

    fflush(stderr);
    sprintf (name,"./budgets/mom_bud_y_tot.out");
    FILE * log_mty = fopen(name,"a");
    fprintf (log_mty, "%8E %8E %8E %8E %8E %8E\n", t, 1.0*i, t-RELEASETIME,
		      vel_m.y, mom_m.y, sur_t.y);
    fclose(log_mty);
    
    fflush(stderr);
    sprintf (name,"./budgets/mom_bud_z_tot.out");
    FILE * log_mtz = fopen(name,"a");
    fprintf (log_mtz, "%8E %8E %8E %8E %8E %8E\n", t, 1.0*i, t-RELEASETIME,
		      vel_m.z, mom_m.z, sur_t.z);
    fclose(log_mtz);

    fflush(stderr);
    sprintf (name,"./budgets/mom_bud_x_air.out");
    FILE * log_max = fopen(name,"a");
    fprintf (log_max, "%8E %8E %8E %8E %8E %8E %8E %8E\n", t, 1.0*i, t-RELEASETIME,
		      vol_a_m, vel_a_m.x, grap_a_int1.x, divu_a_int1.x, divc_a_int1.x);
    fclose(log_max);

    fflush(stderr);
    sprintf (name,"./budgets/mom_bud_y_air.out");
    FILE * log_may = fopen(name,"a");
    fprintf (log_may, "%8E %8E %8E %8E %8E %8E %8E %8E\n", t, 1.0*i, t-RELEASETIME,
		      vol_a_m, vel_a_m.y, grap_a_int1.y, divu_a_int1.y, divc_a_int1.y);
    fclose(log_may);
    
    fflush(stderr);
    sprintf (name,"./budgets/mom_bud_z_air.out");
    FILE * log_maz = fopen(name,"a");
    fprintf (log_maz, "%8E %8E %8E %8E %8E %8E %8E %8E\n", t, 1.0*i, t-RELEASETIME,
		      vol_a_m, vel_a_m.z, grap_a_int1.z, divu_a_int1.z, divc_a_int1.z);
    fclose(log_maz);

    fflush(stderr);
    sprintf (name,"./budgets/mom_bud_x_water.out");
    FILE * log_mwx = fopen(name,"a");
    fprintf (log_mwx, "%8E %8E %8E %8E %8E %8E %8E %8E\n", t, 1.0*i, t-RELEASETIME,
		      vol_w_m, vel_w_m.x, grap_w_int1.x, divu_w_int1.x, divc_w_int1.x);
    fclose(log_mwx);

    fflush(stderr);
    sprintf (name,"./budgets/mom_bud_y_water.out");
    FILE * log_mwy = fopen(name,"a");
    fprintf (log_mwy, "%8E %8E %8E %8E %8E %8E %8E %8E\n", t, 1.0*i, t-RELEASETIME,
		      vol_w_m, vel_w_m.y, grap_w_int1.y, divu_w_int1.y, divc_w_int1.y);
    fclose(log_mwy);
    
    fflush(stderr);
    sprintf (name,"./budgets/mom_bud_z_water.out");
    FILE * log_mwz = fopen(name,"a");
    fprintf (log_mwz, "%8E %8E %8E %8E %8E %8E %8E %8E\n", t, 1.0*i, t-RELEASETIME,
		      vol_w_m, vel_w_m.z, grap_w_int1.z, divu_w_int1.z, divc_w_int1.z);
    fclose(log_mwz);

    fflush(stderr);
    sprintf (name,"./budgets/ke_bud_air.out");
    FILE * fpair = fopen(name,"a");
    fprintf (fpair, "%8E %8E %8E %8E %8E %8E %8E %8E\n", t, 1.0*i, t-RELEASETIME,
		    keAir/2., gpeAir, grap_pa_int, taup_pa_int, dissAir);
    fclose(fpair);

    fflush(stderr);
    sprintf (name,"./budgets/ke_bud_water.out");
    FILE * fpwater = fopen(name,"a");
    fprintf (fpwater, "%8E %8E %8E %8E %8E %8E %8E %8E\n", t, 1.0*i, t-RELEASETIME,
		      keWat/2., gpeWat-gpe_baseWat, grap_pw_int, taup_pw_int, dissWater); 
    fclose(fpwater);
    
  }
  //return 1;

}

event tagging (t = RELEASETIME; t <= T0_*end_sim; t += T0_/tout_tag_my) {

  /*
  scalar f2[];
  scalar_clone (f2, f);
  foreach() {
    f2[] = f[];
  }
#if TREE
  f2.refine = f2.prolongation = fraction_refine;
  restriction({f2});
#endif

  if(do_remove) {
    remove_droplets (f2,min_size,threshold,false); // first remove droplets
    remove_droplets (f2,min_size,threshold,true);  // then remove bubbles
  }
  */

  fprintf(stderr, "I output the tagging every t=%.10e\n", T0_/tout_tag_my);

  if(do_tagging == 1) {

    char name_1[100], name_2[100], name_3[100];
    
    fprintf(stderr, "I do f1-tagging every t=%.10e\n", T0_/tout_glo_my);
    sprintf (name_1, "./tagging/num_bubbles_droplets_f1.out");
    sprintf (name_2, "./tagging/tagging_f1/droplets_f1_%09d.out", i);
    sprintf (name_3, "./tagging/tagging_f1/bubbles_f1_%09d.out", i);
    countDropsBubble (name_1, name_2, name_3, t, RELEASETIME, i, f);
 
    /*   
    fprintf(stderr, "I do f2-tagging every t=%.10e\n", T0_/tout_glo_my);
    sprintf (name_1, "./tagging/num_bubbles_droplets_f2.out");
    sprintf (name_2, "./tagging/tagging_f2/droplets_f2_%09d.out", i);
    sprintf (name_3, "./tagging/tagging_f2/bubbles_f2_%09d.out", i);
    countDropsBubble (name_1, name_2, name_3, t, RELEASETIME, i, f2);
    */

  }
  //return 1;
  
}

event global_obs (t = RELEASETIME; t <= T0_*end_sim; t += T0_/tout_glo_my) {

  scalar f2[];
  scalar_clone (f2, f);
  foreach() {
    f2[] = f[];
  }
#if TREE
  f2.refine = f2.prolongation = fraction_refine;
  restriction({f2});
#endif

  if(do_remove) {
    remove_droplets (f2,min_size,threshold,false); // first remove droplets
    remove_droplets (f2,min_size,threshold,true);  // then remove bubbles
  }

  fprintf(stderr, "I output the global observables file every t=%.10e\n", T0_/tout_glo_my);

  char name_0[100];

  coord stp_eta_0 = {0.,0.,0.}; // it corresponds to 0*Delta 
  double stp_pos_0 = 0.0; // it corresponds to 0*Delta 
  fprintf(stderr, "stp_eta0 for gl. obs. = %.10e\n", stp_pos_0);
  sprintf (name_0, "./budgets/global_obs_ptot_0.out");
  output_global_obs_1 (name_0, i, MAXLEVEL, t, RELEASETIME, eta_m0, cirp_th, k_,
		       f2 , p, stp_eta_0, stp_pos_0);

  coord stp_eta_1 = {0.,0.,0.}; // it corresponds to 0*Delta 
  double stp_pos_1 = my_stp_eta_f2s; // it corresponds to 4*Delta 
  fprintf(stderr, "stp_eta1 for gl. obs. = %.10e\n", stp_pos_1);
  sprintf (name_0, "./budgets/global_obs_ptot_1.out");
  output_global_obs_1 (name_0, i, MAXLEVEL, t, RELEASETIME, eta_m0, cirp_th, k_, 
		       f2s, p, stp_eta_1, stp_pos_1);

}

event eta_loc (t = RELEASETIME; t <= T0_*end_sim; t += T0_/tout_eta_my) {

  if (do_eta_loc == 1) {

    fprintf(stderr, "I output the eta_loc file every t=%.10e\n", T0_/tout_eta_my);
  
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
    //scalar * list_s = {u.x,u.y,u.z,Sxx,Syy,Szz,Sxy,Sxz,Syz,p};
    scalar * list_s = {u.x,u.y,u.z,Sxx,Syy,Szz,Sxy,Sxz,Syz,p_hd};

    fflush(stderr);
    char eta_out[100];
    sprintf (eta_out, "./eta/eta_loc/eta_loc_t%09d.bin", i);

    ///*
    coord stp_eta = {0.,0.,0.};
    double stp_pos = my_stp_eta_f2s; // it corresponds to 4*Delta on 1024**3
    output_int_qtn (eta_out, i, MAXLEVEL, t, RELEASETIME, f2s, list_s, stp_eta, stp_pos);
    //*/

    /* 
    // To revert back and uncomment point.level == MAXLEVEL
    double stp_eta = my_stp_eta_f2s; // it corresponds to 4*Delta on 1024**3
    double stp_pos = 0.0; // it corresponds to 0*Delta
    output_int_qtn (eta_out, i, MAXLEVEL, t, RELEASETIME, f, u, p, stp_eta, stp_pos);
    */
 
  }

}

/** 
   Dump every tout_res period */

event dumpstep (t = RELEASETIME; t <= T0_*end_sim; t += T0_/tout_res_my) {

  fprintf(stderr, "I output the restarting files every t=%.10e\n", T0_/tout_res_my), fflush (stderr);
  
  if(counter < counter_max) {
    counter++;
  }
  else {
    counter = 1;
  }

  char dname[100];
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
    long int res = ftell(fp);
    fclose(fp);

    fflush(stderr);
    char name_1[80];
    sprintf (name_1,"log_restart.out");
    FILE * log_sim = fopen(name_1,"a");
    fprintf (log_sim, "%.10e %.10e %.10e %.10e %.10e\n", t, 1.0*i, t-RELEASETIME, 1.0*counter, 1.0*res);
    fclose(log_sim);

  }

}

event dump_backup (t = RELEASETIME; t <= T0_*end_sim; t += T0_/tout_rbk_my) {
 
  fprintf(stderr, "I output the restarting files for backup every t=%.10e\n", T0_/tout_rbk_my), fflush (stderr);
  
  char dname[100];
  sprintf (dname, "./restart_bk/dump_%09d.bin", i);
  dump (dname);
  
  /** 
     Log restarting info and size of the bin */

  if (pid () == 0) {

    fflush(stderr);
    FILE * fp = fopen(dname, "r");
    fseek(fp, 0L, SEEK_END);
    long int res = ftell(fp);
    fclose(fp);
    
    fflush(stderr);
    char name_1[80];
    sprintf (name_1,"./restart_bk/log_restart_bk.out");
    FILE * log_sim = fopen(name_1,"a");
    fprintf (log_sim, "%.10e %.10e %.10e %.10e\n", t, 1.0*i, t-RELEASETIME, 1.0*res);
    fclose(log_sim);

  }

}

/**
## End 

   We want to run up end_sim physical time. */

event end (t = end_sim*T0_) {

  dump ("end.bin");

  if ( pid() == 0 ) {

    char comm[80];
    sprintf (comm, "ln -sf end.bin restart.bin");
    system(comm);

  }

  return 1; // exit

}

/** 
   Adaptive function. We need a more strict criteria for water speed once the waves starts moving, i.e. t >= RELEASETIME. */ 
   
#if TREE
event adapt (i++) {

  if (t < RELEASETIME) {
    adapt_wavelet ({f,u}, (double[]){femax,uemax,uemax,uemax}, MAXLEVEL, MINLEVEL);
  }

  if (t >= RELEASETIME) {
    foreach () {
      foreach_dimension ()
	u_water.x[] = u.x[]*f[];
    }
    adapt_wavelet ({f,u,u_water}, (double[]){femax,uemax,uemax,uemax,uwemax,uwemax,uwemax}, MAXLEVEL, MINLEVEL);
  }

}
#endif
