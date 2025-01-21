/** 
This is the wave wind interaction simulation with logarithmic wind profile and adaptive grid. This is a turbulent configuration. */

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"  // reduced gravity
#include "view.h"
#include "sandbox/colormap.h"
#include "sandbox/maxruntime.h"  
#include "sandbox/profile6.h"  

/** 
   Input parameters of the simulations. 
   (we preinitialize them to some values, which will be overwritten by
    the ones passed on the command line) */

double d_istep = 1;                    // double istep
double Re_ast = 720;                   // Friction Reynolds number
double BO = 200;                       // Bond number
double Re_wave = 2.55e4;               // Wave Reynolds number
double UstarRATIO = 0.5;               // Friction velocity over wave speed
double r_L0lam = 4.0;                  // ratio between the box size and one wavelength
double rho_r = 1.225/1000.0;           // reference density (air)
double mu_r  = 2.2471881948940954e-06; // reference dynamic viscosity (air)
int MAXLEVEL = 10;                     // max level of the simulation
int MINLEVEL = 6;                      // min level of the simulation

/**
   We define these values: the wave number, fluid depth, wave period, gravity acceleration,
   the wave speed, the friction velocity, the wavelength, alter_MU. */

int istep = 1;    // we change later
double k_  = 1.0; // we change later 
double h_  = 1.0; // we change later 
double T0_ = 1.0; // we change later
double g_  = 1.0; // we change later
double c_  = 1.0; // we change later
double Ustar = 1.0; // we change later
double lam_ = 1.0; // we change later

/**
   The density and viscosity ratios are those of air and water. 
   note: the ratio is air properties over water properties */

double RHO_RATIO = 1.225/1000.;
double MU_RATIO  = 18.31e-6/10.0e-4;

/**
   Input paramters are passed in from the command line. */

int main(int argc, char *argv[]) {

  /**
     Provide the inputs */

  maxruntime (&argc, argv);
  
  if (argc > 1)
    d_istep = atof (argv[1]);
  if (argc > 2)
    Re_ast  = atof(argv[2]);
  if (argc > 3)
    BO = atof(argv[3]);
  if (argc > 4)
    Re_wave = atof(argv[4]);
  if (argc > 5)
    UstarRATIO = atof(argv[5]);
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
  
  if (argc < 11) {

    fprintf(ferr, "Lack of command line arguments. Check! Need %d more arguments\n", 11-argc);
    return 1;

  }

  fprintf(stderr, "************************\n"), fflush (stderr);
  fprintf(stderr, "maximum runtime = %.10e seconds\n", _maxruntime), fflush (stderr);
  fprintf(stderr, "Check of input parameters only\n"), fflush (stderr);
  fprintf(stderr, " istep = %.10e\n ", d_istep), fflush (stderr);
  fprintf(stderr, " Re_ast = %.10e\n ", Re_ast), fflush (stderr);
  fprintf(stderr, " Bo = %.10e\n ", BO), fflush (stderr);
  fprintf(stderr, " Re_wave = %.10e\n ", Re_wave), fflush (stderr);
  fprintf(stderr, " UstarRATIO = %.10e\n ", UstarRATIO), fflush (stderr);
  fprintf(stderr, " r_L0lam = %.10e\n ", r_L0lam), fflush (stderr);
  fprintf(stderr, " Reference density = %8E\n ", rho_r), fflush (stderr);
  fprintf(stderr, " Reference dynamic viscosity = %8E\n ", mu_r), fflush (stderr);
  fprintf(stderr, " MAX LEVEL = %d\n ", MAXLEVEL), fflush (stderr);
  fprintf(stderr, " MIN LEVEL = %d\n ", MINLEVEL), fflush (stderr);
  fprintf(stderr, "************************\n"), fflush (stderr);

  fprintf(stderr, "double istep = %09E\n", d_istep), fflush (stderr);
  istep = (int) d_istep;
  fprintf(stderr, "istep = %09d\n", istep), fflush (stderr);

  /**
     The domain is a cubic box centered on the origin and of length
     $L_0=2*\pi$, periodic in the x- and z-directions. */

  L0   = 2*pi;
  h_   = L0/(2.0*pi); // Water depth set L0/(2*pi) following Wu, Popinet & Deike JFM2022
  lam_ = L0/r_L0lam;
  k_   = 2.0*pi/lam_;
  
  origin (-L0/2., 0, -L0/2.);
  size(L0);

  /**
     Set parameters. */

  rho2     = rho_r;
  mu2      = mu_r;
  rho1     = rho2/RHO_RATIO;
  Ustar    = (mu2*Re_ast)/(rho2*(L0-h_));
  c_       = Ustar/UstarRATIO; 
  g_       = k_*sq(c_)/( 1.0+(rho1-rho2)/(rho1*BO) ); // tune g_ based on c_ and the Bond number 
  f.sigma  = g_*(rho1-rho2)/(BO*sq(k_));
  G.y      = -g_;
  mu1      = rho1*c_*(2.*pi/k_)/Re_wave;
  T0_      = 2.*pi/sqrt(g_*k_ + f.sigma*k_*k_*k_/rho1);

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
     We run! */

  run();

}

void vorticity3D (vector u, vector omg) {

  foreach() {
    omg.x[] = ( (u.z[0,1,0]-u.z[0,-1,0]) - (u.y[0,0,1]-u.y[0,0,-1]) )/(2.0*Delta);
    omg.y[] = ( (u.x[0,0,1]-u.x[0,0,-1]) - (u.z[1,0,0]-u.z[-1,0,0]) )/(2.0*Delta);
    omg.z[] = ( (u.y[1,0,0]-u.y[-1,0,0]) - (u.x[0,1,0]-u.x[0,-1,0]) )/(2.0*Delta);
  }
}

void sliceXY_new (char * fname, scalar s, double zp, int maxlevel, bool do_linear, bool print_bin) {

  boundary ({s}); // must be kept since we use interpolate_linear
  int nn = (1<<maxlevel);
  double stp = 0.999999*(L0+X0-X0)/(double)nn; // to avoid interpolated point coincides with fine grid boundary

  //fprintf(stderr, "I am here 1\n"), fflush (stderr);
  double ** field = (double **) matrix_new (nn, nn, sizeof(double));
  for (int i = 0; i < nn; i++) {
    double xp = stp*i + X0 + stp/2.;
    for (int j = 0; j < nn; j++) {
      double yp = stp*j + Y0 + stp/2.;
      field[i][j] = 0.0;
      for (int k = 0; k < nn; k++) {
        double zp = stp*k + Z0 + stp/2.;
        double val = nodata;
        foreach_point (xp, yp, zp, serial) {
          val = do_linear ? interpolate_linear (point, s, xp, yp, zp) : s[];
	}
	field[i][j] = val != nodata ? val : 0.0; 
	/*
	if (val != nodata) {
          field[i][j] += val/(double)nn;
	}
	*/
      }
    }
  }

  //fprintf(stderr, "I am here 2\n"), fflush (stderr);
  FILE * fpver = fopen (fname,"w");
  if (pid() == 0) { // master

  // reduce it to the first pid()
#if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], sq(nn), MPI_DOUBLE, MPI_SUM, 0,
		MPI_COMM_WORLD);
#endif

    // print
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
    //fprintf(stderr, "I am here 4\n"), fflush (stderr);

  }
#if _MPI
  else // slave
    MPI_Reduce (field[0], NULL, sq(nn), MPI_DOUBLE, MPI_SUM, 0,
		MPI_COMM_WORLD);
#endif
  //fprintf(stderr, "I am here 5\n"), fflush (stderr);

  matrix_free (field);
  fclose (fpver); // we close at the end

}


void output_2d_span_avg (char * fname, scalar s, int maxlevel, bool do_linear, bool print_bin) {

  boundary ({s}); // must be kept since we use interpolate_linear
  int nn = (1<<maxlevel);
  double stp = 0.999999*(L0+X0-X0)/(double)nn; // to avoid interpolated point coincides with fine grid boundary

  //fprintf(stderr, "I am here 1\n"), fflush (stderr);
  double ** field = (double **) matrix_new (nn, nn, sizeof(double));
  for (int i = 0; i < nn; i++) {
    double xp = stp*i + X0 + stp/2.;
    for (int j = 0; j < nn; j++) {
      double yp = stp*j + Y0 + stp/2.;
      field[i][j] = 0.0;
      for (int k = 0; k < nn; k++) {
        double zp = stp*k + Z0 + stp/2.;
        double val = nodata;
        foreach_point (xp, yp, zp, serial) {
          val = do_linear ? interpolate_linear (point, s, xp, yp, zp) : s[];
	}
	field[i][j] += val != nodata ? val/(double)nn : 0.0; 
	/*
	if (val != nodata) {
          field[i][j] += val/(double)nn;
	}
	*/
      }
    }
  }

  //fprintf(stderr, "I am here 2\n"), fflush (stderr);
  FILE * fpver = fopen (fname,"w");
  if (pid() == 0) { // master

  // reduce it to the first pid()
#if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], sq(nn), MPI_DOUBLE, MPI_SUM, 0,
		MPI_COMM_WORLD);
#endif

    // print
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
    //fprintf(stderr, "I am here 4\n"), fflush (stderr);

  }
#if _MPI
  else // slave
    MPI_Reduce (field[0], NULL, sq(nn), MPI_DOUBLE, MPI_SUM, 0,
		MPI_COMM_WORLD);
#endif
  //fprintf(stderr, "I am here 5\n"), fflush (stderr);

  matrix_free (field);
  fclose (fpver); // we close at the end

}

void profile_output (vector u,  
		     vertex scalar phi, int istep, int MAXLEVEL, char * name_pf) {

  scalar omega = new scalar;
  vorticity (u, omega);

  scalar uv = new scalar; scalar duy = new scalar;
  scalar diss = new scalar; scalar ke = new scalar;
  scalar vke = new scalar; scalar dkdy = new scalar;
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
    vke[]  = u.y[]*ke[];
  }
  foreach() {
    dkdy[] = (ke[0,1]-ke[0,-1])/(2.*Delta);
  }

  profiles ({u.x, u.y, u.z, omega, uv, duy, diss, ke, vke, dkdy, f}, phi, rf = 1.0, fname = name_pf, n = 1 << MAXLEVEL);
  delete({omega,uv,duy,diss,ke,vke,dkdy});

}

event init (i = 0) {

  /**
     Print relevant simulation parameters */

  fprintf(stderr, "************************\n"), fflush (stderr);
  fprintf(stderr, "A-posteriori check of simulation parameters\n"), fflush (stderr);
  fprintf(stderr, " Re_ast = %.10e\n ", rho2*Ustar*(L0-h_)/mu2), fflush (stderr);
  fprintf(stderr, " Bo = %.10e\n ", g_*(rho1-rho2)/(f.sigma*sq(k_))), fflush (stderr);
  fprintf(stderr, " Re_wave = %.10e\n ", rho1*c_*(2*pi/k_)/mu1), fflush (stderr);
  fprintf(stderr, " UstarRATIO = %.10e\n ", Ustar/c_), fflush (stderr);
  fprintf(stderr, " r_L0lam = %.10e\n ", r_L0lam), fflush (stderr);
  fprintf(stderr, " Reference density = %8E\n ", rho_r), fflush (stderr);
  fprintf(stderr, " Reference dynamic viscosity = %8E\n ", mu_r), fflush (stderr);
  fprintf(stderr, " MAX LEVEL = %d\n ", MAXLEVEL), fflush (stderr);
  fprintf(stderr, " MIN LEVEL = %d\n ", MINLEVEL), fflush (stderr);
  fprintf(stderr, "************************\n"), fflush (stderr);
  
  /**
     Create directories! */

  fprintf(stderr, "Create directories if needed\n"), fflush (stderr);

  char comm[80];
  sprintf (comm, "mkdir -p slices");
  system(comm);
  sprintf (comm, "mkdir -p field");
  system(comm);
  sprintf (comm, "mkdir -p profiles");
  system(comm);

  /**
     Restart from the latest dump files or initialize a new simulation */

  if (restore ("restart.bin")) {
    vector omg[];
    vorticity3D (u, omg);
    scalar u_x[], omgy[], omgz[], omg_mod[];
    foreach() {
      u_x[]  = u.x[]/Ustar;
      omgz[] = omg.z[]/(2.0*pi/T0_);
      omg_mod[] = 0;
      foreach_dimension() {
        omg_mod[] += sq(omg.x[]/(2.0*pi/T0_));
      }
      omg_mod[] = sqrt(omg_mod[]);
    }
    boundary ({u_x,omgz,omg.x,omg.y,omg.z}); // must be kept?
    boundary ({omg_mod}); // must be kept?
    fprintf(stderr, "Min omg_mod = %8E\n", statsf(omg_mod).min), fflush (stderr);
    fprintf(stderr, "Max omg_mod = %8E\n", statsf(omg_mod).max), fflush (stderr);
    
    fprintf(stderr, "Simulation restarts from a dumped file (TWO-PHASE)\n"), fflush (stderr);
    clear();
    view (fov = 27.5,
          theta = pi/4.0, phi = pi/50.0, psi = 0.0, ty = -0.5,
          width = 1300, height = 1300, bg = {1,1,1}, samples = 4);
    box(false, lc = {1,1,1}, lw = 0.1);
    draw_vof ("f", color = "u.x");
    squares ("u_x" , linear = true, min = -2.0 , max = 20.0, n = {0,0,1}, alpha = -L0/2.0, map = rain_cm);
    squares ("omgz", linear = true, min = -60.0, max = 30.0, n = {1,0,0}, alpha = -L0/2.0, map = balance_cm);
    //squares ("omg_mod", linear = true, min = 0.0, max = 100.0, n = {1,0,0}, alpha = -L0/2.0, map = balance_cm);
    char name_1[80];
    sprintf (name_1, "3D_%d.ppm", istep);
    {
      static FILE * fp = fopen (name_1, "a");
      save (fp = fp);
      fflush (fp);
    }

    fprintf(stderr, "I computed the dissipation\n"), fflush (stderr);
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
      double sqterm = (sq(SDeformxx) + sq(SDeformxy) + sq(SDeformxz) +
          	       sq(SDeformyx) + sq(SDeformyy) + sq(SDeformyz) +
          	       sq(SDeformzx) + sq(SDeformzy) + sq(SDeformzz));
      diss[] = sqterm;
    }

    char filename[100];
    int res        = 9;
    bool do_linear = true;
    bool print_bin = true;
   
    ///*
    int len = 2;
    double pos[len]; 
    for (int i = 0; i <= len; i++) {
      double stp = L0/4.0;
      pos[i] = -L0/4.0+i*stp;
    }
   
    fprintf(stderr, "I make slices\n"), fflush (stderr);
    for (int s = 0; s <= len; s++) {
      double stp = pos[s];
      sprintf (filename, "./slices/di_2d_%09d_%03d.bin", istep, s); // dissipation
      sliceXY_new (filename, diss, stp, res, do_linear, print_bin);
      sprintf (filename, "./slices/fv_2d_%09d_%03d.bin", istep, s); // volume of fluid
      sliceXY_new (filename, f   , stp, res, do_linear, print_bin);
      sprintf (filename, "./slices/ox_2d_%09d_%03d.bin", istep, s); // vorticity - x
      sliceXY_new (filename, omg.x, stp, res, do_linear, print_bin);
      sprintf (filename, "./slices/oy_2d_%09d_%03d.bin", istep, s); // vorticity - y
      sliceXY_new (filename, omg.y, stp, res, do_linear, print_bin);
      sprintf (filename, "./slices/oz_2d_%09d_%03d.bin", istep, s); // vorticity - z
      sliceXY_new (filename, omg.z, stp, res, do_linear, print_bin);
    }
    //*/
    
    fprintf(stderr, "I make slices spanwise-averaged\n"), fflush (stderr);
    sprintf (filename, "./field/di_2d_%09d.bin", istep); // dissipation
    output_2d_span_avg (filename, diss , res, do_linear, print_bin);
    sprintf (filename, "./field/fv_2d_%09d.bin", istep); // volume of fluid
    output_2d_span_avg (filename, f    , res, do_linear, print_bin);
    sprintf (filename, "./field/ox_2d_%09d.bin", istep); // vorticity - x
    output_2d_span_avg (filename, omg.x, res, do_linear, print_bin);
    sprintf (filename, "./field/oy_2d_%09d.bin", istep); // vorticity - y
    output_2d_span_avg (filename, omg.y, res, do_linear, print_bin);
    sprintf (filename, "./field/oz_2d_%09d.bin", istep); // vorticity - z
    output_2d_span_avg (filename, omg.z, res, do_linear, print_bin);


    // We compute the vertical coordinate
    vertex scalar phi_v1[];
    foreach_vertex() {
      phi_v1[] = y;
    }

    char file[99];

    fprintf(stderr, "I plot the profiles\n"), fflush (stderr);
    sprintf (file, "./profiles/prof_%09d_v1.out", istep); // only vertical coordinate
    profile_output (u, phi_v1, istep, MAXLEVEL, file);
    fflush (stderr);
   
    /*
    if (pid() == 0) {
    
      char name_1[80];
      sprintf (name_1,"./profiles/log_pro.out");
      FILE * log_sim = fopen(name_1,"a");
      fprintf (log_sim, "%.10e %.10e\n", t, 1.0*istep);
      fclose(log_sim);
    
    }
    */

    //return 1;

  }
  else {
    fprintf(stderr, "Restarting files absent! Please provide them!\n"), fflush (stderr);
    return 1;
  }
   
}

event end (i = 0; t = t) {

  // To print the time information
  if(pid() == 0) { 
    fflush(stderr);
    char name_1[80];
    sprintf (name_1,"info_restart.out");
    FILE * log_sim = fopen(name_1,"a");
    fprintf (log_sim, "%.10e %d\n", t, istep);
    fclose(log_sim);
  }

  return 1;

}
