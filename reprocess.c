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
#include "sandbox/LS_reinit.h"
#include "sandbox/iso3D.h"  

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
int prt_res = 9;                       // printing resolution
double RELEASETIME = 1000.;            // releasetime

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
   We define some auxiliary fields:
   --> the shifted VoFs and the associated variables */

scalar f2s[];
double my_stp_eta_f2s = 0.1; // we will change later

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
  if (argc > 11)
    prt_res = atoi(argv[11]);
  if (argc > 12)
    RELEASETIME = atof(argv[12]);
  
  if (argc < 13) {

    fprintf(ferr, "Lack of command line arguments. Check! Need %d more arguments\n", 13-argc);
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
  fprintf(stderr, " prt_res = %d\n ", prt_res), fflush (stderr);
  fprintf(stderr, " RELEASETIME = %8E\n ", RELEASETIME), fflush (stderr);
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

/** 
   Creation of a 3D matrix. */

void * matrix_new_3d (int nx, int ny, int nz, size_t size) {
  void ** m = qmalloc(nx, void *);  //Define a pointer that points to every x coordinate.
  char * a  = qmalloc(nx*ny*nz*size, char);
  for (int i=0; i<nx; i++) {
    m[i] = a+i*ny*nz*size;
  }
  return m;
}

/** 
   Output a 3d field in a linearly interpolated uniform grid */

void output_3d (char * fname, scalar s, int maxlevel, bool do_linear, bool print_bin) {

  boundary ({s}); // must be kept since we use interpolate_linear
  int nn = (1<<maxlevel);
  double stp = 0.999999*(L0+X0-X0)/(double)nn; // to avoid interpolated point coincides with fine grid boundary

  double ** field = (double **) matrix_new_3d (nn, nn, nn, sizeof(double));
  for (int i = 0; i < nn; i++) {
    double xp = stp*i + X0 + stp/2.;
    for (int j = 0; j < nn; j++) {
      double yp = stp*j + Y0 + stp/2.;
      for (int k = 0; k < nn; k++) {
        double zp = stp*k + Z0 + stp/2.;
	double val = nodata;
        foreach_point (xp, yp, zp, serial) {
          val = do_linear ? interpolate_linear (point, s, xp, yp, zp) : s[];
	}
        field[i][j*nn+k] = val != nodata ? val : 0.0;
      }
    }
  }

  FILE * fpver = fopen (fname,"w");
  if (pid() == 0) { // master
#if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], cube(nn), MPI_DOUBLE, MPI_SUM, 0,
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
            fprintf(fpver, "%.10e", field[i][j*nn+k]);
            fputc('\n', fpver);  // not double quotation""
          }
        }
      }
    }
    fflush (fpver);
    fclose (fpver); // we close at the end
  }
#if _MPI
  else // slave
    MPI_Reduce (field[0], NULL, cube(nn), MPI_DOUBLE, MPI_SUM, 0,
		MPI_COMM_WORLD);
#endif
  matrix_free (field);

}

// ---------- SAFER XY SLICE (z = zp) ----------
void sliceXY_new (char * fname, scalar s, double zp, int maxlevel, bool do_linear, bool print_bin) {
  boundary ({s}); // required if we may interpolate
  const int nn  = (1 << maxlevel);
  const double stp = 0.999999 * L0 / (double) nn;

  // Wrap zp into [Z0, Z0+L0) in case of periodic z
  while (zp <  Z0)      zp += L0;
  while (zp >= Z0+L0)   zp -= L0;

  double **field = (double **) matrix_new (nn, nn, sizeof(double));

  for (int i = 0; i < nn; i++) {
    const double xp = X0 + (i + 0.5)*stp;
    for (int j = 0; j < nn; j++) {
      const double yp = Y0 + (j + 0.5)*stp;

      double val = nodata;
      foreach_point (xp, yp, zp, serial) {
        if (do_linear && is_leaf(cell))
          val = interpolate_linear(point, s, xp, yp, zp);
        else
          val = s[];
      }
      field[i][j] = (val != nodata) ? val : 0.0;
    }
  }

  FILE *fpver = fopen(fname, "w");
  if (!fpver) { perror("fopen sliceXY_new"); }
  if (pid() == 0 && fpver) {
#if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], sq(nn), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
    if (print_bin) {
      for (int i = 0; i < nn; i++)
        for (int j = 0; j < nn; j++)
          fwrite(&field[i][j], sizeof(double), 1, fpver);
    } else {
      for (int i = 0; i < nn; i++)
        for (int j = 0; j < nn; j++)
          fprintf(fpver, "% .8E\n", field[i][j]);
    }
    fflush(fpver);
    fclose(fpver);
  }
#if _MPI
  else {
    MPI_Reduce (field[0], NULL, sq(nn), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  }
#endif

  matrix_free(field);
}


void sliceYZ_new (char * fname, scalar s, double xp, int maxlevel, bool do_linear, bool print_bin) {

  boundary ({s}); // must be kept since we use interpolate_linear
  int nn = (1<<maxlevel);
  
  double stp = 0.999999*(L0+X0-X0)/(double)nn; // to avoid interpolated point coincides with fine grid boundary

  //fprintf(stderr, "I am here 1\n"), fflush (stderr);
  double ** field = (double **) matrix_new (nn, nn, sizeof(double));
  for (int i = 0; i < nn; i++) {
    double yp = stp*i + Y0 + stp/2.;
    for (int j = 0; j < nn; j++) {
      double zp = stp*j + Z0 + stp/2.;
      //field[i][j] = 0.0;
      //for (int k = 0; k < nn; k++) {
        //double xp = stp*k + X0 + stp/2.;
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
    fflush (fpver);
    fclose (fpver); // we close at the end
    //fprintf(stderr, "I am here 4\n"), fflush (stderr);

  }
#if _MPI
  else // slave
    MPI_Reduce (field[0], NULL, sq(nn), MPI_DOUBLE, MPI_SUM, 0,
                MPI_COMM_WORLD);
#endif
  //fprintf(stderr, "I am here 5\n"), fflush (stderr);

  matrix_free (field);

}

void sliceYZ_avgx_masked (char * fname, scalar s, scalar f,
                          int maxlevel, bool do_linear, bool print_bin,
                          int nx_samples, double f_thresh, char * cov_fname)
{
  boundary ({s});
  if (f.i) boundary ({f}); // mask exists

  const int nn = (1 << maxlevel);
  const double stp_y = 0.999999 * L0 / (double) nn;
  const double stp_z = 0.999999 * L0 / (double) nn;

  if (nx_samples <= 0) nx_samples = (1 << maxlevel); // be mindful: res=10 -> 1024
  const double stp_x = L0 / (double) nx_samples;

  // Accumulators: mean = sum(s*w), sumw = sum(w)
  double **mean = (double **) matrix_new (nn, nn, sizeof(double));
  double **sumw = (double **) matrix_new (nn, nn, sizeof(double));
  for (int i = 0; i < nn; i++)
    for (int j = 0; j < nn; j++)
      mean[i][j] = sumw[i][j] = 0.0;

  // Loop over YZ pixels
  for (int i = 0; i < nn; i++) {
    const double yp = Y0 + (i + 0.5) * stp_y;
    for (int j = 0; j < nn; j++) {
      const double zp = Z0 + (j + 0.5) * stp_z;

      double acc = 0.0, wacc = 0.0;

      // Sample along x
      for (int kx = 0; kx < nx_samples; kx++) {
        double xq = X0 + (kx + 0.5) * stp_x;

        double sval = nodata;
        double fval = 1.0; // default: no mask

        foreach_point (xq, yp, zp, serial) {
          // Interpolate s if requested; avoid linear interp on non-leaf cells for robustness
          if (do_linear) {
            if (is_leaf (cell))
              sval = interpolate_linear (point, s, xq, yp, zp);
            else
              sval = s[];
          } else {
            sval = s[];
          }
          // NEVER interpolate VOF (discontinuous); use cell value
          if (f.i) fval = f[];
        }
        if (sval == nodata) continue; // outside coverage on this rank

        // Clamp f to [0,1] defensively
        if (f.i) {
          if (fval < 0.0) fval = 0.0;
          if (fval > 1.0) fval = 1.0;
        }

        // Weight from VOF
        double w = 1.0;
        if (f.i) {
          if (f_thresh < 0.0) {
            // Continuous weight: proportional contribution of mixed cells
            if      (fval <= 0.0) w = 0.0;
            else if (fval >= 1.0) w = 1.0;
            else                  w = fval;
          } else {
            // Hard threshold: exclude interface band
            if      (fval >= 1.0 - f_thresh) w = 1.0; // water
            else if (fval <= f_thresh)       w = 0.0; // air
            else continue;                              // skip interface
          }
        }

        if (w > 0.0) { acc += sval * w; wacc += w; }
      }

      mean[i][j] = acc;   // sum(s*w)
      sumw[i][j] = wacc;  // sum(w)
    }
  }

// --- MPI reduction of accumulators ---
#if _MPI
if (pid() == 0) {
  // Root rank: reduce in-place into mean/sumw
  MPI_Reduce(MPI_IN_PLACE, mean[0], sq(nn), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(MPI_IN_PLACE, sumw[0], sq(nn), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
} else {
  // Non-root ranks: send buffers; recvbuf must be NULL
  MPI_Reduce(mean[0], NULL,    sq(nn), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(sumw[0], NULL,    sq(nn), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}
#endif

if (pid() == 0) {
  // Convert accumulators -> averaged field (in-place on 'mean')
  for (int i = 0; i < nn; i++)
    for (int j = 0; j < nn; j++)
      mean[i][j] = (sumw[i][j] > 0.0) ? (mean[i][j] / sumw[i][j]) : 0.0;

  // Write averaged field
  FILE *fp = fopen (fname, "w");
  if (!fp) perror("fopen avg file");
  else {
    if (print_bin) {
      for (int i = 0; i < nn; i++)
        for (int j = 0; j < nn; j++)
          fwrite (&mean[i][j], sizeof(double), 1, fp);
    } else {
      for (int i = 0; i < nn; i++)
        for (int j = 0; j < nn; j++)
          fprintf (fp, "% .8E\n", mean[i][j]);
    }
    fflush (fp);
    fclose (fp);
  }

  // Optionally write coverage (sum of weights along x)
  if (cov_fname) {
    FILE *fc = fopen (cov_fname, "w");
    if (!fc) perror("fopen coverage file");
    else {
      if (print_bin) {
        for (int i = 0; i < nn; i++)
          for (int j = 0; j < nn; j++)
            fwrite (&sumw[i][j], sizeof(double), 1, fc);
      } else {
        for (int i = 0; i < nn; i++)
          for (int j = 0; j < nn; j++)
            fprintf (fc, "% .8E\n", sumw[i][j]);
      }
      fflush (fc);
      fclose (fc);
    }
  }
}

matrix_free (mean);
matrix_free (sumw);

}

double my_interpolation(Point point, scalar s, double xp, double yp, double zp) {

  //boundary({s}); // we need to avoid this call since we use my_interpolation
  //                  only if point.level > 0 (i.e. only on some processors)
  //                  The code can get stuck forever here (only some processors enter here)
  //                  Manual boundary conditions must be imposed on any scalar field you want to
  //                  pass to my_interpolate

  // Note: we do not risk to have nodata since
  //       we check point.level > 0 before calling "my_interpolation"

#if dimension == 2
  double xpo = (xp - x)/Delta - s.d.x/2.0;
  double ypo = (yp - y)/Delta - s.d.y/2.0;

  int i = sign(xpo), j = sign(ypo);
  xpo = fabs(xpo), ypo = fabs(ypo);

  return ((s[]*(1. - xpo) + s[i]*xpo)*(1. - ypo) +
          (s[0,j]*(1. - xpo) + s[i,j]*xpo)*ypo);
#else
  double xpo = (xp - x)/Delta - s.d.x/2.0;
  double ypo = (yp - y)/Delta - s.d.y/2.0;
  double zpo = (zp - z)/Delta - s.d.z/2.0;

  int i = sign(xpo), j = sign(ypo), k = sign(zpo);
  xpo = fabs(xpo), ypo = fabs(ypo), zpo = fabs(zpo);

  return (((s[]*(1. - xpo) + s[i]*xpo)*(1. - ypo) +
	   (s[0,j]*(1. - xpo) + s[i,j]*xpo)*ypo)*(1. - zpo) +
	   ((s[0,0,k]*(1. - xpo) + s[i,0,k]*xpo)*(1. - ypo) +
	   (s[0,j,k]*(1. - xpo) + s[i,j,k]*xpo)*ypo)*zpo);
#endif

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

void vof_advection2 (scalar f2, double dt2, face vector uf2, int i) {
  
  face vector uf_stored = uf;
  double dt_stored = dt;
  uf = uf2;
  dt = dt2;
  vof_advection ({f2}, i);
  dt = dt_stored;
  uf = uf_stored;

}

/** 
   Function employed to shift the VoF vertically of my_stp_eta_fs */

void shift_vof (scalar f, coord ufv_fs, int i, double my_stp_eta_fs, double L0, scalar fs) {
  
  // We initialize fs with the latest available f
  scalar_clone(fs,f);
  foreach() {
    fs[] = f[];
  }
  fs.nodump = true; // we need to put here
#if TREE
  fs.refine = fs.prolongation = fraction_refine;
  restriction({fs});
#endif

  // We define the shifting velocity
  face vector uuf_fs[];
  foreach_face() {
    uuf_fs.x[] = ufv_fs.x;
  }

  // We shift fs of my_stp_eta_fs
  double pdt_max = 0.50*CFL*L0/(1 << grid->maxdepth); // CFL = 0.8 and 0.50*CFL = 0.4, which is good for VoF
  unsigned int num_min = my_stp_eta_fs/(fabs(ufv_fs.y)*pdt_max);
  double pdt = my_stp_eta_fs/(1.0*(num_min+1));
  double yp = 0;
  int ps_i = 0;
  while (fabs(yp) < my_stp_eta_fs) {
    vof_advection2 (fs, pdt, uuf_fs, i);
    yp += pdt*ufv_fs.y;
    ps_i++;
    fprintf(stderr, "translation of the VoF along y. Tmp pos: %g. Pseudo it: %d\n", yp, ps_i);
  }

}

/**
   We want to compute some quantities at the interface */

void output_int_qtn (char * fname, int istep, int MAXLEVEL, double time, double RELEASETIME, 
		     scalar c, scalar * list, coord stp_eta, double stp_pos) {

  // Preliminary operation!
  for (scalar s in list) {
    boundary ({s}); // must be kept since we use interpolate_linear
  }

  // We first loop over all the interfacial points 
  // and we count them (per processor)
 
  /* 
  int int_pt = 0;
  foreach(serial) { 
    if (interfacial (point, c)) {
      //if (point.level == MAXLEVEL) {
        coord n = interface_normal(point, c), pp;
        double alpha1 = plane_alpha (c[], n);
        plane_area_center(n, alpha1, &pp);
        double xc = x + Delta*pp.x;
        double yc = y + Delta*pp.y + stp_eta;
        double zc = z + Delta*pp.z;
	zc *= (dimension - 2);
	Point point1 = locate (pc.x, pc.y, pc.z);
	if (point1.level > 0) { // best case
	  POINT_VARIABLES;
	  int_pt++;
	}
      //}
    }
  }
  */

  int int_pt = 0;
  foreach(serial) { 
    if (interfacial (point, c)) {
      int_pt++;
    }
  }

  int tot_column = 8+list_len(list);
  double t_mat[int_pt][tot_column];

  for (int j = 0; j < tot_column; j++) {
    for (int i = 0; i < int_pt; i++) {
      t_mat[i][j] = 0;
    }
  }
  fprintf(stderr, "First pass over interfacial points\n");

  scalar pos[];
#if dimension > 2
  coord G_e = {0.,1.,0.}, Z_e = {0.,0.,0.};
#else
  coord G_e = {0.,1.}, Z_e = {0.,0.};
#endif
  position (c, pos, G_e, Z_e);

  foreach() {
    if(pos[] != nodata) {
      pos[] -= stp_pos; // we correct pos[] to account for the shift in c;
    }
  }

  scalar curv[];
  curvature (c, curv);

  /*
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
  boundary ({Sxx,Syy,Szz,Sxy,Sxz,Syz}); // must be kept
  boundary ({p_a}); // must be kept
  boundary ((scalar *){u}); // must be kept 
  */

  /*
  int count = 0;
  foreach(serial) {
    if (interfacial (point, c)) {
      //if (point.level == MAXLEVEL) {
        coord n = interface_normal(point, c), pp;
        double alpha1 = plane_alpha (c[], n);
        plane_area_center(n, alpha1, &pp);
	double eta = pos[];
	double cur = curv[];
	double norm_2 = sqrt(sq(n.x) + sq(n.y) + sq(n.z));
        double xc = x + Delta*pp.x;
        double yc = y + Delta*pp.y + stp_eta;
        double zc = z + Delta*pp.z;
        zc *= (dimension - 2);
	Point point1 = locate (pc.x, pc.y, pc.z);
	if (point1.level > 0) { // best case
	 
	  POINT_VARIABLES;

	  t_mat[count][0]  = xc;
	  t_mat[count][1]  = zc;
	  t_mat[count][2]  = interpolate_linear(point, p_a, pc.x, pc.y, pc.z);
	  t_mat[count][3]  = interpolate_linear(point, Sxx, pc.x, pc.y, pc.z);
	  t_mat[count][4]  = interpolate_linear(point, Syy, pc.x, pc.y, pc.z);
	  t_mat[count][5]  = interpolate_linear(point, Szz, pc.x, pc.y, pc.z);
	  t_mat[count][6]  = interpolate_linear(point, Sxy, pc.x, pc.y, pc.z);
	  t_mat[count][7]  = interpolate_linear(point, Sxz, pc.x, pc.y, pc.z);
	  t_mat[count][8]  = interpolate_linear(point, Syz, pc.x, pc.y, pc.z);
	  t_mat[count][9]  = interpolate_linear(point, u.x, pc.x, pc.y, pc.z);
	  t_mat[count][10] = interpolate_linear(point, u.y, pc.x, pc.y, pc.z);
#if dimension > 2
	  t_mat[count][11] = interpolate_linear(point, u.z, pc.x, pc.y, pc.z);
#else
	  t_mat[count][11] = 0;
#endif

	  t_mat[count][12] = eta;           // already defined at the interface
	  t_mat[count][13] = cur;           // already defined at the interface
	  t_mat[count][14] = n.x/norm_2;    // already defined at the interface
	  t_mat[count][15] = n.y/norm_2;    // already defined at the interface
	  t_mat[count][16] = n.z/norm_2;    // already defined at the interface
	  t_mat[count][17] = Delta;         // no interpolation is meaningful
	  count++;

	}
      //}
    }
  }
  */

  int count = 0;
  foreach(serial) {
    if (interfacial (point, c)) {
      //if (point.level == MAXLEVEL) {
        coord n = interface_normal(point, c), pp;
        double alpha1 = plane_alpha (c[], n);
        plane_area_center(n, alpha1, &pp);
        coord pc = {0.,0.,0.}, o = {x,y,z};
	foreach_dimension() {
          pc.x = o.x + Delta*pp.x + stp_eta.x;
	}
	t_mat[count][0] = pc.x;
	t_mat[count][1] = pc.z;
	int q = 2;
	for (scalar s in list) {
	  t_mat[count][q++] = interpolate_linear(point, s, pc.x, pc.y, pc.z);
	}
	t_mat[count][12] = pos[];  // already defined at the interface
	t_mat[count][13] = curv[]; // already defined at the interface
	t_mat[count][14] = n.x;    // already defined at the interface
	t_mat[count][15] = n.y;    // already defined at the interface
	t_mat[count][16] = n.z;    // already defined at the interface
	t_mat[count][17] = Delta;  // no interpolation is meaningful
	count++;
      //}
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
  if( pid() == 0 ) {
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
    int tot_el_p0[nproc], disp0[nproc];
    for (int i = 0; i < nproc; i++) {
      tot_el_p0[i] = tot_el_p[i];
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
    fprintf (global_int, "%.10e %.10e %09d %09d %09d\n", 
		          time, time-RELEASETIME, istep, tot_int_pt, tot_column);
    fclose(global_int);

  }
  fprintf(stderr, "Print\n");

}

/** 
   Define some variables for the profile output */

void profile_output (vector u, scalar p, scalar p_hd, scalar f,
		     vertex scalar phi, int fsign, double fthr, int set_fthr, 
		     double ymax, double ymin,
		     int istep, int MAXLEVEL, char * name_pf) {

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

  scalar * list_s = {u.x, u.y, u.z, p, p_hd, omega, uv, duy, diss, ke, vke, dkdy, f};
  profiles (list_s, phi, fname = name_pf, ymax, ymin, f, fsign, fthr, set_fthr,
	    rf = 1.0, fname = name_pf, n = 1 << MAXLEVEL);
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
  sprintf (comm, "mkdir -p field_3d");
  system(comm);
  sprintf (comm, "mkdir -p profiles");
  system(comm);
  sprintf (comm, "mkdir -p eta");
  system(comm);
  sprintf (comm, "mkdir -p eta/eta_loc");
  system(comm);

  /**
     Restore the file! */

  if (restore ("restart.bin")) {
    fprintf(stderr, "File restored!\n"), fflush (stderr);
  }
  else {
    fprintf(stderr, "Restarting files absent! Please provide them!\n"), fflush (stderr);
    return 1;
  }
  
  /**
    Specify the surface elevation shift */
  
  my_stp_eta_f2s = 0.1/k_; // following Wu et al. (JFM2022)

}

event end (i = 0; t = t) {

  /**
   * Shift the VOF as requested (positive y-direction here)
   */
  coord ufv_f2s = {0., +1., 0.};
  shift_vof (f, ufv_f2s, i, my_stp_eta_f2s, L0, f2s);

  /**
   * Hydrostatic reference position (phi_1) based on VOF
   */
  scalar phi_1[];
#if dimension > 2
  coord Z_hd = {0.,0.,0.};
#else
  coord Z_hd = {0.,0.};
#endif
  position (f, phi_1, G, Z_hd, add = false);

  /**
   * Quick diagnostics / views
   */
  vector omg[];
  vorticity3D (u, omg);

  scalar u_x[], omgy[], omgz[], omg_mod[];
  foreach() {
    u_x[]   = u.x[]/Ustar;
    omgz[]  = omg.z[]/(2.0*pi/T0_);
    omg_mod[] = 0.;
    foreach_dimension()
      omg_mod[] += sq(omg.x[]/(2.0*pi/T0_));
    omg_mod[] = sqrt(omg_mod[]);
  }
  boundary ({u_x, omgz, omg.x, omg.y, omg.z});
  boundary ({omg_mod});
  fprintf(stderr, "Min omg_mod = %8E\n", statsf(omg_mod).min), fflush(stderr);
  fprintf(stderr, "Max omg_mod = %8E\n", statsf(omg_mod).max), fflush(stderr);

  fprintf(stderr, "Simulation restarts from a dumped file (TWO-PHASE)\n"), fflush(stderr);
  clear();
  view (fov = 27.5, theta = pi/4.0, phi = pi/50.0, psi = 0.0, ty = -0.5,
        width = 1300, height = 1300, bg = {1,1,1}, samples = 4);
  box(false, lc = {1,1,1}, lw = 0.1);
  draw_vof ("f", color = "u.x");
  squares ("u_x" , linear = true, min = -2.0 , max = 20.0, n = {0,0,1}, alpha = -L0/2.0, map = rain_cm);
  squares ("omgz", linear = true, min = -60.0, max = 30.0, n = {1,0,0}, alpha = -L0/2.0, map = balance_cm);
  {
    char name_1[80];
    sprintf (name_1, "3D_%d.ppm", istep);
    static FILE * fp = fopen (name_1, "a");
    save (fp = fp);
    fflush (fp);
  }

  fprintf(stderr, "I computed the dissipation\n"), fflush(stderr);
  clear();
  view (fov = 10.5, camera = "front", ty = -0.2,
        width = 1300, height = 1300, bg = {1,1,1}, samples = 4);
  box(false, lc = {1,1,1}, lw = 0.1);
  isoline2 ("f", val = 0.5, alpha = 1e-4, lc = {1,0,1}, lw = 3);
  squares ("u_x" , linear = true, min = -2.0 , max = 20.0, n = {0,0,1}, alpha = 1e-4, map = rain_cm);
  cells(n={0,0,1}, alpha = 1e-4);
  {
    char name_1[80];
    sprintf (name_1, "zoom_vel_%d.ppm", istep);
    static FILE * fp = fopen (name_1, "a");
    save (fp = fp);
    fflush (fp);
  }

  /**
   * Simple strain/“dissipation” proxy and stresses (if needed later)
   */
  fprintf(stderr, "I computed the dissipation\n"), fflush(stderr);
  scalar diss[];
  scalar Sxx[], Syy[], Szz[], Sxy[], Sxz[], Syz[];
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
    double sqterm = ( sq(SDeformxx) + sq(SDeformxy) + sq(SDeformxz)
                    + sq(SDeformyx) + sq(SDeformyy) + sq(SDeformyz)
                    + sq(SDeformzx) + sq(SDeformzy) + sq(SDeformzz) );
    diss[] = sqterm;

    Sxx[] = 2.*dudx; 
    Syy[] = 2.*dvdy;
    Szz[] = 2.*dwdz;
    Sxy[] = dudy+dvdx;
    Sxz[] = dudz+dwdx;
    Syz[] = dvdz+dwdy;
  }

  /**
   * Field outputs
   */
  char filename[100], covname[256];
  const int res        = prt_res;
  const bool print_bin = true;

  // Update halos before sampling
  boundary ({f});
  boundary ((scalar *){u});

  // --- Thin 2D slices (YZ at x = 0; XY at z = 0) ---
// VOF: never interpolate (robust near interface)
fprintf(stderr, "[thin] fv YZ start\n"), fflush(stderr);
sprintf (filename, "./slices/fv_2d_%09d.bin", istep);
sliceYZ_new (filename, f,   /*xp=*/0., res, /*do_linear=*/false, print_bin);
fprintf(stderr, "[thin] fv YZ done\n"), fflush(stderr);

// Velocities: linear interpolation is fine
fprintf(stderr, "[thin] ux YZ start\n"), fflush(stderr);
sprintf (filename, "./slices/ux_2d_%09d.bin", istep);
sliceYZ_new (filename, u.x, /*xp=*/0., res, /*do_linear=*/true,  print_bin);
fprintf(stderr, "[thin] ux YZ done\n"), fflush(stderr);

fprintf(stderr, "[thin] uy YZ start\n"), fflush(stderr);
sprintf (filename, "./slices/uy_2d_%09d.bin", istep);
sliceYZ_new (filename, u.y, /*xp=*/0., res, /*do_linear=*/true,  print_bin);
fprintf(stderr, "[thin] uy YZ done\n"), fflush(stderr);


fprintf(stderr, "[thin] uz XY start\n"), fflush(stderr);
sprintf (filename, "./slices/uz_2d_%09d.bin", istep);
// prueba: do_linear=false aquí
sliceYZ_new (filename, u.z, /*zp=*/0., res, /*do_linear=*/false,  print_bin);
fprintf(stderr, "[thin] uz XY done\n"), fflush(stderr);


// --- Spanwise averages along x for YZ planes ---
fprintf(stderr, "I output YZ <·>_x\n"), fflush(stderr);
const int NX_SAMPLES = 256; // choose 64–256; avoid 0 (which would default to 2^res)

// VOF averaged along x (no interpolation)
fprintf(stderr, "[avgx] fv start\n"), fflush(stderr);
sprintf(filename, "./slices/fv_avgx_yz_%09d.bin", istep);
sprintf(covname , "./slices/fv_cov_yz_%09d.bin",  istep);
sliceYZ_avgx_masked(filename, /*s=*/f, /*mask=*/f, res,
                    /*do_linear=*/false, /*print_bin=*/true,
                    /*nx_samples=*/NX_SAMPLES, /*f_thresh=*/-1.0,
                    /*cov_fname=*/covname);
fprintf(stderr, "[avgx] fv done\n"), fflush(stderr);

// u.x averaged along x (masked by water)
fprintf(stderr, "[avgx] ux start\n"), fflush(stderr);
sprintf(filename, "./slices/ux_avgx_yz_%09d.bin", istep);
sprintf(covname , "./slices/ux_cov_yz_%09d.bin",  istep);
sliceYZ_avgx_masked(filename, u.x, f, res,
                    /*do_linear=*/true,  /*print_bin=*/true,
                    /*nx_samples=*/NX_SAMPLES, /*f_thresh=*/-1.0,
                    /*cov_fname=*/covname);
fprintf(stderr, "[avgx] ux done\n"), fflush(stderr);

// u.y averaged along x (masked by water)
fprintf(stderr, "[avgx] uy start\n"), fflush(stderr);
sprintf(filename, "./slices/uy_avgx_yz_%09d.bin", istep);
sprintf(covname , "./slices/uy_cov_yz_%09d.bin",  istep);
sliceYZ_avgx_masked(filename, u.y, f, res,
                    /*do_linear=*/true,  /*print_bin=*/true,
                    /*nx_samples=*/NX_SAMPLES, /*f_thresh=*/-1.0,
                    /*cov_fname=*/covname);
fprintf(stderr, "[avgx] uy done\n"), fflush(stderr);

// u.z averaged along x (masked by water)
fprintf(stderr, "[avgx] uz start\n"), fflush(stderr);
sprintf(filename, "./slices/uz_avgx_yz_%09d.bin", istep);
sprintf(covname , "./slices/uz_cov_yz_%09d.bin",  istep);
sliceYZ_avgx_masked(filename, u.z, f, res,
                    /*do_linear=*/true,  /*print_bin=*/true,
                    /*nx_samples=*/NX_SAMPLES, /*f_thresh=*/-1.0,
                    /*cov_fname=*/covname);
fprintf(stderr, "[avgx] uz done\n"), fflush(stderr);

  /**
   * Simple run logs
   */
  if (pid() == 0) {
    char name_1[80];
    sprintf (name_1,"./profiles/log_pro.out");
    FILE * log_sim = fopen(name_1,"a");
    fprintf (log_sim, "%.10e %.10e %d\n", t, t-RELEASETIME, istep);
    fclose(log_sim);
  }

  if (pid() == 0) { 
    char name_1[80];
    sprintf (name_1,"info_restart.out");
    FILE * log_sim = fopen(name_1, "a");
    fprintf (log_sim, "%.10e %.10e %d\n", t, t-RELEASETIME, istep);
    fclose(log_sim);
  }

  // Return 1 to stop after this postprocess event (single-shot run)
  return 1;
}