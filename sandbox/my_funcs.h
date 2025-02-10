
/** 
   Auxiliary function employed to shift the VoF of a fixed quantity. */

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
   Three components of the vorticity vector */

void vorticity3D (vector u, vector omg) {

  foreach() {
    omg.x[] = ( (u.z[0,1,0]-u.z[0,-1,0]) - (u.y[0,0,1]-u.y[0,0,-1]) )/(2.0*Delta);
    omg.y[] = ( (u.x[0,0,1]-u.x[0,0,-1]) - (u.z[1,0,0]-u.z[-1,0,0]) )/(2.0*Delta);
    omg.z[] = ( (u.y[1,0,0]-u.y[-1,0,0]) - (u.x[0,1,0]-u.x[0,-1,0]) )/(2.0*Delta);
  }

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
   Function to convert the VoF of a LS. WIP */

//void vof_to_ls (scalar f, vertex scalar ls, int imax) {
void vof_to_ls (scalar f, scalar ls_s, int imax) {

  double CFL_ps    = 0.25;
  double delta_min = L0/(1 << grid->maxdepth);
  double dt_max    = CFL_ps*delta_min; // dt_min = CFL*delta_min/(mod(n)), with mod(n)=1
  //scalar ls_s[];
  foreach() {
    //ls_s[] = -(2.*f[] - 1.)*delta_min*0.75;
    ls_s[] = (2.*f[] - 1.)*delta_min*0.75;
  }
#if TREE
  restriction({ls_s});
#endif
  redistance (ls_s, imax, dt_max);
  /*
  foreach_vertex () {
#if dimension > 2
    ls[] = (ls_s[0] + ls_s[-1] + ls_s[0,-1] + ls_s[-1,-1] +
	    ls_s[0,0,-1] + ls_s[-1,0,-1] + ls_s[0,-1,-1] + ls_s[-1,-1,-1])/8.;
#else
    ls[] = (ls_s[] + ls_s[-1] + ls_s[0,-1] + ls_s[-1,-1])/4.;
#endif
  }
  */

}

/** 
   Function to compute the mean value of a generic field s 
   at a specific location y = yp. WIP */

double cmpt_savg_yp (scalar s, double yp, int maxlevel, bool do_linear) {

  int nn = (1<<maxlevel);
  double stp = 0.999999*(L0+X0-X0)/(double)nn; // to avoid interpolated point coincides with fine grid boundary

  double val = 0;
  for (int i = 0; i < nn; i++) {
    double xp = stp*i + X0 + stp/2.;
    for (int k = 0; k < nn; k++) {
      double zp = stp*k + Z0 + stp/2.;
      val += interpolate(s,xp,yp,zp,do_linear)*sq(stp);
    }
  }
  val /= sq(L0);

  return val;

}

/** 
   Output a 2D field at a fixed zp, defined by the user. */

void sliceXY (char * fname, scalar s, double zp, int maxlevel, bool do_linear, bool print_bin) {

  boundary ({s}); // must be kept since we use interpolate_linear
  int nn = (1<<maxlevel);
  double stp = 0.999999*(L0+X0-X0)/(double)nn; // to avoid interpolated point coincides with fine grid boundary

  //fprintf(stderr, "I am here 1\n"), fflush (stderr);
  double ** field = (double **) matrix_new (nn, nn, sizeof(double));
  for (int i = 0; i < nn; i++) {
    double xp = stp*i + X0 + stp/2.;
    for (int j = 0; j < nn; j++) {
      double yp = stp*j + Y0 + stp/2.;
      double val = nodata;
      foreach_point (xp, yp, zp, serial) {
        val = do_linear ? interpolate_linear (point, s, xp, yp, zp) : s[];
      }
      field[i][j] = val != nodata ? val : 0.0; 
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

/** 
   Output a 2D field at a fixed zp, defined by the user (with list) */

void sliceXY_ls (char **fname, scalar * list, double zp, int maxlevel, 
		 bool do_linear, bool print_bin, int pos, int istep) {

  // Preliminary operation!
  for (scalar s in list) {
    boundary ({s}); // must be kept since we use interpolate_linear
  }
  int len = list_len(list);
  int nn = (1<<maxlevel);
  double stp = 0.999999*(L0+X0-X0)/(double)nn; // to avoid interpolated point coincides with fine grid boundary
  
  // Fill the matrix with the entire scalar list
  double ** field = (double **) matrix_new (nn, nn, len*sizeof(double));
  for (int i = 0; i < nn; i++) {
    double xp = stp*i + X0 + stp/2.;
    for (int j = 0; j < nn; j++) {
      double yp = stp*j + Y0 + stp/2.;
      double val = nodata;
      int q = 0;
      for (scalar s in list) {
        foreach_point (xp, yp, zp, serial) {
          val = do_linear ? interpolate_linear (point, s, xp, yp, zp) : s[];
        }
        field[i][len*j+q++] = val != nodata ? val : 0.0; 
      }
    }  
  }

  // Reduction operation
  if (pid() == 0) { // master
#if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], len*sq(nn), MPI_DOUBLE, MPI_SUM, 0,
  	        MPI_COMM_WORLD);
#endif
  }
#if _MPI
  else { // slave
    MPI_Reduce (field[0], NULL, len*sq(nn), MPI_DOUBLE, MPI_SUM, 0,
	        MPI_COMM_WORLD);
  }
#endif

  // Print
  if (pid() == 0) { // only the master prints!
    for (int k = 0; k < len; k++) {
      char filename[100];
      sprintf (filename, "./slices/%s_2d_%03d_%09d.bin", fname[k], pos, istep);
      FILE * fpver = fopen (filename,"w");
      for (int i = 0; i < nn; i++) {
        for (int j = 0; j < nn; j++) {
          fwrite ( &field[i][len*j+k], sizeof(double), 1, fpver );
        }
      }
      fflush (fpver);
      fclose (fpver); // we close at the end
    }
  }
  matrix_free (field); // we deallocate here field

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
   Perform the spanwise averaging in a memory efficient way (with list) */

void output_2d_span_avg_ls (char **fname, scalar * list, int maxlevel, 
		            bool do_linear, bool print_bin, int istep) {

  // Preliminary operation!
  for (scalar s in list) {
    boundary ({s}); // must be kept since we use interpolate_linear
  }
  int len = list_len(list);
  int nn = (1<<maxlevel);
  double stp = 0.999999*(L0+X0-X0)/(double)nn; // to avoid interpolated point coincides with fine grid boundary
  
  // Fill the matrix with the entire scalar list
  double ** field = (double **) matrix_new (nn, nn, len*sizeof(double));
  for (int i = 0; i < nn; i++) {
    double xp = stp*i + X0 + stp/2.;
    for (int j = 0; j < nn; j++) {
      double yp = stp*j + Y0 + stp/2.;
      field[i][j] = 0.0;
      for (int k = 0; k < nn; k++) {
        double zp = stp*k + Z0 + stp/2.;
        double val = nodata;
        int q = 0;
        for (scalar s in list) {
          foreach_point (xp, yp, zp, serial) {
            val = do_linear ? interpolate_linear (point, s, xp, yp, zp) : s[];
          }
          field[i][len*j+q++] += val != nodata ? val/(double)nn : 0.0;
        }	
      }
    }  
  }

  // Reduction operation
  if (pid() == 0) { // master
#if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], len*sq(nn), MPI_DOUBLE, MPI_SUM, 0,
  	        MPI_COMM_WORLD);
#endif
  }
#if _MPI
  else { // slave
    MPI_Reduce (field[0], NULL, len*sq(nn), MPI_DOUBLE, MPI_SUM, 0,
	        MPI_COMM_WORLD);
  }
#endif

  // Print
  if (pid() == 0) { // only the master prints!
    for (int k = 0; k < len; k++) {
      char filename[100];
      sprintf (filename, "./fields/%s_2d_avg_%09d.bin", fname[k], istep);
      FILE * fpver = fopen (filename,"w");
      for (int i = 0; i < nn; i++) {
        for (int j = 0; j < nn; j++) {
          fwrite ( &field[i][len*j+k], sizeof(double), 1, fpver );
        }
      }
      fflush (fpver);
      fclose (fpver); // we close at the end
    }
  }
  matrix_free (field); // we deallocate here field

}

/** 
   Perform the spanwise averaging in a memory efficient way */

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
    fflush(fpver);
  }
#if _MPI
  else // slave
    MPI_Reduce (field[0], NULL, cube(nn), MPI_DOUBLE, MPI_SUM, 0,
		MPI_COMM_WORLD);
#endif
  matrix_free (field);
  fclose (fpver); // we close at the end

}

/** 
   Helper function for fast data sorting */

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
   Define some variables for the profile output */

void profile_output (vector u, scalar p, scalar p_hd,
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

  profiles ({u.x, u.y, u.z, p, p_hd, omega, uv, duy, diss, ke, vke, dkdy, f}, phi, rf = 1.0, fname = name_pf, n = 1 << MAXLEVEL);
  delete({omega,uv,duy,diss,ke,vke,dkdy});

}

/**
   We also want to count the drops and bubbles in the flow */

void countDropsBubble (char * name_1, char * name_2, char * name_3, 
		       double time, double RELEASETIME, int istep, scalar c) {

  scalar m1 = new scalar; // droplets
  scalar m2 = new scalar; // bubbles
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
  double b1x[n1], b1y[n1], b1z[n1];  // droplet
  double v2[n2]; // bubble
  double b2x[n2], b2y[n2], b2z[n2];  // bubble

  /**
  We initialize to zero */

  for (int j=0; j<n1; j++) {
    v1[j] = b1x[j] = b1y[j] = b1z[j] = 0.0;
  }
  for (int j=0; j<n2; j++) {
    v2[j] = b2x[j] = b2y[j] = b2z[j] = 0.0;
  }

  /**
  We proceed with calculation. */

  fprintf(stderr, "I enter the for loop for reduction\n");
  foreach(reduction(+:v1[:n1]),reduction(+:v2[:n2]),
	  reduction(+:b1x[:n1]),reduction(+:b1y[:n1]),reduction(+:b1z[:n1])
	  reduction(+:b2x[:n2]),reduction(+:b2y[:n2]),reduction(+:b2z[:n2])) {

    // droplets
    if (m1[] > 0) {
      int j = m1[] - 1;
      v1[j]  += dv()*c[]; 
      b1x[j] += dv()*c[]*x; 
      b1y[j] += dv()*c[]*y; 
      b1z[j] += dv()*c[]*z; 
    }

    // bubbles
    if (m2[] > 0) {
      int j = m2[] - 1;
      v2[j]  += dv()*(1.0-c[]);
      b2x[j] += dv()*(1.0-c[])*x; 
      b2y[j] += dv()*(1.0-c[])*y; 
      b2z[j] += dv()*(1.0-c[])*z; 
    }

  }
  delete({m1,m2});
  fprintf(stderr, "I go out of the for loop for reduction\n");

  if (pid() == 0) {

    fflush(stderr);

    /**
    We first output the number of droplets and bubbles. */

    FILE * ftot = fopen(name_1,"a");
    fprintf (ftot, "%.10e %.10e %.10e %.10e %.10e\n", time, 1.0*istep, time-RELEASETIME, 1.0*(n1-1), 1.0*(n2-1)); // we remove the main region
    fclose(ftot);

    /**
    We output separately the volume and position of each droplet and bubble to file. */

    FILE * fdrop = fopen(name_2,"a");
    for (int j=0; j<n1; j++) {
      fprintf (fdrop, "%d %.10e %.10e %.10e %.10e\n",
                         j, v1[j], b1x[j]/v1[j], b1y[j]/v1[j], b1z[j]/v1[j]);
    }
    fclose(fdrop);
    fflush(fdrop);

    FILE * fbubb = fopen(name_3,"a");
    for (int j=0; j<n2; j++) {
      fprintf (fbubb, "%d %.10e %.10e %.10e %.10e\n",
                        j, v2[j], b2x[j]/v2[j], b2y[j]/v2[j], b2z[j]/v2[j]);
    }
    fclose(fbubb);
    fflush(fbubb);

  }

}

/** 
   Compute the bulk dissipation in air and water */

int dissipation_rate (double mu1, double mu2, vector u, scalar f1s, scalar f2s, double* rates) {
	
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
			     sq(SDeformzx) + sq(SDeformzy) + sq(SDeformzz));
    rateWater += mu1*(0.0+f1s[])*sqterm; // water
    rateAir   += mu2*(1.0-f2s[])*sqterm; // air
  }
  rates[0] = rateWater;
  rates[1] = rateAir;
  return 0;

}

/**
   We want to compute some quantities at the interface */

void output_int_qtn (char * fname, int istep, int MAXLEVEL, double time, double RELEASETIME, 
		     scalar c, vector u, scalar p_a, double stp_eta, double stp_pos) {

  // We first loop over all the interfacial points 
  // and we count them (per processor)
  
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
	Point point1 = locate (xc, yc, zc);
	if (point1.level > 0) { // best case
	  POINT_VARIABLES;
	  int_pt++;
	}
      //}
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
	Point point1 = locate (xc, yc, zc);
	if (point1.level > 0) { // best case
	 
	  POINT_VARIABLES;

	  t_mat[count][0]  = xc;
	  t_mat[count][1]  = zc;
	  t_mat[count][2]  = my_interpolation(point, p_a, xc, yc, zc);
	  t_mat[count][3]  = my_interpolation(point, Sxx, xc, yc, zc);
	  t_mat[count][4]  = my_interpolation(point, Syy, xc, yc, zc);
	  t_mat[count][5]  = my_interpolation(point, Szz, xc, yc, zc);
	  t_mat[count][6]  = my_interpolation(point, Sxy, xc, yc, zc);
	  t_mat[count][7]  = my_interpolation(point, Sxz, xc, yc, zc);
	  t_mat[count][8]  = my_interpolation(point, Syz, xc, yc, zc);
	  t_mat[count][9]  = my_interpolation(point, u.x, xc, yc, zc);
	  t_mat[count][10] = my_interpolation(point, u.y, xc, yc, zc);
#if dimension > 2
	  t_mat[count][11] = my_interpolation(point, u.z, xc, yc, zc);
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
   We want to compute some global observables */

void output_global_obs_1 (char * fname, int istep, int MAXLEVEL, double time, double RELEASETIME,
		          double eta_m0, double cirp_th, double k_,  
		          scalar c, scalar p_a, double stp_eta, double stp_pos) {

  // Preliminary calculations
  
  // --> position
  scalar pos[];
#if dimension > 2
  coord G = {0.,1.,0.}, Z = {0.,0.,0.};
#else
  coord G = {0.,1.}, Z = {0.,0.};
#endif
  position (c, pos, G, Z);

  foreach() {
    if(pos[] != nodata) {
      pos[] -= stp_pos; // we correct pos[] to account for the shift in c;
    }
  }

  // --> quantities to be probed at the interface
  // n.b.: for the stress we can use just u since it contains derivatives
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

  // --> my diagnosis of eta_my, area_my, amp_my (without offsets)
  double area_my = 0; double eta_my = 0; 
  foreach(reduction(+:area_my), reduction(+:eta_my)) {
    if (interfacial (point, c)) {
      //if (point.level == MAXLEVEL) {
      if (point.level == MAXLEVEL && abs(pos[]-eta_m0) < cirp_th) {
      //if (abs(pos[]-eta_m0) < cirp_th) {
        coord n      = interface_normal(point, c), pp;
        double alpha = plane_alpha (c[], n);
        double ar    = pow(Delta, dimension - 1)*plane_area_center (n, alpha, &pp);
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
      //if (abs(pos[]-eta_m0) < cirp_th) {
        coord n = interface_normal(point, c), pp;
        double alpha = plane_alpha (c[], n);
        double ar = pow(Delta, dimension - 1)*plane_area_center (n, alpha, &pp);
        amp2_my += ar*2.0*sq(pos[]-eta_my);
      }
    }
  }
  amp2_my = amp2_my/area_my;

  // --> compute some mean quantities to be used later (with offset)
  double area = 0; double eta_m = 0; double pr_m = 0;
  double u_xi = 0; double u_yi  = 0; double u_zi = 0;
  foreach(reduction(+:area), reduction(+:eta_m), reduction(+:pr_m),
	  reduction(+:u_xi), reduction(+:u_yi) , reduction(+:u_zi)) {
    if (interfacial (point, c)) {
      //if (point.level == MAXLEVEL) {
      if (point.level == MAXLEVEL && abs(pos[]-eta_m0) < cirp_th) {
      //if (abs(pos[]-eta_m0) < cirp_th) {
	coord n      = interface_normal(point, c), pp;
        double alpha = plane_alpha (c[], n);
        double ar    = pow(Delta, dimension - 1)*plane_area_center (n, alpha, &pp);
	double eta   = pos[]; // must be here since here it is defined
        double xc = x + Delta*pp.x;
        double yc = y + Delta*pp.y + stp_eta;
        double zc = z + Delta*pp.z;
        zc *= (dimension - 2);
	Point point1 = locate (xc, yc, zc);
	if (point1.level > 0) {
        
	  POINT_VARIABLES;
	
	  area  += ar;
	  eta_m += ar*eta; // already defined at the interface

          pr_m += ar*my_interpolation(point, p_a, xc, yc, zc);
          u_xi += ar*my_interpolation(point, u.x, xc, yc, zc);
          u_yi += ar*my_interpolation(point, u.y, xc, yc, zc);
#if dimension > 2
          u_zi += ar*my_interpolation(point, u.z, xc, yc, zc);
#else
	  u_zi = 0.;
#endif

	}
      }
    }
  }
  area  += 1.0e-12; // to prevent arithmetic failure if area = 0
  eta_m  = eta_m/area; // we need the mean surface elevation
  pr_m   = pr_m/area;  // we need the mean pressure
  u_xi   = u_xi/area; u_yi = u_yi/area; u_zi = u_zi/area;
  //fprintf(stderr, "area %8E, eta %8E, pr_m %8E\n", area, eta_m, pr_m);
  //fprintf(stderr, "u.x %8E, u.y %8E, u.z %8E\n", u_xi, u_yi, u_zi);

  // --> compute the amplitude, stress, velocity at the interface, energy fluxes
  double amp2  = 0.0;
  double mf_px = 0.0; double mf_py = 0.0; double mf_pz = 0.0; // mom. flux pressure
  double mp_px = 0.0; double mp_py = 0.0; double mp_pz = 0.0; // mean pressure contribution
  double mf_vx = 0.0; double mf_vy = 0.0; double mf_vz = 0.0; // mom. flux viscous stress
  double ef_pn = 0.0; double ef_pp = 0.0; // en. flux pressure (without-with subtraction)
  double ef_fn = 0.0; double ef_fp = 0.0; // en. flux pressure (no mean pressure)
  double ef_vn = 0.0; double ef_vp = 0.0; // en. flux viscous stress
  foreach(reduction(+:amp2), // amplitude 
	  reduction(+:mf_px), reduction(+:mf_py), reduction(+:mf_pz), // pressure - momentum flux (x,y,z)
	  reduction(+:mp_px), reduction(+:mp_py), reduction(+:mp_pz), // pressure - momentum flux (x,y,z) - (mean pressure)
	  reduction(+:mf_vx), reduction(+:mf_vy), reduction(+:mf_vz), // viscous  - momentum flux (x,y,z)
	  reduction(+:ef_pn), reduction(+:ef_pp), // pressure en flux (without-with subtraction)
	  reduction(+:ef_fn), reduction(+:ef_fp), // pressure en flux (without the mean pressure)
	  reduction(+:ef_vn), reduction(+:ef_vp)) { // viscous en flux (without-with subtraction)
    if (interfacial (point, c)) {
      //if (point.level == MAXLEVEL) {
      if (point.level == MAXLEVEL && abs(pos[]-eta_m0) < cirp_th) {
      //if (abs(pos[]-eta_m0) < cirp_th) {
	coord n      = interface_normal(point, c), pp;
        double alpha = plane_alpha (c[], n);
        double ar    = pow(Delta, dimension - 1)*plane_area_center (n, alpha, &pp);
	double eta   = pos[]; // must be here since here it is defined
	normalize (&n); // |n| = 1, we should normalize since we use the normals for online calculations
        double xc = x + Delta*pp.x;
        double yc = y + Delta*pp.y + stp_eta;
        double zc = z + Delta*pp.z;
        zc *= (dimension - 2);
	Point point1 = locate (xc, yc, zc);
	if (point1.level > 0) {
        
	  POINT_VARIABLES;
	
	  amp2 += ar*2.0*sq(eta-eta_m); // amplitude
 
          double pr_int = 0;   
	  double Sxx_int = 0; double Syy_int = 0; double Szz_int = 0;
	  double Sxy_int = 0; double Sxz_int = 0; double Syz_int = 0;
	  double ux_int  = 0; double uy_int  = 0; double uz_int  = 0;

          pr_int  = my_interpolation(point, p_a, xc, yc, zc);
	  Sxx_int = my_interpolation(point, Sxx, xc, yc, zc);
	  Syy_int = my_interpolation(point, Syy, xc, yc, zc);
	  Szz_int = my_interpolation(point, Szz, xc, yc, zc);
	  Sxy_int = my_interpolation(point, Sxy, xc, yc, zc);
	  Sxz_int = my_interpolation(point, Sxz, xc, yc, zc);
	  Syz_int = my_interpolation(point, Syz, xc, yc, zc);
	  ux_int  = my_interpolation(point, u.x, xc, yc, zc);
	  uy_int  = my_interpolation(point, u.y, xc, yc, zc);
#if dimension > 2
          uz_int  = my_interpolation(point, u.z, xc, yc, zc);
#else
	  uz_int  = 0.;
#endif

          // momentum flux --> pressure	  
	  mf_px += ar*( -(pr_int)*n.x ); // x-Mom. flux due to pressure work
	  mf_py += ar*( -(pr_int)*n.y ); // y-Mom. flux due to pressure work
	  mf_pz += ar*( -(pr_int)*n.z ); // z-Mom. flux due to pressure work
	  
	  mp_px += ar*( +(pr_m)*n.x ); // mean pressure contribution - x
	  mp_py += ar*( +(pr_m)*n.y ); // mean pressure contribution - y
	  mp_pz += ar*( +(pr_m)*n.z ); // mean pressure contribution - z
	  
	  // momentum flux --> viscous stress
	  mf_vx += ar*mu2*(n.x*Sxx_int + n.y*Sxy_int + n.z*Sxz_int); // x-Mom. flux due to viscous dissipation
          mf_vy += ar*mu2*(n.x*Sxy_int + n.y*Syy_int + n.z*Syz_int); // y-Mom. flux due to viscous dissipation
          mf_vz += ar*mu2*(n.x*Sxz_int + n.y*Syz_int + n.z*Szz_int); // z-Mom. flux due to viscous dissipation
          
	  // energy flux --> pressure
	  ef_pn += ar*( -(pr_int-pr_m)*( (ux_int)*n.x+(uy_int)*n.y+(uz_int)*n.z) ); // En. flux due to pressure (no vel. subtraction)
	  ef_pp += ar*( -(pr_int-pr_m)*( (ux_int-u_xi)*n.x+(uy_int-u_yi)*n.y+(uz_int-u_zi)*n.z) ); // En. flux due to pressure
	  ef_fn += ar*( -(pr_int)*( (ux_int)*n.x+(uy_int)*n.y+(uz_int)*n.z) ); // En. flux due to pressure (no vel. subtraction)
	  ef_fp += ar*( -(pr_int)*( (ux_int-u_xi)*n.x+(uy_int-u_yi)*n.y+(uz_int-u_zi)*n.z) ); // En. flux due to pressure

	  // energy flux --> viscous stress
	  ef_vn += ar*mu2*( ( Sxx_int*n.x + Sxy_int*n.y + Sxz_int*n.z )*(ux_int) +
	                    ( Sxy_int*n.x + Syy_int*n.y + Syz_int*n.z )*(uy_int) +
	                    ( Sxz_int*n.x + Syz_int*n.y + Szz_int*n.z )*(uz_int) ); // En. flux due to viscous dissipation (no vel. subtraction)
          ef_vp += ar*mu2*( ( Sxx_int*n.x + Sxy_int*n.y + Sxz_int*n.z )*(ux_int-u_xi) +
	                    ( Sxy_int*n.x + Syy_int*n.y + Syz_int*n.z )*(uy_int-u_yi) +
	                    ( Sxz_int*n.x + Syz_int*n.y + Szz_int*n.z )*(uz_int-u_zi) ); // En. flux due to viscous dissipation
	  
	}
      }
    }
  }
  amp2 = amp2/area;

  // --> compute quantities at the top/bottom boundary (for y-mom in air and water)
  /*
  double pre_a_m = 0.0, pre_w_m = 0.0;
  double vol_a_m = 0.0, vol_w_m = 0.0;
  foreach(reduction(+:pre_a_m) reduction(+:pre_w_m)
          reduction(+:vol_a_m) reduction(+:vol_w_m)) {
    pre_a_m = pre_a_m + (1.0-f[])*p[]*dv();
    pre_w_m = pre_w_m + (0.0+f[])*p[]*dv();
    vol_a_m = vol_a_m + (1.0-f[])*dv();
    vol_w_m = vol_w_m + (0.0+f[])*dv();
  }
  pre_a_m /= vol_a_m; pre_w_m /= vol_w_m;
  */

  double pre_a_m = 0.0;
  foreach_boundary(top,reduction(+:pre_a_m)) {
    double ff = 0.5*(f[0,1,0]+f[]);
    double pf = 0.5*(p[0,1,0]+p[]);
    pre_a_m = pre_a_m + (1.0-ff)*pf*sq(Delta);
  }
  pre_a_m = pre_a_m/sq(L0);

  double pre_w_m = 0.0;
  foreach_boundary(bottom,reduction(+:pre_w_m)) {
    double ff = 0.5*(f[0,-1,0]+f[]);
    double pf = 0.5*(p[0,-1,0]+p[]);
    pre_w_m = pre_w_m + (0.0+ff)*pf*sq(Delta);
  }
  pre_w_m = pre_w_m/sq(L0);

  double py_a_flux = 0.0, tx_a_flux = 0.0, ty_a_flux = 0.0;
  foreach_boundary (top,reduction(+:py_a_flux),reduction(+:tx_a_flux)
		    reduction(+:ty_a_flux)) {
    double pre_int = 0.5*(p[0,1,0]+p[]); // we interpolate to the face
    py_a_flux += sq(Delta)*( -(pre_int) );  
    tx_a_flux += sq(Delta)*mu2*(u.x[0,1,0]-u.x[])/Delta;  
    ty_a_flux += sq(Delta)*mu2*(u.y[0,1,0]-u.y[])/Delta;  
  }
  double py_w_flux = 0.0, ty_w_flux = 0.0;
  foreach_boundary (bottom,reduction(+:py_w_flux),reduction(+:ty_w_flux)) {
    double pre_int = 0.5*(p[0,-1,0]+p[]); // we interpolate to the face
    py_w_flux += sq(Delta)*( -(pre_int) );  
    ty_w_flux += sq(Delta)*mu1*(u.x[]-u.x[0,-1,0])/Delta;  
  }

  if (pid() == 0) {
  
    fflush(stderr);
    FILE * global_obs = fopen (fname, "a");
    fprintf (global_obs, "%8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E\n", 
                         time, 1.0*istep, time-RELEASETIME, area, eta_m, k_*sqrt(amp2), area_my, eta_my, k_*sqrt(amp2_my), pr_m,  
			 mf_px, mf_py, mf_pz, 
			 mp_px, mp_py, mp_pz, 
			 mf_vx, mf_vy, mf_vz, 
			 u_xi, u_yi, u_zi, 
			 ef_pn, ef_pp,
			 ef_fn, ef_fp,
			 ef_vn, ef_vp,
			 f.sigma*area, 
			 pre_a_m, pre_w_m, 
			 py_a_flux,tx_a_flux,ty_a_flux,py_w_flux,ty_w_flux);
    fclose(global_obs);

  }

}
