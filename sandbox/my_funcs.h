/**
# My functions */

void vof_advection2 (scalar f2, double dt2, face vector uf2, int i) {
  face vector uf_stored = uf;
  double dt_stored = dt;
  uf = uf2;
  dt = dt2;
  vof_advection ({f2}, i);
  dt = dt_stored;
  uf = uf_stored;
}

void sliceXY(char * fname, scalar s, double zp, int maxlevel, bool do_linear) {

  FILE * fpver = fopen (fname,"w");
  int nn = (1<<maxlevel);
  double ** field = matrix_new (nn, nn, sizeof(double));
  double stp = L0/(double)nn;
  for (int i = 0; i < nn; i++) {
    double xp = stp*i + X0 + stp/2.;
    for (int j = 0; j < nn; j++) {
      double yp = stp*j + Y0 + stp/2.;
      if(do_linear) {
        field[i][j] = interpolate(s,xp,yp,zp);
      }
      else {
        Point point = locate (xp,yp,zp);
        field[i][j] = point.level >= 0 ? s[] : nodata;
      }
    }
  }
  if (pid() == 0) { // master
#if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], sq(nn), MPI_DOUBLE, MPI_MIN, 0,
                MPI_COMM_WORLD);
#endif
    for (int i = 0; i < nn; i++) {
      for (int j = 0; j < nn; j++) {
        fwrite ( &field[i][j], sizeof(double), 1, fpver );
      }
    }
    fflush (fpver);
  }
#if _MPI
  else // slave
    MPI_Reduce (field[0], NULL, nn*nn, MPI_DOUBLE, MPI_MIN, 0,
                MPI_COMM_WORLD);
#endif
  matrix_free (field);
  fclose (fpver); // we close at the end

}

void * matrix_new_3d (int nx, int ny, int nz, size_t size) {
  void ** m = qmalloc(nx, void *);  //Define a pointer that points to every x coordinate.
  char * a  = qmalloc(nx*ny*nz*size, char);
  for (int i=0; i<nx; i++) {
    m[i] = a+i*ny*nz*size;
  }
  return m;
}

void output_2d_span_avg(char * fname, scalar s, int maxlevel, bool do_linear, bool print_bin) {

  int nn = (1<<maxlevel);
  double stp = 0.999999*(L0+X0-X0)/(double)nn; // to avoid interpolated point coincides with fine grid boundary

  //fprintf(stderr, "I am here 1\n"), fflush (stderr);
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

  //fprintf(stderr, "I am here 2\n"), fflush (stderr);
  FILE * fpver = fopen (fname,"w");
  if (pid() == 0) { // master

    // reduce it to the first pid()
#if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], cube(nn), MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
#endif

    //fprintf(stderr, "I am here 2b\n"), fflush (stderr);
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
    //fprintf(stderr, "I am here 3\n"), fflush (stderr);

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
    //fprintf(stderr, "I am here 4\n"), fflush (stderr);

  }
#if _MPI
  else // slave
    MPI_Reduce (field[0], NULL, cube(nn), MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
#endif
  //fprintf(stderr, "I am here 5\n"), fflush (stderr);

  matrix_free (field);
  fclose (fpver); // we close at the end

}

double my_interpolation(Point point, scalar s, double xp, double yp, double zp) {

  //boundary({s}); // we need to avoid this call since we use my_interpolation 
  //                  only if point.level > 0 (i.e., only on some processors)
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
    MPI_Reduce (field[0], NULL, cube(nn), MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
#endif
  matrix_free (field);
  fclose (fpver); // we close at the end

}

/**
   Compute the bulk dissipation in air and water. */

void profile_output (vector u, int istep, int MAXLEVEL) {

  char file[99];
  sprintf (file, "./profiles/prof_%09d.out", istep);

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
    ke[]   = sq(u.x[])+sq(u.y[])+sq(u.z[]);
    vke[]  = u.y[]*ke[];
  }
  foreach() {
    dkdy[] = (ke[0,1]-ke[0,-1])/(2.*Delta);
  }

  vertex scalar phi[];
  foreach_vertex() {
    phi[] = y;
  }

  profiles ({u.x, u.y, u.z, p, omega, uv, duy, diss, ke, vke, dkdy}, phi, rf = 1.0, fname = file, n = 1 << MAXLEVEL);
  delete({omega,uv,duy,diss,ke,vke,dkdy});

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
