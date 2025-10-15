/**
# Spectrum read-in and wave field initialization functions */

double *F_kxky_, *phase_;
double *kx_, *ky_;
double dkx_, dky_;

/**
   Allocate arrays based on N_mode_ */

void allocate_arrays(int N_mode_) {
  size_t N2D = N_mode_ * (N_mode_ + 1);
  size_t N1D = N_mode_;

  F_kxky_ = (double *) malloc(N2D * sizeof(double));
  phase_  = (double *) malloc(N2D * sizeof(double));
  kx_     = (double *) malloc(N1D * sizeof(double));
  ky_     = (double *) malloc((N1D + 1) * sizeof(double));

  if (!F_kxky_ || !phase_ || !kx_ || !ky_) {
    fprintf(stderr, "Memory allocation failed!\n");
    exit(EXIT_FAILURE);
  }
}

/**
   Deallocate them */

void free_arrays() {
  if (F_kxky_) free(F_kxky_);
  if (phase_)  free(phase_);
  if (kx_)     free(kx_);
  if (ky_)     free(ky_);

  F_kxky_ = NULL;
  phase_  = NULL;
  kx_     = NULL;
  ky_     = NULL;
}

/**
   Function to read it */

void load_array(char *filename, double *target, int expected_count, char *label) {
  float *buffer = (float *) malloc(expected_count * sizeof(float));
  if (!buffer) {
    fprintf(stderr, "Fatal error: malloc failed for buffer loading '%s'\n", filename);
#if _MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
#else
    exit(1);
#endif
  }
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Fatal error: cannot open file '%s'\n", filename);
#if _MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
#else
    exit(1);
#endif
  }

  fseek(fp, 0, SEEK_END);
  long file_size = ftell(fp);
  rewind(fp);

  long expected_size = expected_count * sizeof(float);
  if (file_size != expected_size) {
    fprintf(stderr,
        "Fatal error: file '%s' has size %ld bytes, expected %ld bytes\n",
        filename, file_size, expected_size);
    fclose(fp);
#if _MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
#else
    exit(1);
#endif
  }

  size_t read_count = fread(buffer, sizeof(float), expected_count, fp);
  fclose(fp);
  if (read_count != (size_t)expected_count) {
    fprintf(stderr,
        "Fatal error: expected %d elements in '%s', but read %zu\n",
        expected_count, filename, read_count);
    free(buffer);
#if _MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
#else
    exit(1);
#endif
  }

  for (int i = 0; i < expected_count; i++) {
    target[i] = (double) buffer[i];
  }
  
  free(buffer);

  fprintf(stderr, "%s loaded from '%s'\n", label, filename);
  fflush(stderr);
}

/**
   Random number generator. */

double randInRange(int min, int max) {
  return min + (rand() / (double) (RAND_MAX) * (max - min + 1));
}

/** 
A MPI compatible function that reads in kx_, ky_, and F_kxky_. For now kx_, ky_ are 1D arrays while F_kxky_ is a 2D array. The root process reads in and then broadcast to all other processes. It looks for files named F_kxky, kx, ky. If they are absent, it returns an error. */

void power_input(int N_mode_) {

  int length1D = N_mode_;
  int length2D = N_mode_*(N_mode_+1);

#if _MPI

  int rank, size;
  int root = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Read in on the zeroth thread
  if (rank == root) {
    load_array("F_kxky", F_kxky_, length2D  , "F_kxky_");
    load_array("kx",     kx_,     length1D  , "kx_");
    load_array("ky",     ky_,     length1D+1, "ky_");
    
    int index = 0;
    srand(MY_RAND); // We can seed it differently for different runs
    fprintf(stderr, "MY_RAND Num. is: %d\n", MY_RAND), fflush (stderr);
    for (int i=0; i<N_mode_; i++) {
      for (int j=0; j<N_mode_+1; j++) {
	index = j*N_mode_ + i;
	phase_[index] = randInRange (0, 2.*pi);
      }
    }
  }

  // Broadcast to other threads
  MPI_Bcast(kx_    , length1D  , MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Bcast(ky_    , length1D+1, MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Bcast(F_kxky_, length2D  , MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Bcast(phase_ , length2D  , MPI_DOUBLE, root, MPI_COMM_WORLD);
  fprintf(stderr, "The MPIBcast is done!\n"), fflush (stderr);

#else
  
  load_array("F_kxky", F_kxky_, length2D  , "F_kxky_");
  load_array("kx",     kx_,     length1D  , "kx_");
  load_array("ky",     ky_,     length1D+1, "ky_");

  int index = 0;
  srand(MY_RAND); // We can seed it differently for different runs
  for (int i=0; i<N_mode_;i++) {
    for (int j=0;j<N_mode_+1;j++) {
      index = j*N_mode_ + i;
      phase_[index] = randInRange (0, 2.*pi);
    }
  }

#endif

}

/**
   Function to initialize the wave field (based on linear wave equations). */

double wave (double x, double y, int N_mode_) {
  double eta = 0;
  for (int i=0; i<N_mode_; i++) {
    for (int j=0; j<N_mode_+1; j++) {
      int index = j*N_mode_ + i;
      double ampl = sqrt(2.*F_kxky_[index]*dkx_*dky_);
      double a = (kx_[i]*x + ky_[j]*y + phase_[index]);
      eta += ampl*cos(a);
    }
  }
  return eta;
}

/**
   Functions that compute the orbital velocity (based on linear wave equations). */

double uin_x (double x, double y, double z, int N_mode_) {
  double uin_x = 0;
  for (int i=0; i<N_mode_; i++) {
    for (int j=0; j<N_mode_+1; j++) {
      int index = j*N_mode_ + i;
      double ampl = sqrt(2.*F_kxky_[index]*dkx_*dky_);
      double z_actual = (z < ampl ? (z) : ampl);
      // fprintf(stderr, "z = %g, ampl = %g, z_actual = %g\n", z, ampl, z_actual);
      double kmod = sqrt(sq(kx_[i]) + sq(ky_[j]));
      double theta = atan(ky_[j]/kx_[i]);
      double a = (kx_[i]*x + ky_[j]*y + phase_[index]);
      uin_x += sqrt(g_*kmod)*ampl*exp(kmod*z_actual)*cos(a)*cos(theta);
    }
  }
  return uin_x;
}

double uin_y (double x, double y, double z, int N_mode_) {
  double uin_y = 0;
  for (int i=0; i<N_mode_; i++) {
    for (int j=0; j<N_mode_+1; j++) {
      int index = j*N_mode_ + i;
      double ampl = sqrt(2.*F_kxky_[index]*dkx_*dky_);
      double z_actual = (z < ampl ? (z) : ampl);
      double kmod = sqrt(sq(kx_[i]) + sq(ky_[j]));
      double theta = atan(ky_[j]/kx_[i]);
      double a = (kx_[i]*x + ky_[j]*y + phase_[index]);
      uin_y += sqrt(g_*kmod)*ampl*exp(kmod*z_actual)*cos(a)*sin(theta);
    }
  }
  return uin_y;
}

double uin_z (double x, double y, double z, int N_mode_) {
  double uin_z = 0;
  for (int i=0; i<N_mode_; i++) {
    for (int j=0; j<N_mode_+1; j++) {
      int index = j*N_mode_ + i;
      double ampl = sqrt(2.*F_kxky_[index]*dkx_*dky_);
      double z_actual = (z < ampl ? (z) : ampl);
      double kmod = sqrt(sq(kx_[i]) + sq(ky_[j]));
      double a = (kx_[i]*x + ky_[j]*y + phase_[index]);
      uin_z += sqrt(g_*kmod)*ampl*exp(kmod*z_actual)*sin(a);
    }
  }
  return uin_z;
}
