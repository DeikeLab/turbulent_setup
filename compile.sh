#
#module load intel-mpi/gcc/2021.3.1
module purge
module load intel-oneapi/2024.2
module load intel-mpi/oneapi/2021.13
module load intel-tbb/2021.13 intel-rt/2024.2

rm -rf windwave_turb && mkdir -p windwave_turb
mkdir -p windwave_turb
#
# compile
#
CC99='mpicc -std=c99' qcc -grid=octree -events -autolink -O2 -Wall -pipe -D_FORTIFY_SOURCE=2 -D_MPI=1 -disable-dimensions -o windwave_turb/windwave_turb curved_fixREtau.c -L$BASILISK/gl -lfb_tiny -lm
#qcc -grid=octree -events -autolink -O2 -Wall -pipe -D_FORTIFY_SOURCE=2 -D_MPI=0 -disable-dimensions -o windwave_turb/windwave_turb curved_fixREtau.c -L$BASILISK/gl -lfb_tiny -lm
