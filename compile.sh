#
mkdir -p windwave_turb
#
# compile
#
CC99='mpicc -std=c99' qcc -grid=octree -events -autolink -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -D_MPI=1 -DTREE -o windwave_turb/windwave_turb curved_fixREtau.c -L$BASILISK/gl -lfb_tiny -lm
