#
mkdir -p windwave_turb_linear_old_eta
#
# compile
#
CC99='mpicc -std=c99' qcc -grid=octree -events -autolink -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -D_MPI=1 -DTREE -o windwave_turb_linear_old_eta/windwave_turb_linear_old_eta curved_fixREtau_linear_old_eta.c -L$BASILISK/gl -lfb_tiny -lm
