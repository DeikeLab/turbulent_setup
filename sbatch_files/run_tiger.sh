#!/bin/bash
#SBATCH --job-name=720_0p3          # create a short name for your job
#SBATCH --nodes=4                   # node count
#SBATCH --ntasks-per-node=40        # number of tasks per node
#SBATCH --cpus-per-task=1           # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G            # memory per cpu-core (4G is default)
#SBATCH --time=23:59:00             # total run time limit (HH:MM:SS)
##SBATCH --time=00:30:00             # total run time limit (HH:MM:SS)
####SBATCH --mail-type=begin        # send email when job begins
####SBATCH --mail-type=end          # send email when job ends
####SBATCH --mail-user=ns8802@princeton.edu
#
module load openmpi/gcc/3.1.5/64
#
ulimit -s unlimited
#
# Input parameters
#
Re_ast=720.0;
BO=200.0;
Re_wave=2.5597571943e+04;
UstarRATIO=0.9;
ak=0.3;
r_L0lam=4.0;
rho_r=0.001225;           
mu_r=2.2471881948940954e-06; 
MAXLEVEL=10;
MINLEVEL=6;
RELEASETIME=3.4324016086e+03;
uemaxRATIO=0.3;
end_sim=1000000.0;
do_eta_loc=1;
do_profile=1;
do_fields=1;
do_tagging=0;
prt_res=9; 
from_pr=1;
#
srun precursor -m 23:59:00 $Re_ast $BO $Re_wave $UstarRATIO $ak $r_L0lam \
	                   $rho_r $mu_r $MAXLEVEL $MINLEVEL $RELEASETIME \
	                   $uemaxRATIO $end_sim \
			   $do_eta_loc $do_profile $do_fields $do_tagging \
	                   $prt_res $from_pr > out.log 2>&1
#
