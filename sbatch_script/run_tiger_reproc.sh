#!/bin/bash
#SBATCH --job-name=reproc            # create a short name for your job
#SBATCH --nodes=4                   # node count
#SBATCH --ntasks-per-node=96        # number of tasks per node
#SBATCH --cpus-per-task=1           # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G            # memory per cpu-core (4G is default)
#SBATCH --time=23:59:00             # total run time limit (HH:MM:SS)
#
#module load openmpi/gcc/4.1.0
#module load intel-mpi/gcc/2021.3.1
module load intel-mpi/gcc/2021.15
#
ulimit -s unlimited
#
Re_ast=720.0;
BO=200.0;
Re_wave=2.5597571943e+04;
UstarRATIO=0.9;
r_L0lam=4.0;
rho_r=0.001225;
mu_r=2.2471881948940954e-06; 
MAXLEVEL=10;
MINLEVEL=6;
prt_res=9;
RELEASETIME=3.2708e+03;
#
DUMP_DIR="restart_bk"
#
for dump_file in "$DUMP_DIR"/dump_*.bin; do
    
    # Create symbolic link
    ln -sf "$dump_file" restart.bin

    file_number=$(basename "$dump_file" | sed -E 's/dump_([0-9]+)\.bin/\1/')

    # Print the current file number
    echo "Running with file number: $file_number"

    # Define log file names
    OUT_LOG="out_${file_number}.log"
    ERR_LOG="err_${file_number}.log"
    
    # Run your code
    srun reprocess -m 23:30:00 "$file_number" "$Re_ast" "$BO" "$Re_wave" "$UstarRATIO" "$r_L0lam" "$rho_r" "$mu_r" "$MAXLEVEL" "$MINLEVEL" "$prt_res" "$RELEASETIME" > "$OUT_LOG" 2> "$ERR_LOG"

    # Optionally check for errors
    if [[ $? -ne 0 ]]; then
        echo "Job with file number $file_number failed."
        exit 1
    fi

done
#
