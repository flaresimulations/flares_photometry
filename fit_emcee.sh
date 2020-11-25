#!/bin/bash
#SBATCH --ntasks 6
#SBATCH -A dp004
#SBATCH -p cosma7
#SBATCH --job-name=get_fit
#SBATCH -t 0-2:30
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=6
#SBATCH -o logs/std_output.%J
#SBATCH -e logs/std_error.%J

module purge
module load gnu_comp/7.3.0 openmpi/3.0.1 hdf5/1.10.3 python/3.6.5

source ../flares_pipeline/venv_fl/bin/activate

mpiexec -n 6 python3 LF_fit.py Schechter


echo "Job done, info follows..."
sacct -j $SLURM_JOBID --format=JobID,JobName,Partition,MaxRSS,Elapsed,ExitCode
exit
