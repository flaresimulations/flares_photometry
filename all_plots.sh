#!/bin/bash
#SBATCH --ntasks 1
#SBATCH -A dp004
#SBATCH -p cosma7
#SBATCH --job-name=get_plots
#SBATCH -t 0-0:30
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH -o logs/std_output.%J
#SBATCH -e logs/std_error.%J

module purge
module load gnu_comp/7.3.0 openmpi/3.0.1 hdf5/1.10.3 python/3.6.5


#export PY_INSTALL=/cosma/home/dp004/dc-love2/.conda/envs/eagle/bin/python

source ../flares_pipeline/venv_fl/bin/activate

python3 gsmf_compare.py
python3 Z_evo.py
python3 LF.py 0
python3 LF.py 1
python3 LF.py 2
python3 fit_evo.py
python3 beta_relation.py 0
python3 att_lum.py 0
python3 att_lum.py 1
python3 att_lum.py 2
python3 UVLF_env.py 0
python3 line_lumEW_stamps.py 0
python3 line_lumEW_stamps.py 1
python3 line_lumEW_stamps.py 2
python3 line_lumOIII_EW.py
python3 line_lumCIII_EW.py
python3 LFOIII_deBarros.py
python3 sfr_obscured.py 0
python3 sfr_obscured.py 1



echo "Job done, info follows..."
sacct -j $SLURM_JOBID --format=JobID,JobName,Partition,MaxRSS,Elapsed,ExitCode
exit
