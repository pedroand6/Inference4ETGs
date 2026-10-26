#!/bin/bash
# Example Slurm job array for the spectral fits: one MultiNest fit per array task, on one core.
# Submit from the folder that holds run_spectral_fits.py, synthetic_model.py, splus-filters/ and
# spectral_fits/jobs/ (see the docstring of run_spectral_fits.py), with:   sbatch slurm_fit.sh
#
# What to edit:
#   --array     0-<N-1>%<P>: N is the number of .npz files in spectral_fits/jobs/ (201 for 33 mocks and 168
#               galaxies, so 0-200); P is how many tasks may run at once (the cores you may use).
#   --time      a fit may take several hours; raise this if tasks are killed (finished fits are kept).
#   LD_LIBRARY_PATH   the folder with libmultinest.so.
#   the activation line   your Python environment with bagpipes and pymultinest.
#
# Resubmitting is safe: a task whose <job_id>.h5 already exists finishes at once.
#SBATCH --job-name=spectral_fit
#SBATCH --array=0-200%20
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=12:00:00
#SBATCH --output=slurm_%A_%a.log

export LD_LIBRARY_PATH=$HOME/MultiNest/lib:$LD_LIBRARY_PATH
source $HOME/venv/bin/activate

# --task-index picks the job at this position in the sorted list of .npz files
python run_spectral_fits.py --task-index $SLURM_ARRAY_TASK_ID
