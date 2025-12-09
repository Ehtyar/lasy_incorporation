#!/usr/bin/env bash

#SBATCH --partition=gpu_a100_low
# necessary to set the account also to the queue name because otherwise access is not allowed at the moment
#SBATCH --account=low
#SBATCH --time=24:00:00
# Sets batch job's name
#SBATCH --job-name=fresnel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mincpus=4
#SBATCH --cpus-per-task=16
#SBATCH --mem=4024000
#SBATCH --chdir=/home/marqua27/jupyter_notebooks/flfoc_fresnel_out
#SBATCH -o stdout
#SBATCH -e stderr

echo 'Running program...'
source /home/marqua27/jupyter.profile

cd /home/marqua27/jupyter_notebooks/

echo 'starting...'
python runflfoc_fresnel.py $1 $2
#python runflfoc.py 98 15
#python runflfoc.py 102
