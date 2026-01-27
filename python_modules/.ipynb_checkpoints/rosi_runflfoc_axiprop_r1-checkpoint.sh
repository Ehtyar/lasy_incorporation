#!/usr/bin/env bash

#SBATCH --partition=cpu-genoa
# necessary to set the account also to the queue name because otherwise access is not allowed at the moment
#SBATCH --time=1-00:00:00
# Sets batch job's name
#SBATCH --job-name=r1axiprop
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mincpus=4
#SBATCH --cpus-per-task=16
#SBATCH --mem=1500000
#SBATCH --chdir=/home/marqua27/jupyter_notebooks/rosi_r1_flfoc_axiprop_out
#SBATCH -o stdout
#SBATCH -e stderr

echo 'Running program...'
source /home/marqua27/jupyter.profile

cd /home/marqua27/jupyter_notebooks/

echo 'starting...'
python runflfoc_axiprop_xyt.py $1 $2 "rosi" "r1" # percent c, num iterations, profile settings

