#!/usr/bin/env bash

#SBATCH --partition=fwkt_v100
# necessary to set the account also to the queue name because otherwise access is not allowed at the moment
#SBATCH --account=fwkt_v100
#SBATCH --time=24:00:00
# Sets batch job's name
#SBATCH --job-name=wtest5
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mincpus=4
#SBATCH --cpus-per-task=16
#SBATCH --mem=320000
#SBATCH --chdir=/home/marqua27/jupyter_notebooks/wtest5_out
#SBATCH -o stdout
#SBATCH -e stderr

echo 'Running program...'
source /home/marqua27/jupyter.profile

cd /home/marqua27/jupyter_notebooks/

echo 'starting...'
python flying_focus.py 4