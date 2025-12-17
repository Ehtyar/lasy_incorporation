#!/usr/bin/env bash

#SBATCH --partition=fwkt_v100
# necessary to set the account also to the queue name because otherwise access is not allowed at the moment
#SBATCH --account=fwkt_v100
#SBATCH --time=24:00:00
# Sets batch job's name
#SBATCH --job-name=direct_flfoc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mincpus=4
#SBATCH --cpus-per-task=16
#SBATCH --mem=378000
#SBATCH --chdir=/home/marqua27/jupyter_notebooks/flfoc_direct_out/ono
#SBATCH -o stdout
#SBATCH -e stderr

echo 'Running program...'
source /home/marqua27/jupyter.profile

cd /home/marqua27/jupyter_notebooks/

echo 'starting...'
python runflfoc_direct.py no 20
#python runflfoc.py 98
#python runflfoc.py 102
