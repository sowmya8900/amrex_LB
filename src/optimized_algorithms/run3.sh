#!/bin/bash
##SBATCH --qos=shared
#SBATCH --qos=regular
##SBATCH --qos=debug
#SBATCH --constraint=cpu
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --mem=1GB
#SBATCH --account=nintern
##SBATCH --mail-type=begin,end,fail
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=ananda@ucsd.edu
#SBATCH --output="%A-lb.out"

./experiment3.sh 