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
#SBATCH --output="%A-lb.out"

./experiment1.sh 