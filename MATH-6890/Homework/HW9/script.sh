#!/bin/sh
#SBATCH --ntasks=2
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --output=op.txt

mpirun ./heat2d.bin -nx=100000
