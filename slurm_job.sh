#!/bin/bash
#SBATCH --job-name=docking_sim        
#SBATCH --nodes=1               
#SBATCH --ntasks=1              
#SBATCH --cpus-per-task=27     
#SBATCH --mem=1G                 
#SBATCH --time=00:10:00
#SBATCH --account=m24028
#SBATCH --partition=mesonet
#SBATCH --output=docking_output.out
#SBATCH --error=docking_error.err

lscpu | grep -i "CPU(s):" | head -n 1
time conda run -n chps --live-stream python generate_cubes.py
