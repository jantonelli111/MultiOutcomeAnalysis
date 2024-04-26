#!/bin/bash
#
#SBATCH -p fasse # partition (queue)
#SBATCH -c 48 # number of cores
#SBATCH --mem 184GB 
#SBATCH -t 7-00:00 # time (D-HH:MM)

singularity exec --cleanenv --env R_LIBS_USER=../RELEASE_3_14 ../ifxrstudio:RELEASE_3_14.sif Rscript MainAnalysis.R

#cp $HOME/R/ifxrstudio/RELEASE_3_14 .
#cp /n/singularity_images/informatics/ifxrstudio/ifxrstudio:RELEASE_3_14.sif .
