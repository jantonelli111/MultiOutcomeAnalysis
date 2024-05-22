#!/bin/bash
#
#SBATCH -p serial_requeue # partition (queue)
#SBATCH -c 8 # number of cores
#SBATCH --mem 48GB 
#SBATCH -t 0-12:00 # time (D-HH:MM)

singularity exec --cleanenv --env R_LIBS_USER=../RELEASE_3_14 ../ifxrstudio:RELEASE_3_14.sif Rscript MainAnalysis.R --analysis test1

#cp $HOME/R/ifxrstudio/RELEASE_3_14 .
#cp /n/singularity_images/informatics/ifxrstudio/ifxrstudio:RELEASE_3_14.sif .
