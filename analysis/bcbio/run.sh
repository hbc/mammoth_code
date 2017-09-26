#!/bin/sh
#BSUB -q priority
#BSUB -J mammoth
#BSUB -oo main.out
#BSUB -n 1
#BSUB -R "rusage[mem=8024]"
#BSUB -W 150:00

date
unset JAVA_HOME

bcbio_nextgen.py ../config/mammoth_vc.yaml -n 48 -t ipython -s lsf -q mcore '-rW=150:00' -r mincores=2 -rminconcores=2 --retries 3 --timeout 180 --tag mammoth

date
