#!/bin/bash
#PBS -N tache       # Nom de la tâche
#PBS -A hii-894-aa  # Identifiant Rap; ID
#PBS -l nodes=16:ppn=8 # Number of nodes and number of processors per node
#PBS -l walltime=30:00:00    # Duration in seconds


#source /software/soft.computecanada.ca.sh
#module load compilers/intel/14.0
#module load mpi/openmpi/1.6.5

cd $PBS_O_WORKDIR

# Commande à exécuter
#mpiexec -n 128 /home/delcampo/bin/raxmlHPC-MPI-SSE3 -m GTRCAT -c 25 -e 0.1 -p 31415 -d -f d -N 1000 -n alex -s stram02.trim.phy