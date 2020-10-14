#!/bin/bash -l
#SBATCH -A SNIC2020-5-300
#SBATCH -t 01-10:00:00
#SBATCH --job-name=Seqvec_PS1
#SBATCH --error=Seqvec_PS1.error
#SBATCH --output=Seqvec_PS1.out

for filename in *.fa;
    do
        if [ -f "${filename%.*}".npz ]; then
                echo "$filename"
        else
                echo "${filename%.*}"
                seqvec -i $filename -o "${filename%.*}".npz --protein True
        fi
    done
~                                                                                                                                                                                                                                                                                 
~                                                                                                                                                                                                                                                                                 
~                                                                                                                                                                                                                                                                                 
~                                                                                                                                                                                                                                                                                 
~               
