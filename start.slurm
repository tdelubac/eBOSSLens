#!/bin/bash 
#SBATCH -p r3n -o my.stdout -e my.stderr --mail-user=romainalexis.meyer@epfl.ch --mail-type=ALL 

python eBOSSLens.py
