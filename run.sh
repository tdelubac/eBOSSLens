#! /bin/bash
#SBATCH -p r3 -o my.stdout -e my.stderr --mail-user=xinlun.cheng@epfl.ch --mail-type=ALL
python main.py CFHTLS_Wide.txt --savedir ../CFHTLS_Wide
python genCatalog.py ../CFHTLS_Wide ../CFHTLS_Wide/CFHTLS_candidates.fits
