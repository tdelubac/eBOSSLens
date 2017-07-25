#! /bin/bash
#SBATCH -p r3 -o my.stdout -e my.stderr --mail-user=xinlun.cheng@epfl.ch --mail-type=ALL
python main.py plate-mjd.txt --savedir ../FullSearch
