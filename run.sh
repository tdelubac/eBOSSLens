#! /bin/bash
#SBATCH -p r3 -o my.stdout -e my.stderr --mail-user=xinlun.cheng@epfl.ch --mail-type=ALL
python main.py plate-mjd.txt -c 3.5 --savedir ../FullSearch/3.5
python main.py plate-mjd.txt -c 3.0 --savedir ../FullSearch/3.0
