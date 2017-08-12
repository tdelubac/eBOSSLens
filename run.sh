#! /bin/bash
#SBATCH -p r4 -o my.stdout -e my.stderr --mail-user=romainalexis.meyer@epfl.ch --mail-type=ALL
python main.py ../RefactorTest/test_list.txt --qso --lya
