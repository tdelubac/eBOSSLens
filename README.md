# eBOSSLens
## Currently under refactoring. Galaxy lensing should be working. 
## Refactoring is done in structure. Optimization is underway.
## Use at own risk



## TO-DO list:
* Full search on all cases
* Blind search parallelization optimization
* Further improve in parameter configuration (i.e. using json config file)
* Write a better document (invoke function descriptions + minimal example)


## Notes

* The user must call the main.py file to launch the code. The main.py must be provided with a txt file containing either plates/mjd in two columns for a full plate search or plates/mjd/fiberid for specific searches. Please note that at least two entries are necessary at the moment. 
* The main search is located in eBOSSLens.eBOSSLens(...)
* The type of background and foreground object is a galaxy by default. To change it, use the options:
 --qso for a foreground QSO
 --jpt for a Jackpot search (two background ELGs)
 --lya for a background LAE
* The user must also provide a saving directory with --savedir

## Minimal examples

* python main.py ../Test/plates_mjd_list.txt -savedir ../Test/ 
* python main.py ../Test/plates_mjd_list.txt -savedir ../Test/ --qso --lya
