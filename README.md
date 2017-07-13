# eBOSSLens
## Currently under refactoring. Galaxy lensing should be working. 
## Use at own risk

## TO-DO list:
* Complete objectifing the code for QSO, LyA, jackpot lensing
* Further improve in parameter configuration (i.e. using json config file)
* Write a better document

## Note for further improver
* Search for `#TODO` in the project for the parts need to be modified
* Parameters/returns of functions might also be modified if necessary
* Use `os.path.join` instead of adding strings to form directory
* `topdir` has been removed since it could be merged into `savedir`
* Main calculation function is `eBOSSLens.eBOSSLens`
* `main.py` is only an example of how to invoke the function
* Document could be written in source code file
