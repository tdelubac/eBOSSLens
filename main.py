# import argparse
import os
import numpy as np
from eBOSSLens import eBOSSLens
from utils import make_sure_path_exists


if __name__ == "__main__":
    '''
    # Required argument
    parser = argparse.ArgumentParser()
    parser.add_argument("plate", type=int, help="Plate number")
    parser.add_argument("mjd", type=int, help="MJD")
    parser.add_argument("fiberid", type=int, help="FiberId")
    # parser.add_argument("chi2", type=float)
    # parser.add_argument("width", type=float)
    # parser.add_argument("sig", type=float)
    # parser.add_argument("sdir", type=str)
    # parser.add_argument("dataVersion", type=str)
    args = parser.parse_args()
    eBOSSLens(args.plate, args.mjd, args.fiberid, False, False, False)
    '''
    plate = 4391
    mjd = 55866
    fiberid = np.arange(1, 1001, 1, dtype=int)
    datav = 'v5_7_0'
    datadir = '../SCRATCH'
    savedir = os.path.join('../FullSearch', str(plate) + "-" + str(mjd))
    make_sure_path_exists(savedir)
    for each in fiberid:
        try:
            eBOSSLens(plate, mjd, each, datav, False, False, False, savedir,
                      datadir)
        except Exception as reason:
            print(str(each) + " " + str(reason))
