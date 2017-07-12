# import argparse
import numpy as np
from eBOSSLens import eBOSSLens


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
    plate = 6739
    mjd = 56393
    fiberid = np.arange(1, 1000, 1, dtype=int)
    for each in fiberid:
        try:
            eBOSSLens(plate, mjd, each, False, False, False)
        except Exception as reason:
            print(str(each) + " " + str(reason))
