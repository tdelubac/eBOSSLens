# import argparse
import os
import multiprocessing as mul
import numpy as np
from eBOSSLens import eBOSSLens
from utils import make_sure_path_exists


def para_return(func, params, num_thread=4):
    res_list = []
    p = mul.Pool(processes=num_thread)
    for each in params:
        def _log_result(res):
            res_list.append(res)
        p.apply_async(func, each, callback=_log_result)
    p.close()
    p.join()
    return res_list


def lensFinder(plate, mjd, fiberid, datav='v5_7_0', datadir='../SCRATCH'):
    savedir = os.path.join('../FullSearch', str(plate) + "-" + str(mjd))
    make_sure_path_exists(savedir)
    try:
        eBOSSLens(plate, mjd, fiberid, datav, False, False, False, savedir,
                  datadir)
    except Exception as reason:
        print(reason)


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
    platemjd = np.loadtxt('plate-mjd.txt', dtype=int)
    for each in platemjd:
        fiberid = np.arange(1, 1001, 1, dtype=int)
        args = []
        for fid in fiberid:
            args.append((each[0], each[1], fid,))
        res = para_return(lensFinder, args, 12)
    # lensFinder(6739, 56393, 272)
