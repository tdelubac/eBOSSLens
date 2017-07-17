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


def lensFinder(plate, mjd, fiberid, datav='v5_7_0', datadir='/SCRATCH'):
    savedir = os.path.join('../FullSearch', str(plate) + "-" + str(mjd))
    make_sure_path_exists(savedir)
    try:
        eBOSSLens(plate, mjd, fiberid, datav, False, False, False, savedir,
                  datadir)
    except Exception as reason:
        text = str(plate) + " " + str(mjd) + " " + str(fiberid) + " " + \
            str(reason)
        print(text)


if __name__ == "__main__":
    platemjd = np.loadtxt('plate-mjd.txt', dtype=int)
    for each in platemjd:
        fiberid = np.arange(1, 1001, 1, dtype=int)
        args = []
        for fid in fiberid:
            args.append((each[0], each[1], fid,))
        res = para_return(lensFinder, args, 12)
    # Uncomment below and comment above to debug
    # lensFinder(4391, 55866, 999, datadir='../SCRATCH')
