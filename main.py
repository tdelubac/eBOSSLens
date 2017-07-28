import argparse
import os
import multiprocessing as mul
import numpy as np
from eBOSSLens import eBOSSLens
from plateStats import plateStats
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


def lensFinder(plate, mjd, fiberid, datav, datadir, savedir, lya, qso, jpt,
               bwidth, bsig, maxchi2, doplot):
    sd = os.path.join(savedir, str(plate) + "-" + str(mjd))
    make_sure_path_exists(sd)
    try:
        eBOSSLens(plate, mjd, fiberid, datav, lya, qso, jpt, sd, datadir,
                  max_chi2=maxchi2, bwidth=bwidth, bsig=bsig, doPlot=doplot)
    except Exception as reason:
        text = str(plate) + " " + str(mjd) + " " + str(fiberid) + " " + \
            str(reason)
        print(text)


def setArg(argument, defValue):
    if argument is not None:
        return argument
    else:
        return defValue


if __name__ == "__main__":
    # Command line argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("pmfile", help="The file for plate and mjd input",
                        type=str)
    parser.add_argument("-c", help="Max accepted chisquare", type=float)
    parser.add_argument("-w", help="Width for peak detection", type=float)
    parser.add_argument("-s", help="Sig for peak detection", type=float)
    parser.add_argument("--dataversion", help="Data version string", type=str)
    parser.add_argument("--datadir", help="Data directory string", type=str)
    parser.add_argument("--savedir", help="Save directory string", type=str)
    parser.add_argument("--lya", help="LyA search", action="store_true")
    parser.add_argument("--qso", help="QSO search", action="store_true")
    parser.add_argument("--jpt", help="Jackpot search", action="store_true")
    args = parser.parse_args()
    # Storing command line argument
    pmfile = args.pmfile
    maxchi2 = setArg(args.c, 2.5)
    bwidth = setArg(args.w, 60.0)
    bsig = setArg(args.s, 1.2)
    dataversion = setArg(args.dataversion, 'v5_7_0')
    datadir = setArg(args.datadir, '/SCRATCH')
    savedir = setArg(args.savedir, '../FullSearch')
    lya = setArg(args.lya, False)
    qso = setArg(args.qso, False)
    jpt = setArg(args.jpt, False)
    # Read and begin
    platemjd = np.loadtxt(pmfile, dtype=int)
    for each in platemjd:
        fiberid = np.arange(1, 1001, 1, dtype=int)
        args = []
        for fid in fiberid:
            args.append((each[0], each[1], fid, dataversion, datadir, savedir,
                         lya, qso, jpt, bwidth, bsig, maxchi2, True))
        res = para_return(lensFinder, args, 12)
        plateStats(each[0], each[1], savedir)
    # Uncomment below and comment above to debug
    # lensFinder(4198, 55480, 14, 'v5_7_0', '../SCRATCH', '../PlotCheck',
    #            False, False, False, 60.0, 1.2, 2.5, True)
