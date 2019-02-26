#!/usr/bin/env python

__author__ = "Xiang Gao, Huaiying Lin, Qunfeng Dong"
__copyright__ = "Copyright 2016, Bayesian-based LCA taxonomic classification"
__credits__ = ["Xiang Gao", "Huaiying Lin","Qunfeng Dong"]
__license__ = "GPL"
__version__ = "2.1"
__maintainer__ = "Huaiying Lin"
__email__ = "hlin2@luc.edu"
__status__ = "Python3"

import argparse
import os
import sys
import pandas as pd
import numpy as np
import time

#### constant #####
RANKNAMES = ["superkingdom","phylum","class","order","family","genus","species"]
headers = [(r,r+"-score") for r in RANKNAMES]
headers = [y for z in headers for y in z]

##### function #####
def float_100range(x):
    x = float(x)
    if x < 0.0 or x > 100.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 100.0]"%(x,))
    return x

def blca_parser(blcaout, level, cutoff, rel):
    na = os.path.splitext(os.path.split(os.path.basename(blcaout))[1])[0]
    try:
        blca = pd.read_csv(blcaout, sep="\t", header=None)
        blca.columns = ["ReadID", "Taxonomy"]
        tx = blca["Taxonomy"].str.split(";|:", expand=True)
        tx.replace("", np.nan, inplace=True)
        tx = tx.dropna(how='all', axis=1)
    except:
        print("Your input file is not in BLCA output format!")
        sys.exit(1)
    tx = tx.drop([x for x in range(0,tx.shape[1],3)],axis=1)
    tx.columns = headers
    total = tx.shape[0]
    tx[[r+"-score" for r in RANKNAMES]] = tx[[r+"-score" for r in RANKNAMES]].astype(float)
    sel_tx = tx.loc[:,[s.startswith(level) for s in tx.columns]]
    sel_tx = sel_tx.loc[sel_tx[level+"-score"] > cutoff,:]
    outtx = sel_tx.groupby(level).agg("count")
    outtx.columns = [na]
    uclRow = pd.Series({na: total - int(outtx.sum(axis=0))}, name="Unclassified")
    outtx = outtx.append(uclRow)
    if rel:
        outtx[na] = outtx[na]/total * 100
    return outtx

def loop_wrapper(list_of_blcaouts,level, cutoff, rel):
    fF = list_of_blcaouts.pop(0)
    merDF = blca_parser(fF, level, cutoff, rel)
    for f in list_of_blcaouts:
        tmpDF = blca_parser(f, level, cutoff, rel)
        merDF = merDF.merge(tmpDF, on=level, how="outer").fillna(0)
    return merDF

def main():
    Nf = len(args.files)
    startTime = time.time()
    outFile = loop_wrapper(args.files, args.tax, args.cutoff, not args.rawCounts)
    outFile.to_csv(args.output,sep='\t')
    endTime = time.time()
    print("Merged %i files" % Nf)
    print("Time elapsed: %i seconds" % (round((endTime - startTime))))


##### parser ######

parser = argparse.ArgumentParser(prog = "generate_abundance_table.py", description=
     '''    Utility script to combine all BLCA outputs and generate abundance table.
    
    Sample ID will be automatically extracted from the file name as the ID used in the combined table.
    
    Multiple files or wild card file pattern are accepted.''',
     epilog="No warrenty comes with this script. Author: hlin2@luc.edu. \nAny suggestions or bugs report are welcomed.",
     add_help=False, formatter_class=argparse.RawTextHelpFormatter)
##### Required arguments #####
required = parser.add_argument_group('positional arguments')
required.add_argument("files", metavar = "input.blca.out", nargs = "+", help="Input BLCA output files.")
##### Reads filter #####
readsoptions = parser.add_argument_group('reads filter arguments')
readsoptions.add_argument("-c", "--cutoff", default=80.0, help="minimum confidence score to be considered as a valid taxonomy assignment. Default: 80", type=float_100range)
##### Other arguments #####
optional = parser.add_argument_group('other arguments')
optional.add_argument("-t", "--tax", default='genus',
                      help="taxonomy level to construct the table. Default: genus", type=str,
                      choices=["phylum","class","order","family","genus","species"])
optional.add_argument("--rawCounts", help="use raw counts instead of relative abundance",
                      action='store_true')
optional.add_argument("-o", "--output",
                      help="output file name. Default: merged.table.<level>.txt", type=str)
optional.add_argument("-h", "--help", help="show this help message and exit", action="help")
##### parse arguments #####
args = parser.parse_args()
##### Output file name ####
if not args.output:
    args.output="merged.table."+args.tax+".txt"

# print(args.files)
if __name__ == "__main__":
    main()