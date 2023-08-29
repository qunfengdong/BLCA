#!/usr/bin/env python

__author__ = "Xiang Gao, Huaiying Lin, Qunfeng Dong"
__copyright__ = "Copyright 2016, Bayesian-based LCA taxonomic classification"
__credits__ = ["Xiang Gao", "Huaiying Lin", "Qunfeng Dong"]
__license__ = "GPL"
__version__ = "1.3"
__maintainer__ = "Qunfeng Dong"
__email__ = "qunfengd@gmail.com"
__status__ = "Python3"

import sys
import math
import glob
import re
import os
import getopt
import subprocess
import argparse

'''

Get a subset of lineage information from NCBI 16S Microbial database.

'''

##### parser ######

parser = argparse.ArgumentParser(description=''' << Bayesian-based LCA taxonomic classification method >>\n\n   Please make sure the following softwares are in your PATH:
         1.muscle (http://www.drive5.com/muscle/downloads.htm), muscle should be the program's name.
         2.ncbi-blast suite (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
         3.clustalo (http://www.clustal.org/omega/), clustalo should be the program's name.
         4.Biopython should be installed locally.''',
                                 epilog="No warrenty comes with this script. Author: hlin2@luc.edu. \nAny suggestions or bugs report are welcomed.",
                                 add_help=False, formatter_class=argparse.RawTextHelpFormatter)
##### Other arguments #####
optional = parser.add_argument_group('optional arguments')
optional.add_argument("--dir", default='db',
                      help="The local directory name where you want to store the formatted database. Default: db", type=str)
optional.add_argument("-d", "--database", default='https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz',
                      help="The database link that you want to download from and format. Default: https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz", type=str)
optional.add_argument("--taxdmp", default='ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip',
                                          help="The taxonomy database dmp link from NCBI. Default: ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip", type=str)
optional.add_argument(
    "-h", "--help", help="show this help message and exit", action="help")
##### parse arguments #####
args = parser.parse_args()

nodes = args.dir + '/nodes.dmp'
names = args.dir + '/names.dmp'
merged = args.dir + '/merged.dmp' #some NCBI taxIDs were merged together by the NCBI taxonomy database
				#we need to use the one that are associated with lineage information

################ Function ##############


def setup_tax_db(url, dbstr):
    os.system("wget -P ./" + dbstr + "/ " + url)
    os.system("tar zxvf ./" + dbstr + "/taxdb.tar.gz -C " + dbstr + "/")
    print("taxdb downloaded!\n######################\n>> Please run the following before the next script:\nexport BLASTDB=" +
          os.getcwd() + "/" + dbstr + "/" + "\n######################")


def lineage(tid, allnode, namesdic, nonexist):
    levels = ["superkingdom", "phylum", "class", "order",
              "family", "genus", "subspecies", "species"]
    if tid not in nonexist:
        r = str()
        name = namesdic[tid][0]
        pa, rank = allnode[tid]
        if rank in levels:
            r += ":".join([rank, name]) + ";"
        if pa in nonexist:
            r += ""
        else:
            r += lineage(pa, allnode, namesdic, nonexist)
        return r
    else:
        return ""


def get_taxdb(url, dbstr):
    os.system("wget -P " + dbstr + "/ " + url)
    os.system("unzip " + dbstr + "/taxdmp.zip -d " + dbstr + "/")
    print(">> NCBI Taxonomy Database downloaded!")


def gi_taxid_nucl(fsa, gi_taxid_dmp):
    fa = open(fsa)
    gilist = {}
    for ln in fa:
        lne = ln.rstrip()
        if lne[0] == ">":
            gilist[lne.split()[0].replace(">", "")] = lne.split("|")[1]
    fa.close()
    glst = gilist.values()
    tagi = {}
    dmp = open(gi_taxid_dmp)
    for l in dmp:
        ln = l.rstrip()
        line = ln.split("\t")
        if line[0] in glst:
            tagi[line[0]] = line[1]
            glst.remove(line[0])
    dmp.close()
    outmap = open(fsa + ".taxid_mapping", 'wa')
    for k, v in gilist.items():
        outmap.write(k + "\t" + tagi[v] + "\n")
    outmap.close()


def check_program(prgname):
    '''Check whether a program has been installed and put in the PATH'''
    import os
    path = os.popen("which " + prgname).read().rstrip()
    if len(path) > 0 and os.path.exists(path):
        print(prgname + " is located in your PATH!")
    else:
        print("ERROR: " + prgname +
              " is NOT in your PATH, please set up " + prgname + "!")
        sys.exit(1)


def get_db_gitaxid(url, dbstr):
    dblist = re.split(r'\.|\/+', url)[-3]
    check_program("blastdbcmd")
    os.system("wget -P " + dbstr + "/ " + url)
    os.system("tar zxvf " + dbstr + "/" + dblist + ".tar.gz -C " + dbstr + "/")
    subprocess.call(['blastdbcmd', '-entry', 'all', '-outfmt', '%a|%T', '-db',
                     dbstr + '/' + dblist, '-out', dbstr + '/' + dblist + '.ACCtaxID'])
    print(">> Accession and TaxIDs from " + dblist + " are extracted!!")

#####################################################################
# End of functions
# Working scripts
####
#####################################################################


os.system("mkdir " + args.dir)
get_taxdb(args.taxdmp, args.dir)
get_db_gitaxid(args.database, args.dir)
dbname = re.split(r'\.|\/+', args.database)[-3]


#Find all the merged taxaID
mergeddict = {}
me = open(merged)
for ln in me:
    ln = ln.rstrip()
    line = re.split('\|', ln)
    key = line[0].strip()
    value = line[1].strip()
    mergeddict[key] = value
me.close()

### Read in TaxID from the 16S Microbial Database ####
tx = open(args.dir + "/" + dbname + ".ACCtaxID")
print("\n>> Loading " + dbname + " TaxID list...")
g2txdic = {}
for ln in tx:
    ln = ln.rstrip().split("|")
    if ln[1] in mergeddict:
        ln[1] = mergeddict[ln[1]]
    g2txdic[ln[0]] = ln[1]
tx.close()

glst = g2txdic.values()
print(">> Loading " + dbname + " TaxID list...DONE!\nTotal",
      len(set(glst)), "TaxID to fetch!")

# Find all the nodes in the subset of database
print(">> Loading nodes.dmp...This will take 10-15 minutes, please wait!")
nolist = list(set(glst))
allnode = {}
nnolist = len(nolist) + 1

while len(nolist) < nnolist:
    #   nodmp=open('subset.nodes.dmp')
    print(">> 1 >> Open nodes.dmp file!")
    nodmp = open(nodes)
    i = 1
    nnolist = len(nolist)
    for ln in nodmp:
        lne = ln.rstrip()
        line = re.split("\t\|\t|\t\|", lne)
        i += 1
        if i % 100000 == 0:
            print("Scanning nodes.dmp line:", i)
        if line[0] in nolist:
            allnode[line[0]] = line[1:3]
            nolist.remove(line[0])
            if line[1] not in nolist:
                nolist.append(line[1])
    nodmp.close()
    print(">> 2 >> Close nodes.dmp file!")
    print(">> Remaining # TaxID to look for:", len(nolist))
#       break
print(">> Loading nodes.dmp...DONE!")

# Find all the names in the subset of database
nalist = allnode.keys()
namesdic = {}
print(">> Loading names.dmp...\nThis may take 20-30 minutes. Please wait!")
na = open(names)
id = 0
for ln in na:
    lne = ln.rstrip()
    line = re.split("\t\|\t|\t\|", lne)
    if line[3] == "scientific name" and line[0] in nalist:
        namesdic[line[0]] = line[1:4]
    id += 1
    if id % 50000 == 0:
        print(">> ", id, "names recorded!")
na.close()
print(">> Loading names.dmp...DONE!")

os.system("rm -f " + args.dir + "/" + dbname + ".ACC.taxonomy")
subtx = open(args.dir + "/" + dbname + ".ACC.taxonomy", 'w')
print(">> Generating a subset of taxonomy file.")
for gi, ti in g2txdic.items():
    subtx.write(gi + "\t" + lineage(ti, allnode, namesdic, nolist) + "\n")
subtx.close()
print(">> Taxonomy file generated!")

# setup_tax_db(taxdb)
