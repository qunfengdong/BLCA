#!/usr/bin/env python

__author__ = "Xiang Gao, Huaiying Lin, Qunfeng Dong"
__copyright__ = "Copyright 2016, Bayesian-based LCA taxonomic classification"
__credits__ = ["Xiang Gao", "Huaiying Lin", "Qunfeng Dong"]
__license__ = "GPL"
__version__ = "1.2"
__maintainer__ = "Huaiying Lin"
__email__ = "hlin2@luc.edu"
__status__ = "Python3"

import sys
import os
import subprocess
import argparse

try:
    from Bio import AlignIO,SeqIO
except ImportError:
    sys.stderr.write("Error! BioPython is not detected!\n")
    sys.exit(1)

'''
Get a subset of lineage information from GreenGene Microbial database.
'''

##### parser ######

parser = argparse.ArgumentParser(description=
                                 ''' << Bayesian-based LCA taxonomic classification method >>\n\n   Please make sure the following softwares are in your PATH:
                                     1.muscle (http://www.drive5.com/muscle/downloads.htm), muscle should be the program's name.
                                     2.ncbi-blast suite (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
                                     3.clustalo (http://www.clustal.org/omega/), clustalo should be the program's name.
                                     4.Biopython should be installed locally.
                            
                                     This is the utility script to format Greengene Database before running the BLCA taxonomy profiling.
                                     >> Please first download the Greengenes fasta and taxonomy files from https://greengenes.secondgenome.com/?prefix=downloads/greengenes_database/gg_13_5/.''',
                                 epilog="No warrenty comes with this script. Author: hlin2@luc.edu. \nAny suggestions or bugs report are welcomed.",
                                 add_help=False, formatter_class=argparse.RawTextHelpFormatter)
##### Other arguments #####
optional = parser.add_argument_group('optional arguments')
optional.add_argument("--dir", default='gg',
                      help="The local directory name where you want to store the formatted database. Default: gg",
                      type=str)
optional.add_argument("--ggfasta", default='gg_13_5.fasta.gz',
                      help="The GreenGene database fasta file. Default: gg_13_5.fasta.gz", type=str)
optional.add_argument("--ggtax", default='gg_13_5_taxonomy.txt.gz',
                      help="The GreenGene database taxonomy file. Default: gg_13_5_taxonomy.txt.gz", type=str)
optional.add_argument("-t","--fulltax", action='store_true',
                      help="Extract a subset of GreenGene with only reads with full taxonomy information. This could take a while.")
optional.add_argument("-h", "--help", help="show this help message and exit", action="help")
##### parse arguments #####
args = parser.parse_args()

cwd = os.getcwd()
args.ggfasta = os.path.join(cwd, args.ggfasta)
args.ggtax = os.path.join(cwd, args.ggtax)
args.dir = os.path.join(cwd, args.dir)

################ Function ##############

def check_program(prgname):
    '''Check whether a program has been installed and put in the PATH'''
    path = os.popen("which " + prgname).read().rstrip()
    if len(path) > 0 and os.path.exists(path):
        print(prgname + " is located in your PATH!")
    else:
        print("ERROR: ", prgname, " is NOT in your PATH, please set up ", "!")
        sys.exit(1)

def format_ext_gg_taxfile(taxfile):
    ''' Read in Greengenes taxonomy file, and format it and write out into BLCA compatiable format '''
    leveldic = {'c': 'class', 'g': 'genus', 'k': 'superkingdom', 'f': 'family', 'o': 'order', 'p': 'phylum',
                's': 'species'}
    taxin = open(taxfile)
    fsqName = []
    outtaxfile = os.path.splitext(taxfile)[0] + ".fullTax.taxonomy"
    if os.path.isfile(outtaxfile):
        os.remove(outtaxfile)
    taxout=open(outtaxfile,'w')
    for l in taxin:
        ln = l.rstrip().replace(" ", "").split("\t")
        if len(ln) == 2:
            tmp = dict(x.split("__") for x in ln[1].split(";"))
            naN = sum([ tmp[k] == "" for k in tmp ])
            if naN == 0:
                fsqName.append(ln[0])
                taxout.write(ln[0] + "\t")
                oldnames = list(tmp.keys())
                for k in oldnames:
                    tmp[leveldic[k]] = tmp.pop(k)
                for key, val in tmp.items():
                    taxout.write(key + ':' + val + ";")
                taxout.write("\n")
    taxin.close()
    taxout.close()
    print(os.path.basename(taxfile), "has been formatted and outputted!")
    return(fsqName)

def format_gg_taxfile(taxfile):
    ''' Read in Greengenes taxonomy file, and format it and write out into BLCA compatiable format '''
    leveldic={'c':'class','g':'genus','k':'superkingdom','f':'family','o':'order','p':'phylum','s':'species'}
    taxin=open(taxfile)
    outtaxfile = os.path.splitext(taxfile)[0] + ".taxonomy"
    if os.path.isfile(outtaxfile):
        os.remove(outtaxfile)
    taxout=open(outtaxfile,'w')
    for l in taxin:
        ln=l.rstrip().replace(" ","").split("\t")
        if len(ln) == 2:
            taxout.write(ln[0] + "\t")
            tmp = dict(x.split("__") for x in ln[1].split(";"))
            oldnames = list(tmp.keys())
            for k in oldnames:
                tmp[leveldic[k]]=tmp.pop(k)
            for key, val in tmp.items():
                taxout.write(key+':'+val+";")
            taxout.write("\n")
    taxin.close()
    taxout.close()
    print(os.path.basename(taxfile), "has been formatted and outputted!")

def extract_faseq(fsa, seqlist):
    ''' Extract seqs according to taxonomy info'''
    outfsa = os.path.splitext(fsa)[0] + ".fullTax.fasta"
    outhandler = open(outfsa,"w")
    record_dict = SeqIO.index(fsa, "fasta")
    for ID in seqlist:
        SeqIO.write(record_dict[ID], outhandler, "fasta")
    outhandler.close()
    print("Full taxonomy sequences have been extracted from",fsa,"!")
    return(outfsa)

def make_blastdb(dbfsa):
    ''' Format downloaded Greengenes sequence fasta files into blast compatible format '''
    check_program('makeblastdb')
    dbname = os.path.splitext(dbfsa)[0]
    subprocess.call(['makeblastdb', '-dbtype', "nucl", '-in', dbfsa, '-parse_seqids', '-out', dbname])
    print(os.path.basename(dbfsa), "has been formatted!")

def ungz(file, folder):
    fname = os.path.splitext(os.path.basename(file))[0]
    os.system("gzip -dc " + file + " > " + folder + "/" + fname)
    print(file, "has be unzipped!")
    return(os.path.join(folder,fname))


######## MAIN ###########

if __name__ == "__main__":
    os.system('rm -fr ' + args.dir)
    os.system('mkdir ' + args.dir)

    dbfsa = ungz(args.ggfasta, args.dir)
    taxfile = ungz(args.ggtax, args.dir)

    if args.fulltax:
        validseqs = format_ext_gg_taxfile(taxfile)
        newdbfsa = extract_faseq(dbfsa,validseqs)
        make_blastdb(newdbfsa)
    else:
        format_gg_taxfile(taxfile)
        make_blastdb(dbfsa)