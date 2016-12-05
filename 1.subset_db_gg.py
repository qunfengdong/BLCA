#!/usr/bin/env python

import sys
import math
import glob
import re
import os
import getopt
import subprocess

__author__ = "Xiang Gao, Huaiying Lin, Qunfeng Dong"
__copyright__ = "Copyright 2016, Bayesian-based LCA taxonomic classification"
__credits__ = ["Xiang Gao", "Huaiying Lin","Qunfeng Dong"]
__license__ = "GPL"
__version__ = "1.2"
__maintainer__ = "Huaiying Lin"
__email__ = "ying.eddi2008@gmail.com"
__status__ = "Production"

'''

Get a subset of lineage information from NCBI 16S Microbial database.

'''

def usage():
        print "\n<< Bayesian-based LCA taxonomic classification method >>\n\n   Please make sure the following softwares are in your PATH:\n\t1.muscle (http://www.drive5.com/muscle/downloads.htm), muscle should be the program's name.\n\t2.ncbi-blast suite (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)\n\t3.Biopython is installed locally.\n"
        print "   This is the utility script to format Greengene Database before running the BLCA taxonomy profiling.\n"
        print 'Usage: python '+sys.argv[0]+'\n'
        print "Arguments:\n - Optional:"
        print "\t-d\t\tThe database file that you want to format. Default: gg_13_5.fasta.gz.\n\t-t\t\tThe taxonomy database link from Greengenes. Default: gg_13_5_taxonomy.txt.gz. \n - Other:"
        print "\t-h\t\tShow program usage and quit"

dbfsafile='gg_13_5.fasta.gz'
dbtaxfile='gg_13_5_taxonomy.txt.gz'


#### Get options ####
opts, args=getopt.getopt(sys.argv[1:],"d:t:u:h",['Database','fsadb','taxdb','help'])
for o,a in opts:
        if o == "-d":
                dbfsafile=a
        elif o == "-t":
                dbtaxfile=a
        elif o in ('-h','--help'):
                print usage()
                sys.exit()
        else:
                assert False, 'unhandle option'

################ Function ##############

def check_program(prgname):
        '''Check whether a program has been installed and put in the PATH'''
        import os
        path=os.popen("which "+prgname).read().rstrip()
        if len(path) > 0 and os.path.exists(path):
                print prgname+" is located in your PATH!"
        else:
                print "ERROR: "+prgname+" is NOT in your PATH, please set up "+prgname+"!"
                sys.exit(1)


def format_gg_taxfile(taxfile,folder):
	''' Read in Greengenes taxonomy file, and format it and write out into BLCA compatiable format '''
	leveldic={'c':'class','g':'genus','k':'superkingdom','f':'family','o':'order','p':'phylum','s':'species'}	
	taxin=open(folder+'/'+taxfile)
	outtaxname=taxfile.split(".")
	outtaxfile=''.join(outtaxname[:len(outtaxname)-1:])+".taxonomy"
	if os.path.isfile(outtaxfile):
		os.remove(outtaxfile)
	taxout=open(folder+'/'+outtaxfile,'aw')
	for l in taxin:
		ln=l.rstrip().replace(" ","").split("\t")
		taxout.write(ln[0]+"\t")
		tmp=dict(x.split("__") for x in ln[1].split(";"))
		for k in tmp.keys():
			tmp[leveldic[k]]=tmp.pop(k)
		for key, val in tmp.items():
			taxout.write(key+':'+val+";")
		taxout.write("\n")
	taxin.close()
	taxout.close()
	print ">> "+taxfile+" has been formatted and outputted!"

def make_blastdb(dbfsa,folder):
	''' Format downloaded Greengenes sequence fasta files into blast compatiable format '''
	check_program('makeblastdb')
	dbname=folder+'/'+dbfsa.rstrip('.fasta')
	subprocess.call(['makeblastdb','-dbtype',"nucl",'-in',folder+'/'+dbfsa,'-parse_seqids','-out',dbname])
	print dbfsa+" has been formatted!"

def ungz(file,folder):
	fname=file.rstrip('.gz')
	os.system("gzip -dc "+file+" > "+folder+"/"+fname)
	print fname+" has been unzipped!"

######## MAIN ###########

foldername='gg'

os.system('rm -fr '+foldername)
os.system('mkdir '+foldername)

ungz(dbfsafile,foldername)
ungz(dbtaxfile,foldername)

taxfile=dbtaxfile.rstrip('.gz')
dbfsa=dbfsafile.rstrip('.gz')

format_gg_taxfile(taxfile,foldername)
make_blastdb(dbfsa,foldername)


