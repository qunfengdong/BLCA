#!/usr/bin/env python

__author__ = "Xiang Gao, Huaiying Lin, Qunfeng Dong"
__copyright__ = "Copyright 2016, Bayesian-based LCA taxonomic classification"
__credits__ = ["Xiang Gao", "Huaiying Lin","Qunfeng Dong"]
__license__ = "GPL"
__version__ = "1.2"
__maintainer__ = "Huaiying Lin"
__email__ = "hlin2@luc.edu"
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
					  help="The local directory name where you want to store the formatted database. Default: gg", type=str)
optional.add_argument("--ggfasta", default='gg_13_5.fasta.gz',
					  help="The GreenGene database fasta file. Default: gg_13_5.fasta.gz", type=str)
optional.add_argument("--ggtax", default='gg_13_5_taxonomy.txt.gz',
					  help="The GreenGene database taxonomy file. Default: gg_13_5_taxonomy.txt.gz", type=str)
optional.add_argument("-h","--help",help="show this help message and exit",action="help")
##### parse arguments #####
args = parser.parse_args()

################ Function ##############

def check_program(prgname):
	'''Check whether a program has been installed and put in the PATH'''
	path=os.popen("which "+prgname).read().rstrip()
	if len(path) > 0 and os.path.exists(path):
		print(prgname+" is located in your PATH!")
	else:
		print("ERROR: "+prgname+" is NOT in your PATH, please set up "+prgname+"!")
		sys.exit(1)


def format_gg_taxfile(taxfile,folder):
	''' Read in Greengenes taxonomy file, and format it and write out into BLCA compatiable format '''
	leveldic={'c':'class','g':'genus','k':'superkingdom','f':'family','o':'order','p':'phylum','s':'species'}	
	taxin=open(folder+'/'+taxfile)
	outtaxname=taxfile.split(".")
	outtaxfile=''.join(outtaxname[:len(outtaxname)-1:])+".taxonomy"
	if os.path.isfile(outtaxfile):
		os.remove(outtaxfile)
	taxout=open(folder+'/'+outtaxfile,'w')
	for l in taxin:
		ln=l.rstrip().replace(" ","").split("\t")
		taxout.write(ln[0]+"\t")
		if len(ln) > 2:
			tmp=dict(x.split("__") for x in ln[1].split(";"))
			for k in tmp:
				tmp[leveldic[k]]=tmp.pop(k)
			for key, val in tmp.items():
				taxout.write(key+':'+val+";")
			taxout.write("\n")
	taxin.close()
	taxout.close()
	print(">> "+taxfile+" has been formatted and outputted!")

def make_blastdb(dbfsa,folder):
	''' Format downloaded Greengenes sequence fasta files into blast compatible format '''
	check_program('makeblastdb')
	dbname=folder+'/'+dbfsa.rstrip('.fasta')
	subprocess.call(['makeblastdb','-dbtype',"nucl",'-in',folder+'/'+dbfsa,'-parse_seqids','-out',dbname])
	print(dbfsa+" has been formatted!")

def ungz(file,folder):
	fname=file.rstrip('.gz')
	os.system("gzip -dc "+file+" > "+folder+"/"+fname)
	print(fname+" has been unzipped!")

######## MAIN ###########

if __name__ == "__main__":
	os.system('rm -fr '+args.dir)
	os.system('mkdir '+args.dir)

	ungz(args.ggfasta,args.dir)
	ungz(args.ggtax,args.dir)

	dbfsa=args.ggfasta.rstrip('.gz')
	taxfile=args.ggtax.rstrip('.gz')

	format_gg_taxfile(taxfile,args.dir)
	make_blastdb(dbfsa,args.dir)


