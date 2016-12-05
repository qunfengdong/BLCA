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
        print "   This is the utility script to format 16S Microbial Database from NCBI before running the BLCA taxonomy profiling. This could be used for other subsets of NCBI formatted database for blast too.\n"
	print 'Usage: python '+sys.argv[0]+'\n'
        print "Arguments:\n - Optional:"
        print "\t-d\t\tThe database link that you want to download from and format. Default: ftp://ftp.ncbi.nlm.nih.gov/blast/db/16SMicrobial.tar.gz.\n\t-t\t\tThe taxonomy database link from NCBI. Default: ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip.\n\t-u\t\tThe taxdb from NCBI. Default: ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz \n - Other:"
	print "\t-h\t\tShow program usage and quit"	

db='ftp://ftp.ncbi.nlm.nih.gov/blast/db/16SMicrobial.tar.gz'
tax='ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip'
taxdb='ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz'

nodes='db/nodes.dmp'
names='db/names.dmp'

#### Get options ####
opts, args=getopt.getopt(sys.argv[1:],"d:t:u:h",['Database','fsadb','taxdb','help'])
for o,a in opts:
        if o == "-d":
                db=a
	elif o == "-t":
		tax=a
	elif o == "-u":
		taxdb=a
        elif o in ('-h','--help'):
                print usage()
                sys.exit()
        else:
                assert False, 'unhandle option'

################ Function ##############
                
def setup_tax_db(url):
	os.system("wget -P ./db/ "+url)
	os.system("tar zxvf ./db/taxdb.tar.gz -C db/")
	print "taxdb downloaded!\n######################\n>> Please run the following before the next script:\nexport BLASTDB="+os.getcwd()+"/db/"+"\n######################"

def lineage(tid,allnode,namesdic,nonexist):
	levels=["superkingdom","phylum","class","order","family","genus","subspecies","species"]
        if tid not in nonexist:
                r=str()
        #       print allnode[tid]
                name=namesdic[tid][0]
                pa, rank=allnode[tid]
                if rank in levels:
                	r+=":".join([rank,name])+";"
        #       print pa, rank
                if pa in nonexist:
                        r+=""
                else:
                        r+=lineage(pa,allnode,namesdic,nonexist)
        #       print r
        #       print "----"
                return r
        else:
                return ""

def get_taxdb(url):
	os.system("wget -P db/ "+url)
	os.system("unzip db/taxdmp.zip -d db/")
	print ">> NCBI Taxonomy Database downloaded!"

def gi_taxid_nucl(fsa,gi_taxid_dmp):
        fa=open(fsa)
        gilist={}
        for ln in fa:
                lne=ln.rstrip()
                if lne[0] == ">":
                        gilist[lne.split()[0].replace(">","")]=lne.split("|")[1]
        fa.close()
        glst=gilist.values()
        tagi={}
        dmp=open(gi_taxid_dmp)
        for l in dmp:
                ln=l.rstrip()
                line=ln.split("\t")
                if line[0] in glst:
                        tagi[line[0]]=line[1]
                        glst.remove(line[0])
        dmp.close()
        outmap=open(fsa+".taxid_mapping",'wa')
        for k,v in gilist.items():
                outmap.write(k+"\t"+tagi[v]+"\n")
        outmap.close()

def check_program(prgname):
	'''Check whether a program has been installed and put in the PATH'''
	import os
	path=os.popen("which "+prgname).read().rstrip()
	if len(path) > 0 and os.path.exists(path):
		print prgname+" is located in your PATH!"
	else:
		print "ERROR: "+prgname+" is NOT in your PATH, please set up "+prgname+"!"
		sys.exit(1)


def get_db_gitaxid(url):
	dblist=re.split(r'\.|\/+',url)[-3]
	check_program("blastdbcmd")
	os.system("wget -P db/ "+url)
	os.system("tar zxvf db/"+dblist+".tar.gz -C db/")
	subprocess.call(['blastdbcmd','-entry','all','-outfmt','%a|%T','-db','db/'+dblist,'-out','db/'+dblist+'.ACCtaxID'])
	print ">> Accession and TaxIDs from "+dblist+" are extracted!!"

#####################################################################
####	End of functions
####	Working scripts
####
#####################################################################

os.system("mkdir db")
get_taxdb(tax)

get_db_gitaxid(db)

dbname=re.split(r'\.|\/+',db)[-3]

### Read in TaxID from the 16S Microbial Database ####
tx=open("db/"+dbname+".ACCtaxID")
print "\n>> Loading "+dbname+" TaxID list..."
g2txdic={}
for ln in tx:
        ln=ln.rstrip().split("|")
        g2txdic[ln[0]]=ln[1]
tx.close()

glst=g2txdic.values()
print ">> Loading "+dbname+" TaxID list...DONE!\nTotal",len(set(glst)),"TaxID to fetch!"

### Find all the nodes in the subset of database
print ">> Loading nodes.dmp...This will take 10-15 minutes, please wait!"
nolist=list(set(glst))
allnode={}
nnolist=len(nolist)+1
while len(nolist) < nnolist:
#	nodmp=open('subset.nodes.dmp')
	print ">> 1 >> Open nodes.dmp file!"
	nodmp=open(nodes)
	i=1	
#       print nolist
	nnolist=len(nolist)
        for ln in nodmp:
                lne=ln.rstrip()
                line=re.split("\t\|\t|\t\|",lne)
		i+=1
		if i%100000 == 0:
			print "Scanning nodes.dmp line:", i
                if line[0] in nolist:
                        allnode[line[0]]=line[1:3]
                        nolist.remove(line[0])
                        if line[1] not in nolist:
                                nolist.append(line[1])
        nodmp.close()
	print ">> 2 >> Close nodes.dmp file!"
	print ">> Remaining # TaxID to look for:",len(nolist)
#       break
print ">> Loading nodes.dmp...DONE!"

### Find all the names in the subset of database
nalist=allnode.keys()
namesdic={}
print ">> Loading names.dmp...\nThis may take 20-30 minutes. Please wait!"
na=open(names)
id=0
for ln in na:
        lne=ln.rstrip()
        line=re.split("\t\|\t|\t\|",lne)
        if line[3]=="scientific name" and line[0] in nalist:
                namesdic[line[0]]=line[1:4]
		id+=1
		if id%2000==0:
			print ">> ", id, "names recorded!"
na.close()
print ">> Loading names.dmp...DONE!"

os.system("rm -f db/"+dbname+".ACC.taxonomy")
subtx=open("db/"+dbname+".ACC.taxonomy",'wa')
print ">> Generating a subset of taxonomy file."
for gi, ti in g2txdic.items():
        subtx.write(gi+"\t"+lineage(ti,allnode,namesdic,nolist)+"\n")
subtx.close()
print ">> Taxonomy file generated!"

# setup_tax_db(taxdb)
