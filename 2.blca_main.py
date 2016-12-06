#!/usr/bin/env python

__author__ = "Xiang Gao, Huaiying Lin, Qunfeng Dong"
__copyright__ = "Copyright 2016, Bayesian-based LCA taxonomic classification"
__credits__ = ["Xiang Gao", "Huaiying Lin","Qunfeng Dong"]
__license__ = "GPL"
__version__ = "0.8"
__maintainer__ = "Huaiying Lin"
__email__ = "ying.eddi2008@gmail.com"
__status__ = "Production"


import sys
import os
import math

try:
	from Bio import AlignIO
except ImportError:
	sys.stderr.write("Error! BioPython is not detected!\n")
	sys.exit(1)

import random
import subprocess
import getopt
# import argparse

'''
BLCA Core annotation tool

'''

def usage():
	print "\n<< Bayesian-based LCA taxonomic classification method >>\n\n   Please make sure the following softwares are in your PATH:\n\t1.muscle (http://www.drive5.com/muscle/downloads.htm), muscle should be the program's name.\n\t2.ncbi-blast suite (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)\n\t3.Biopython should be installed locally.\n"
	print 'Usage: python '+sys.argv[0]+' -i <fasta file> [option]\n'
	print " \nArguments:\n - Required:"
	print "\t-i\t\tInput fasta file.\n - Taxonomy Profiling Options [filtering of hits]:"
	print "\t-n\t\tNumber of times to bootstrap. Default: 100"
	print "\t-t\t\tExtra number of nucleotides to include at the beginning and end of the hits. Default: 10"
	print "\t-d\t\tProportion of hits to include from top hit. Default: 0.1 [0-1]"
	print "\t-e\t\tMinimum evalue to include for blastn. Default: 0.1"
	print "\t-a\t\tMinimum bitscore to include for blastn hits. Default: 100"
	print "\t-c\t\tMinimum coverage to include. Default: 0.95 [0-1]"
	print "\t-b\t\tMinimum identity score to include. Default: 95 [0-100]"
	print "\t-r\t\tReference Taxonomy file for the Database. Default: db/16SMicrobial.ACC.taxonomy"
	print "\t-q\t\tRefernece blast database. Default: db/16SMicrobial"
	print "\t-o\t\tOutput file name. Default: <fasta>.blca.out\n - Alignment Options:"
	print "\t-m\t\tAlignment match score. Default: 1"
	print "\t-f\t\tAlignment mismatch penalty. Default: -2.5"
	print "\t-g\t\tAlignment gap penalty. Default: -2\n - Other:"
	print "\t-h\t\tShow program usage and quit" 

####set up default parameters #######
### bootstrap times ###
nper=100            # number of bootstrap to permute
### Filter hits per query ###
iset=float(95.0)    # identify threshold
eset=float(0.1)     # evalue threshold
gap=10              # number of nucleotides to include at the beginning and end of the hits
topper=float(0.1)   # top percentage to include hits
cvrset=float(0.95)  # coverage to include
bset=float(100)     # minimum bitscore to include
### Alignment options ###
ngap=-2             # gap penalty
match=1             # match score
mismatch=-2.5       # mismatch penalty
### pre-formatted taxonomy and blastn database ###
tax='./db/16SMicrobial.ACC.taxonomy'
db='./db/16SMicrobial'

opts, args=getopt.getopt(sys.argv[1:],"b:c:d:e:f:g:i:lm:n:o:r:q:t:h",['Minimum Identity','Minimum Coverage','Top Proportion','Minimum evale','Alignment Mismatch Penalty','Alignment Gap Penalty','Input File','Long Output','Alignment Match Score','Number of Bootstrap to Run','Output File','Taxonomy File of Database','Reference Blastn Database','Number of nt Length of Sequence','help'])
for o,a in opts:
	if o == "-i":
		fsa=a
		outfile=fsa+".blca.out"
	elif o == "-n":
		nper=float(a)
	elif o == "-t":
		gap=float(a)
	elif o == "-c":
            cvrset=a
	elif o == "-l":
		longout=True
	elif o == "-d":
		topper=float(a)
	elif o == "-e":
		eset=float(a)
	elif o == "-f":
		mismatch=float(a)
	elif o == "-g":
		ngap=float(a)
	elif o == "-q":
		db=a
	elif o == "-m":
		match=float(a)
	elif o == "-b":
		iset=float(a)
	elif o== "-o":
		outfile=a
	elif o== "-r":
		tax=a
	elif o== "-a":
		bset=float(a)	
	elif o in ('-h','--help'):
		print usage()
		sys.exit(1)
	else:
		assert False, 'unhandle option'

def check_taxdb():
	''' Check whether BLASTDB environmental variable has been setup'''
	if 'BLASTDB' not in os.environ.keys():
		print "ERROR: taxdb location has not been set up as the environmental variable BLASTDB. Please set up as \n\texport BLASTDB=/location/of/taxdb.bti/and/taxdb.btd/"
		sys.exit(1)

def check_program(prgname):
	'''Check whether a program has been installed and put in the PATH'''
	import os
	path=os.popen("which "+prgname).read().rstrip()
	if len(path) > 0 and os.path.exists(path):
		print prgname+" is located in your PATH!"
	else:
		print "ERROR: "+prgname+" is NOT in your PATH, please set up "+prgname+"!"
		sys.exit(1)

def run_blastn(fsa):
	'''Running blastn with customized output format'''
	check_program("blastn")
	import subprocess
	print ">> Running blast!!"
	subprocess.call(['blastn','-query',fsa,'-evalue',str(eset),'-dust','no','-soft_masking','false','-db',db,'-num_threads','4','-outfmt','6 std score sstrand slen qlen','-out',fsa+'.blastn'])
	print ">> Blastn Finished!!"

def get_dic_from_aln(aln):
	'''Read in alignment and convert it into a dictionary'''
	alignment=AlignIO.read(aln,"clustal")
	alndic={}
	for r in alignment:
		alndic[r.id]=list(r.seq)
	return alndic

def pairwise_score(alndic,query,match,mismatch,ngap):
	'''Calculate pairwise alignment score given a query'''
	nt=["A","C","T","G","g","a","c","t"]
	hitscore={}
	for k,v in alndic.items():
		if k != query:
			hitscore[k]=0
			for i in range(len(v)):
				if (alndic[query][i] in nt) and (v[i] in nt) and (alndic[query][i] == v[i]):
					hitscore[k]+=float(match)
				elif (alndic[query][i] not in nt) and (v[i] not in nt) and (alndic[query][i] == v[i]):
					hitscore[k]+=float(0)
				elif ((alndic[query][i] not in nt) or (v[i] not in nt)) and (alndic[query][i] != v[i]):
					hitscore[k]+=float(mismatch)
				elif (alndic[query][i] in nt) and (v[i] in nt) and (alndic[query][i] != v[i]):
					hitscore[k]+=float(ngap)		
	total=sum(hitscore.values())
	for k,v in hitscore.items():
		hitscore[k]=v/total
	return hitscore

def random_aln_score(alndic,query,match,mismatch,ngap):
	'''Randomize the alignment, and calculate the score'''
	nt=["A","C","T","G","g","a","c","t"]
	idx=[]
	for i in range(len(alndic.values()[0])):
		idx.append(random.choice(range(len(alndic.values()[0]))))
	hitscore={}
	for k,v in alndic.items():
		if k != query:
			hitscore[k]=0
			for i in idx:
                                if (alndic[query][i] in nt) and (v[i] in nt) and (alndic[query][i] == v[i]):
                                        hitscore[k]+=float(match)
                                elif (alndic[query][i] not in nt) and (v[i] not in nt) and (alndic[query][i] == v[i]):
                                        hitscore[k]+=float(0)
                                elif ((alndic[query][i] not in nt) or (v[i] not in nt)) and (alndic[query][i] != v[i]):
                                        hitscore[k]+=float(mismatch)
                                elif (alndic[query][i] in nt) and (v[i] in nt) and (alndic[query][i] != v[i]):
                                        hitscore[k]+=float(ngap)           
        return hitscore				

def get_gap_pos(query,alndic):
	'''Get the gap position in the alignment'''
	for i in range(len(alndic[query])):
		if alndic[query][i] != "-":
			start=i
			break
	for i in range(len(alndic[query])-1,0,-1):
		if alndic[query][i] != "-":
			end=i
			break
	return start,end

def cut_gap(alndic,start,end):
	'''Given a start and end gap position, truncate the alignmnet'''
	trunc_alndic={}
	for k_truc,v_truc in alndic.items():
		trunc_alndic[k_truc]=v_truc[start:end]
	return trunc_alndic

def read_tax_acc(taxfile):
	tx=open(taxfile)
	acctax={}
	print "> 3 > Read in taxonomy information!"
	for l in tx:
        	lne=l.rstrip().strip(";").split("\t")
		if len(lne)==2:
	                acctax[lne[0].split(".")[0]]=dict( x.split(":",1) for x in lne[1].split(";"))
	tx.close()
	return acctax

################################################################
##
## 	Running Script Start
##
################################################################

## check taxdb
# check_taxdb()

## check whether blastdbcmd is located in the path
check_program("blastdbcmd")

## check whether muscle is located in the path
check_program("muscle")

## Run blastn and output fas.blastn output
run_blastn(fsa)

## read in blastn output
b=open(fsa+'.blastn')
qtosdic={}
giinfo={}
sublist=[]
for ln in b:
	ln=ln.rstrip()
	line=ln.split("\t")
	sstart=int(line[8])
	send=int(line[9])
	sstrand=line[13]
	slen=int(line[14])
	bltscore=float(line[12])
	evalue=float(line[10])
	identity=float(line[2])
    ### extend sstart and send by gap ###
	if (sstart-gap)<1:
		sstart=1
	else:
		sstart-=gap
	if (send+gap)>slen:
		send=slen
	else:
		send+=gap
    ### calculate coverage ###
	coverage=float(line[3])/float(line[15])
    ### filter hits ###
	if evalue < eset and identity > iset and coverage >= cvrset and line[11]>bset:
		giinfo[line[1]+":"+str(sstart)+"-"+str(send)]=[sstart,send,sstrand,bltscore,evalue,identity]
		if line[0].replace("|","_") not in sublist:
			topsc=float(line[11])
			qtosdic[line[0].replace("|","_")]=[line[1]+":"+str(sstart)+"-"+str(send)]
			sublist.append(line[0].replace("|","_"))
		else:
			bitsc=float(line[11])
			if (topsc-bitsc)/topsc < topper: 
				qtosdic[line[0].replace("|","_")].append(line[1]+":"+str(sstart)+"-"+str(send))

b.close()
print "> 1 > Read in blast output!"

### read in pre-formatted lineage information ###
acc2tax=read_tax_acc(tax)

### read in input fasta file ###
f=open(fsa)
fsaln=f.readlines()
f.close()
fsadic={}
for k in range(0,len(fsaln),2):
	name=fsaln[k].rstrip().split(" ")[0].replace(">","").replace("|","_")
	seq=fsaln[k+1].rstrip()
	fsadic[name]=seq	
#print "> 4 > Fasta file read in!!"


levels=["superkingdom","phylum","class","order","family","genus","subspecies","species"]
for k1,v1 in qtosdic.items():
	os.system("rm -f "+k1+".dblist")
	os.system("rm -f "+k1+".hitdb.fsa")
	### Get all the hits list belong to the same query ###
	msafsa=open(k1+".dblist",'aw')
	for g in v1:
		msafsa.write(g.split(":")[0]+" "+str(giinfo[g][0])+"-"+str(giinfo[g][1])+" "+giinfo[g][2]+"\n")
	msafsa.close()
	os.system("blastdbcmd -db "+db+" -entry_batch "+k1+".dblist -outfmt %f > "+k1+".hitdb.fsa")
	### Add query fasta sequence to extracted hit fasta ###
	fifsa=open(k1+".hitdb.fsa",'aw')
	fifsa.write(">"+k1+"\n"+fsadic[k1])
	fifsa.close()
	os.system("rm "+k1+".dblist")
	### Run muscle ###
	os.system("muscle -quiet -clw -in "+k1+".hitdb.fsa -out "+k1+".muscle")
	alndic=get_dic_from_aln(k1+".muscle")
	os.system("rm "+k1+".hitdb.fsa")
	os.system("rm "+k1+".muscle")
    	### get gap position and truncate the alignment###
	start,end=get_gap_pos(k1,alndic)
	trunc_alndic=cut_gap(alndic,start,end)
	orgscore=pairwise_score(trunc_alndic,k1,match,mismatch,ngap)
	### start bootstrap ###
	perdict={}	# record alignmet score for each iteration
	pervote={}	# record vote after nper bootstrap
	### If any equal score, average the vote ###
	for j in range(nper):
		tmpdic=random_aln_score(trunc_alndic,k1,match,mismatch,ngap)
		perdict[j]=tmpdic
		mx=max(tmpdic.values())
		mxk=[k3 for k3,v3 in tmpdic.items() if v3 == mx]
		if len(mxk)==1:
			if pervote.has_key(max(tmpdic,key=tmpdic.get)):
				pervote[max(tmpdic,key=tmpdic.get)]+=float(1)
			else:
				pervote[max(tmpdic,key=tmpdic.get)]=float(1)
		else:
			por=float(1)/len(mxk)
			for h in mxk:
				if pervote.has_key(h):
                                	pervote[h]+=por
                        	else:
                                	pervote[h]=por
	### normalize vote by total votes ###
	ttlvote=sum(pervote.values())
	for k4,v4 in pervote.items():
		pervote[k4]=v4/ttlvote*100	
	### 
	hitstax={}
	outout=open(outfile,'aw')
	for k5,v5 in orgscore.items():
		shortk5=k5.split(".")[0]
		if acc2tax.has_key(shortk5):
			k2tax=acc2tax[shortk5]
			mistx=list(set(levels)-set(k2tax.keys()))
			for mis in mistx:
				k2tax[mis]="Not Available"
			if not pervote.has_key(k5):
				pervote.update({k5:0})
			hitstax[k5]=[k2tax,pervote[k5],orgscore[k5]]
	voten=[y[1] for y in hitstax.values()]
	prbn=[z[2] for z in hitstax.values()]
	outout.write(k1+"\t")
	for le in levels:
		lex=[x[0][le] for x in hitstax.values()]
		lexsum={}
		prosum={}
		for d in range(len(lex)):
#			print voten[d]
			if lexsum.has_key(lex[d]):
				lexsum[lex[d]]+=voten[d]
				prosum[lex[d]]+=prbn[d]
			else:
				lexsum[lex[d]]=voten[d]
				prosum[lex[d]]=prbn[d]
		outout.write(le+":"+max(lexsum,key=lexsum.get)+";"+str(max(lexsum.values()))+";")
	outout.write("\n")
	outout.close()
print ">> Taxonomy file generated!!"
