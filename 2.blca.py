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
import os
# import argparse

'''


'''

def usage():
	print "\n<< Bayesian-based LCA taxonomic classification method >>\n\n   Please make sure the following softwares are in your PATH:\n\t1.muscle (http://www.drive5.com/muscle/downloads.htm)\n\t2.ncbi-blast suite (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)\n\t3.Biopython should be installed locally.\n"
	print 'Usage: python '+sys.argv[0]+' -i <fasta file> [option]\n'
	print " \nArguments:\n - Required:"
	print "\t-i\t\tInput fasta file.\n - Taxonomy Profiling Options [filitering of hits]:"
	print "\t-n\t\tNumber of times to bootstrap. Default: 100"
	print "\t-t\t\tExtra number of nucleotides to include at the beginning and end of the hits. Default: 10"
	print "\t-d\t\tProportion of hits to include from top hit. Default: 0.1 [0-1]"
	print "\t-e\t\tMinimum evalue to include for blastn. Default: 0.1"
	print "\t-a\t\tMinimum bitscore to include for blastn hits. Default: 100"
	print "\t-c\t\tMinimum coverage to include. Default: 0.95 [0-1]"
	print "\t-b\t\tMinimum identity score to include. Default: 95 [0-100]"
	print "\t-r\t\tReference Taxonomy file for the Database. Default: db/16SMicrobial.taxID.taxonomy"
	print "\t-q\t\tRefernece blast database. Default: db/16SMicrobial"
	print "\t-o\t\tOutput file name. Default: <fasta>.blca.out\n - Alignment Options:"
	print "\t-m\t\tAlignment match score. Default: 1"
	print "\t-f\t\tAlignment mismatch penalty. Default: -2.5"
	print "\t-g\t\tAlignment gap penalty. Default: -2\n - Other:"
	print "\t-h\t\tShow program usage and quit" 

iset=float(95.0)
eset=float(0.1)
nper=100
gap=10
topper=float(0.1)
cvrset=float(0.95)
ngap=-2
match=1
mismatch=-2.5
bset=float(100)
tax='./db/16SMicrobial.taxID.taxonomy'
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
	subprocess.call(['blastn','-query',fsa,'-evalue',str(eset),'-dust','no','-soft_masking','false','-db',db,'-num_threads','4','-outfmt','6 std score sstrand slen sseq staxids sscinames qlen','-out',fsa+'.blastn'])
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


################################################################
##
## 	Running Script Start
##
################################################################

## check taxdb
check_taxdb()

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
gi2taxid={}
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
	sname=line[17]
	if (sstart-gap)<1:
		sstart=1
	else:
		sstart-=gap
	if (send+gap)>slen:
		send=slen
	else:
		send+=gap
	coverage=float(line[3])/float(line[18])
	tids=line[16].split(";")
	if evalue < eset and identity > iset and tids[0] != '1161941' and coverage >= cvrset and line[11]>bset:
		giinfo[line[1]+":"+str(sstart)+"-"+str(send)]=[sstart,send,sstrand,bltscore,evalue,identity,sname]
               	gi2taxid[line[1].split("|")[1]]=tids[0]
		if line[0].replace("|","_") not in sublist:
			topsc=float(line[11])
			qtosdic[line[0].replace("|","_")]=[line[1]+":"+str(sstart)+"-"+str(send)]
			sublist.append(line[0].replace("|","_"))
		else:
			bitsc=float(line[11])
#			print ">>query:",line[0],">>topscore:",topsc,">>bitscore:", bitsc
			if (topsc-bitsc)/topsc < topper: 
				qtosdic[line[0].replace("|","_")].append(line[1]+":"+str(sstart)+"-"+str(send))

b.close()
print "> 1 > Read in blast output!"
#print qtosdic

# print "> 2 > Writing GI to TAXID file!!"
# os.system("rm "+fsa+".GItaxID")
# gtout=open(fsa+".GItaxID",'wa')
# for k, v in gi2taxid.items():
#         gtout.write(str(k)+"\t"+str(v)+"\n")
# gtout.close()
# print "> 2 > GI to TAXID file outputted!!"

tx=open(tax)
tidtax={}
print "> 3 > Read in taxonomy information!"
for l in tx:
	l=l.rstrip()
	l=l.strip(";")
	lne=l.split("\t")
	if len(lne)==2:
		tidtax[lne[0]]=dict( x.split(":",1) for x in lne[1].split(";"))
tx.close()

f=open(fsa)
fsaln=f.readlines()
f.close()
fsadic={}
for k in range(0,len(fsaln),2):
	name=fsaln[k].rstrip().split(" ")[0].replace(">","").replace("|","_")
	seq=fsaln[k+1].rstrip()
	fsadic[name]=seq	
#print "> 4 > Fasta file read in!!"


levels=["superkingdom","phylum","class","order","family","genus","species"]
for k,v in qtosdic.items():
#	os.system("rm "+k+".dblist")
	msafsa=open(k+".dblist",'aw')
	for g in v:
		msafsa.write(g.split("|")[1]+" "+str(giinfo[g][0])+"-"+str(giinfo[g][1])+" "+giinfo[g][2]+"\n")
	msafsa.close()
	os.system("blastdbcmd -db "+db+" -entry_batch "+k+".dblist -outfmt %f > "+k+".hitdb.fsa")
	fifsa=open(k+".hitdb.fsa",'aw')
	fifsa.write(">"+k+"\n"+fsadic[k])
	fifsa.close()
#	os.system("clustalw2 -INFILE="+k+".hitdb.fsa -OUTFILE="+k+".clustalw")
#	os.system("t_coffee "+k+".hitdb.fsa -gapopen 0 -gapext -2 -mode quickaln -outfile "+k+".t_coffee -tree "+k+".t_coffee.tre -output aln")
	os.system("muscle -quiet -clw -in "+k+".hitdb.fsa -out "+k+".muscle")
	alndic=get_dic_from_aln(k+".muscle")
	os.system("rm "+k+".muscle")
        os.system("rm "+k+".hitdb.fsa")
        os.system("rm "+k+".dblist")
	start,end=get_gap_pos(k,alndic)
	trunc_alndic=cut_gap(alndic,start,end)
	orgscore=pairwise_score(trunc_alndic,k,match,mismatch,ngap)
	perdict={}
	pervote={}
	for j in range(nper):
		tmpdic=random_aln_score(trunc_alndic,k,match,mismatch,ngap)
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
	ttlvote=sum(pervote.values())
	for k4,v4 in pervote.items():
		pervote[k4]=v4/ttlvote*100	
#	print "For read "+k+"\tThe original score:\tThe permutated vote:\tAll permutation score[list]:"
#	print "################################"
#	print "@Read_"+k+" Blastn_rawscore\tEvalue\tIdentity\tMSA_orgscore\tVote\tSuperkingdom\tPhylum\tClass\tOrder\tFamily\tGenus"
	hitstax={}
	outout=open(outfile,'aw')
	for k1,v1 in orgscore.items():
		k2=[x for x in giinfo.keys() if k1 in x][0]
		k2tax=tidtax[gi2taxid[k2.split("|")[1]]]
		mistx=list(set(levels)-set(k2tax.keys()))
		for mx in mistx:
			k2tax[mx]="Not Available"		
#		print k2, giinfo[k2][3],giinfo[k2][4]
#		if pervote.has_key(k1):
#			print k1+"\t"+str(giinfo[k2][3])+"\t"+str(giinfo[k2][4])+"\t"+str(giinfo[k2][5])+"\t"+str(v1)+"\t"+str(pervote[k1])+"\t",[perdict.values()[i][k1] for i in range(nper)]
#			print "> table > "+k2+"\t"+str(giinfo[k2][3])+"\t"+str(giinfo[k2][4])+"\t"+str(giinfo[k2][5])+"\t"+str(v1)+"\t"+str(pervote[k1])+"\t"+k2tax['superkingdom']+"\t"+k2tax['phylum']+"\t"+k2tax['class']+"\t"+k2tax['order']+"\t"+k2tax['family']+"\t"+k2tax['genus']
#		else:
#			print k1+"\t"+str(giinfo[k2][3])+"\t"+str(giinfo[k2][4])+"\t"+str(giinfo[k2][5])+"\t"+str(v1)+"\t"+str(0)+"\t",[perdict.values()[i][k1] for i in range(nper)]
#			print "> table > "+k2+"\t"+str(giinfo[k2][3])+"\t"+str(giinfo[k2][4])+"\t"+str(giinfo[k2][5])+"\t"+str(v1)+"\t"+str(0)+"\t"+k2tax['superkingdom']+"\t"+k2tax['phylum']+"\t"+k2tax['class']+"\t"+k2tax['order']+"\t"+k2tax['family']+"\t"+k2tax['genus']
		if not pervote.has_key(k1):
			pervote.update({k1:0})
		hitstax[k2]=[k2tax,pervote[k1],orgscore[k1]]
	voten=[y[1] for y in hitstax.values()]
	prbn=[z[2] for z in hitstax.values()]
#	print "---------------------------------------"
#	print "> result > @"+k+" Taxonomy Prediction with Confidence."
	outout.write(k+"\t")
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
#		print "> result > "+le+":"
#		for k5,v5 in lexsum.items():
#			print "> result > \t"+k5+":\t"+str(v5)+"\t("+str(prosum[k5])+")"
#		for k5,v5 in lexsum.items():
#			if float(v5)>=cutoff:
#				print le+":"+k5+"; "+str(v5)+";",
		outout.write(le+":"+max(lexsum,key=lexsum.get)+";"+str(max(lexsum.values()))+";")
	relgi = [x for x in giinfo.keys() if max(pervote,key=pervote.get) in x][0]
	outout.write("strain:"+giinfo[relgi][6]+";"+str(max(pervote.values()))+";\n")
	outout.close()
print ">> Taxonomy file generated!!"
