#!/usr/bin/env python

__author__ = "Xiang Gao, Huaiying Lin, Qunfeng Dong"
__copyright__ = "Copyright 2016, Bayesian-based LCA taxonomic classification"
__credits__ = ["Xiang Gao", "Huaiying Lin", "Qunfeng Dong"]
__license__ = "GPL"
__version__ = "2.1"
__maintainer__ = "Huaiying Lin"
__email__ = "hlin2@luc.edu"
__status__ = "Python3"

import os
from sys import exit, stderr
from collections import defaultdict

try:
    from Bio import AlignIO, SeqIO
except ImportError:
    stderr.write("Error! BioPython is not detected!\n")
    exit(1)

import random
import argparse

import time

startTime = time.time()

'''
BLCA Core annotation tool

'''


def float_1range(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]" % (x,))
    return x

def float_100range(x):
    x = float(x)
    if x < 0.0 or x > 100.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 100.0]" % (x,))
    return x


def check_taxdb():
    ''' Check whether BLASTDB environmental variable has been setup'''
    if 'BLASTDB' not in os.environ.keys():
        print(
            "ERROR: taxdb location has not been set up as the environmental variable BLASTDB. Please set up as \n\texport BLASTDB=/location/of/taxdb.bti/and/taxdb.btd/")
        exit(1)


def check_program(prgname):
    '''Check whether a program has been installed and put in the PATH'''
    import os
    path = os.popen("which " + prgname).read().rstrip()
    if len(path) > 0 and os.path.exists(path):
        print(prgname + " is located in your PATH!")
    else:
        print("ERROR: " + prgname + " is NOT in your PATH, please set up " + prgname + "!")
        exit(1)


def run_blastn(fsa, eset, nsub, db, proc):
    '''Running blastn with customized output format'''
    check_program("blastn")
    import subprocess
    print("> > Running blast!!")
    # subprocess.call(['blastn','-query',fsa,'-evalue',str(eset),'-dust','no','-soft_masking','false','-db',db,'-num_threads',str(proc),'-outfmt','6 std score sstrand slen qlen','-max_target_seqs',str(nsub),'-out',fsa+'.blastn'])
    subprocess.call(
        ['blastn', '-query', fsa, '-evalue', str(eset), '-perc_identity', '90', '-soft_masking', 'false', '-db', db,
         '-num_threads', str(proc), '-outfmt', '6 std score sstrand slen qlen', '-max_target_seqs', str(nsub), '-out',
         fsa + '.blastn'])
    print("> > Blastn Finished!!")


def get_dic_from_aln(aln):
    '''Read in alignment and convert it into a dictionary'''
    alignment = AlignIO.read(aln, "clustal")
    alndic = {}
    for r in alignment:
        alndic[r.id] = list(r.seq)
    return alndic


def pairwise_score(alnmtdic, query, match, mismatch, ngap):
    '''Calculate pairwise alignment score given a query'''
    nt = ["A", "C", "T", "G", "g", "a", "c", "t"]
    hitscore = defaultdict(float)
    for k, v in alnmtdic.items():
        if k != query:
            for i in range(len(v)):
                if (alnmtdic[query][i] in nt) and (v[i] in nt) and (alnmtdic[query][i] == v[i]):
                    hitscore[k] += match
                elif (alnmtdic[query][i] not in nt) and (v[i] not in nt) and (alnmtdic[query][i] == v[i]):
                    hitscore[k] += 0
                elif ((alnmtdic[query][i] not in nt) or (v[i] not in nt)) and (alnmtdic[query][i] != v[i]):
                    hitscore[k] += mismatch
                elif (alnmtdic[query][i] in nt) and (v[i] in nt) and (alnmtdic[query][i] != v[i]):
                    hitscore[k] += ngap
    total = float(sum(hitscore.values()))
    if total <= 0:
        total = 1
    for k, v in hitscore.items():
        hitscore[k] = v / total
    return hitscore


def random_aln_score(alnmtdic, query, match, mismatch, ngap):
    '''Randomize the alignment, and calculate the score'''
    nt = ["A", "C", 'T', "G", "g", "a", "c", "t"]
    idx = []
    alnlist = list(alnmtdic.values())
    for i in range(len(alnlist[0])):  # 'dict_values' object does not support indexing in python 3
        idx.append(random.choice(range(len(alnlist[0]))))
    hitscore = defaultdict(float)
    for k, v in alnmtdic.items():
        if k != query:
            for i in idx:
                if (alnmtdic[query][i] in nt) and (v[i] in nt) and (alnmtdic[query][i] == v[i]):
                    hitscore[k] += match
                elif (alnmtdic[query][i] not in nt) and (v[i] not in nt) and (alnmtdic[query][i] == v[i]):
                    hitscore[k] += 0
                elif ((alnmtdic[query][i] not in nt) or (v[i] not in nt)) and (alnmtdic[query][i] != v[i]):
                    hitscore[k] += mismatch
                elif (alnmtdic[query][i] in nt) and (v[i] in nt) and (alnmtdic[query][i] != v[i]):
                    hitscore[k] += ngap
    return hitscore


def get_gap_pos(query, alnmtdic):
    '''Get the gap position in the alignment'''
    for i in range(len(alnmtdic[query])):
        if alnmtdic[query][i] != "-":
            start = i
            break
    for i in range(len(alnmtdic[query]) - 1, 0, -1):
        if alnmtdic[query][i] != "-":
            end = i
            break
    return start, end


def cut_gap(alndic, start, end):
    '''Given a start and end gap position, truncate the alignmnet'''
    trunc_alndic = {}
    for k_truc, v_truc in alndic.items():
        trunc_alndic[k_truc] = v_truc[start:end]
    return trunc_alndic


def read_tax_acc(taxfile, IDlen):
    acctax = {}
    print(">  > Reading in taxonomy information! ....")
    with open(taxfile) as tx:
        for l in tx:
            lne = l.rstrip().strip(";").split("\t")
            if len(lne) == 2:
                if len(lne[0]) > IDlen:
                    print("ERROR: Your reference sequence ID length longer than 32, please shorten your ID length and try again!! ")
                    exit(1)
                else:
                    acctax[lne[0].split(".")[0]] = dict(x.split(":", 1) for x in lne[1].split(";"))
    if len(acctax) == 0:
        print("ERROR: No taxonmy record was read in. Please examine the format of your taxonomy file according to the instruction on github.")
        exit(1)
    return acctax


##### parser ######

parser = argparse.ArgumentParser(description=
                                 ''' << Bayesian-based LCA taxonomic classification method >>\n\n   Please make sure the following softwares are in your PATH:
    1.muscle (http://www.drive5.com/muscle/downloads.htm), muscle should be the program's name.
    2.ncbi-blast suite (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
    3.clustalo (http://www.clustal.org/omega/), clustalo should be the program's name.
    4.Biopython should be installed locally.''',
                                 epilog="No warrenty comes with this script. Author: hlin2@luc.edu. \nAny suggestions or bugs report are welcomed.",
                                 add_help=False, formatter_class=argparse.RawTextHelpFormatter)
##### Required arguments #####
required = parser.add_argument_group('required arguments')
required.add_argument("-i", "--fsa", help="Input fasta file", type=str, required=True)
##### Taxonomy filtering arguments #####
taxoptions = parser.add_argument_group('taxonomy profiling options [filtering of hits]')
taxoptions.add_argument("-x", "--skipblast", help="skip blastn. Default: blastn is not skipped", action='store_true')
taxoptions.add_argument("-n", "--nper", help="number of times to bootstrap. Default: 100", type=int, default=100)
taxoptions.add_argument("-j", "--nsub", help="maximum number of subjects to include for each query reads. Default: 50",
                        type=int, default=50)
taxoptions.add_argument("-d", "--topper", help="proportion of hits to include from top hit. Default: 0.1 [0-1]",
                        type=float_1range, default=0.1)
taxoptions.add_argument("-e", "--eset", help="minimum evalue to include for blastn. Default: 0.1", type=float,
                        default=0.1)
taxoptions.add_argument("-b", "--bset", help="minimum bitscore to include for blastn hits. Default: 100", type=int,
                        default=100)
taxoptions.add_argument("-c", "--cvrset", help="minimum coverage to include. Default: 0.85 [0-1]", type=float_1range,
                        default=0.85)
taxoptions.add_argument("--iset", help="minimum identity score to include. Default: 90 [0-100]", type=float_100range,
                        default=90.0)
##### Alignment control arguments #####
alignoptions = parser.add_argument_group('alignment control arguments')
alignoptions.add_argument("-a", "--align", default='clustalo',
                          help="alignment tool: clustal omega or muscle. Default: clustalo", type=str,
                          choices=["clustalo", "muscle"])
alignoptions.add_argument("-m", "--match", default=1.0, help="alignment match score. Default: 1", type=float)
alignoptions.add_argument("-f", "--mismatch", default=-2.5, help="alignment mismatch penalty. Default: -2.5",
                          type=float)
alignoptions.add_argument("-g", "--ngap", default=-2.0, help="alignment gap penalty. Default: -2", type=float)
##### Other arguments #####
optional = parser.add_argument_group('other arguments')
optional.add_argument("-r", "--tax", default='./db/16S_ribosomal_RNA.ACC.taxonomy',
                      help="reference taxonomy file for the Database. Default: db/16S_ribosomal_RNA.taxonomy", type=str)
optional.add_argument("-q", "--db", default='./db/16S_ribosomal_RNA',
                      help="refernece blast database. Default: db/16S_ribosomal_RNA", type=str)
optional.add_argument("-t", "--gap", default=10,
                      help="extra number of nucleotides to include at the beginning and end of the hits. Default: 10",
                      type=int)
optional.add_argument("-o", "--outfile", help="output file name. Default: <fasta>.blca.out", type=str)
optional.add_argument("-p", "--proc", default=2,
                      help="how many processors are used in blastn step. Default: 2 processors", type=int)
optional.add_argument("-h", "--help", help="show this help message and exit", action="help")
##### parse arguments #####
args = parser.parse_args()
##### Output file name ####
if not args.outfile:
    args.outfile = args.fsa + ".blca.out"

################################################################
##
##  Running Script Start
##
################################################################

## check taxdb

## check whether blastdbcmd is located in the path
##check_program("blastdbcmd")

## check whether alignment tool is located in the path
if args.align == "clustalo":
    check_program("clustalo")
    IDlenallow = 100
elif args.align == "muscle":
    check_program("muscle")
    IDlenallow = 32

### read in input fasta file ###
fsadic = {}
with open(args.fsa) as f:
    for r in SeqIO.parse(f, "fasta"):
        if len(r.id) > IDlenallow:
            print("ERROR: Your query ID length is longer than " + str(
                IDlenallow) + ", please shorten your ID length and try again!!")
            exit(1)
        else:
            fsadic[r.id] = str(r.seq)
print(">  > Fasta file read in!")

### read in pre-formatted lineage information ###
acc2tax = read_tax_acc(args.tax, IDlen=IDlenallow)

## Run blastn and output fas.blastn output
if not args.skipblast:
    run_blastn(args.fsa, args.eset, args.nsub, args.db, args.proc)
elif args.skipblast:
    print(">  > blastn is skipped")

## read in blastn output
print(">  > read in blast file...")
b = open(args.fsa + '.blastn', 'r')
print(">  > blastn file opened")
qtosdic = {}
giinfo = {}
## place holder for tracking input reads
sublist = set()  # major speed boost when list replaced by set
for ln in b:
    ln = ln.rstrip()
    line = ln.split("\t")
    sstart = int(line[8])
    send = int(line[9])
    sstrand = line[13]
    slen = int(line[14])
    bltscore = float(line[12])
    evalue = float(line[10])
    identity = float(line[2])
    ### extend sstart and send by gap ###
    if (sstart - args.gap) < 1:
        sstart = 1
    else:
        sstart -= args.gap
    if (send + args.gap) > slen:
        send = slen
    else:
        send += args.gap
    ### calculate coverage ###
    coverage = float(int(line[7]) - int(line[6]) + 1) / float(line[15])
    ### filter hits ###
    if evalue <= args.eset and identity >= args.iset and coverage >= args.cvrset and float(
            line[11]) >= args.bset:  # minimum values to include for blastn hits can be equal to user defined values
        giinfo[line[1] + ":" + str(sstart) + "-" + str(send)] = [sstart, send, sstrand, bltscore, evalue, identity]
        if line[0].replace("|", "_") not in sublist:
            subcount = 1
            topsc = float(line[11])
            qtosdic[line[0].replace("|", "_")] = [line[1] + ":" + str(sstart) + "-" + str(send)]
            sublist.add(line[0].replace("|", "_"))
        else:
            bitsc = float(line[11])
            if (topsc - bitsc) / topsc < args.topper and subcount <= args.nsub:
                qtosdic[line[0].replace("|", "_")].append(line[1] + ":" + str(sstart) + "-" + str(send))
                subcount += 1
b.close()
print(">  > blast output read in")

print(">  > Start aligning reads...")
# print("Cut offs")
# print("evalue:",args.eset,"identity:", args.iset, "coverage:", args.cvrset)
# print qtosdic.values()
os.system("rm -f " + args.outfile)
outout = open(args.outfile, 'w+')
levels = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
for seqn in fsadic.keys():
    # for k1,v1 in qtosdic.items():
    k1 = seqn
    if k1 in qtosdic and k1 not in acc2tax:
        v1 = qtosdic[k1]
        os.system("rm -f " + k1 + ".dblist")
        os.system("rm -f " + k1 + ".hitdb.fsa")
        ### Get all the hits list belong to the same query ###
        msafsa = open(k1 + ".dblist", 'w')
        for g in v1:
            # if blastn returns hits on minus strand, reverse sstart and send positions in -range arg
            if giinfo[g][1] < giinfo[g][0]:
                msafsa.write(g.split(":")[0] + " " + str(giinfo[g][1]) + "-" + str(giinfo[g][0]) + " " + giinfo[g][2] + "\n")
            else:
                msafsa.write(g.split(":")[0] + " " + str(giinfo[g][0]) + "-" + str(giinfo[g][1]) + " " + giinfo[g][2] + "\n")
        msafsa.close()
        os.system("blastdbcmd -db " + args.db + " -entry_batch " + k1 + ".dblist -outfmt %f > " + k1 + ".hitdb.fsa")
        ### Add query fasta sequence to extracted hit fasta ###
        fifsa = open(k1 + ".hitdb.fsa", 'a')
        fifsa.write(">" + k1 + "\n" + fsadic[k1])
        fifsa.close()
        os.system("rm " + k1 + ".dblist")
        ### Run alignment ###
        if args.align == "clustalo":
            os.system("clustalo --threads=" + str(args.proc) + " --outfmt=clu -i " + k1 + ".hitdb.fsa -o " + k1 + ".aln")
        if args.align == "muscle":
            os.system("muscle -quiet -clw -in " + k1 + ".hitdb.fsa -out " + k1 + ".aln")
        alndic = get_dic_from_aln(k1 + ".aln")
        os.system("rm " + k1 + ".hitdb.fsa")
        os.system("rm " + k1 + ".aln")
        # print "Processing:",k1
        ### get gap position and truncate the alignment###
        start, end = get_gap_pos(k1, alndic)
        trunc_alndic = cut_gap(alndic, start, end)
        orgscore = pairwise_score(trunc_alndic, k1, args.match, args.mismatch, args.ngap)
        ### start bootstrap ###
        perdict = {}  # record alignmet score for each iteration
        pervote = {}  # record vote after nper bootstrap
        ### If any equal score, average the vote ###
        for j in range(args.nper):
            tmpdic = random_aln_score(trunc_alndic, k1, args.match, args.mismatch, args.ngap)
            # print("iter: ", j, ":", tmpdic) # for checking purpose
            perdict[j] = tmpdic
            mx = max(tmpdic.values())
            mxk = [k3 for k3, v3 in tmpdic.items() if v3 == mx]
            ### check if there are multiple max hits ###
            if len(mxk) == 1:
                if max(tmpdic, key=tmpdic.get) in pervote:
                    pervote[max(tmpdic, key=tmpdic.get)] += float(1)
                else:
                    pervote[max(tmpdic, key=tmpdic.get)] = float(1)
            else:
                por = float(1) / len(mxk)
                for h in mxk:
                    if h in pervote:
                        pervote[h] += por
                    else:
                        pervote[h] = por
        ### normalize vote by total votes ###
        ttlvote = sum(pervote.values())
        for k4, v4 in pervote.items():
            pervote[k4] = v4 / ttlvote * 100
        ### hit taxnomy parsing ###
        hitstax = {}
        for k5, v5 in orgscore.items():
            shortk5 = k5.split(":")[0].split(".")[0]
            if shortk5 in acc2tax:
                k2tax = acc2tax[shortk5]
                mistx = list(set(levels) - set(k2tax.keys()))
                for mis in mistx:
                    k2tax[mis] = "NA"
                if k5 not in pervote:
                    pervote.update({k5: 0})
                hitstax[k5] = [k2tax, pervote[k5], orgscore[k5]]
        voten = [y[1] for y in hitstax.values()]
        prbn = [z[2] for z in hitstax.values()]
        outout.write(k1 + "\t")
        for le in levels:
            lex = [x[0][le] for x in hitstax.values()]
            lexsum = {}
            prosum = {}
            for d in range(len(lex)):
                # print voten[d]
                if lex[d] in lexsum:
                    lexsum[lex[d]] += voten[d]
                    prosum[lex[d]] += prbn[d]
                else:
                    lexsum[lex[d]] = voten[d]
                    prosum[lex[d]] = prbn[d]
            outout.write(le + ":" + max(lexsum, key=lexsum.get) + ";" + str(max(lexsum.values())) + ";")
        outout.write("\n")
    elif k1 in acc2tax:
        print("[WARNING] Your sequence " + k1 + " has the same ID as the reference database! Please correct it!")
        print("...Skipping sequence " + k1 + " ......")
        outout.write(k1 + "\tSkipped\n")
    else:
        outout.write(k1 + "\tUnclassified\n")

outout.close()
print(">  > Taxonomy file generated!!")

endTime = time.time()
print("Time elapsed: %i minutes" % (round((endTime - startTime) / 60)))
