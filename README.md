Bayesian-based LCA taxonomic classification method
--------------------------------------------------
Bayesian LCA-based Taxonomic Classification Method (BLCA) is a Bayesian-based method that provides a solid probabilistic basis for evaluating the taxonomic assignments for the query sequences with bootstrap confidence scores, which is based on Bayesian posterior probability that quantitatively weigh each database hit sequence according to its similarity to the query sequence - the more similar database hit sequence to the query, the more its contribution to the taxonomic assignment of the query. 

We implemented the above algorithm as a simple python script here.

## Update
* **Sep 16 2019** Minor update of 2.blca_main.py to fix minus strand range from Blastn output (credit to [Carter Hoffman](mailto:hoffmanc@ohsu.edu)).
* **Jul 8 2019** Minor update of 2.blca_main.py to fix clustalo's compatibility issue with blast 2.9.0.
* **Jun 3 2019** Minor update of 2.blca_main.py to fix the hidden 100% confidence score bug.
* **May 9 2019** Minor update to 1.subset_db_gg.py to include a new function to extract only sequences with full taxonomy information.
* **Feb 26 2019 update** One utility script (generate_abundance_table.py) for merging multiple BLCA output is available in the utils folder.
* **Feb 21 2019 update** The entire package has been updated to python 3.
* **Nov 15 2018 update** Thanks to Kristjan's contribution, now we incorporated the use of clustalo as alignment software. Also, now BLCA main program is based on **python 3**.
* **May 11 2017 update** to be compatible for the latest blastn v2.5 and added a new parameter -j to limit the accepted hits number to 50. After another round of testing, we've decided to change the default value of coverage and identify filter to 0.80 and 90 respectively.

## Important Note -- Please do read
* BLCA has migrated to **Python 3**. If you'd like to use python2.7, please install from release (https://github.com/qunfengdong/BLCA/releases).
* BLCA currently is compatible with **blast 2.9.0+**, please make sure you have the latest blast: version 2.9.0 or above. 
* There should **NOT** be any "|" (pipe) present in the sequence ID of input fasta, database fasta and taxonomy files.

## Prerequisities
* Python 3
* Linux
* Biopython

### **The following programs should be in your PATH:**

* BLAST binary (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/)
* clustalo (http://www.clustal.org/omega/)
* MUSCLE (http://www.drive5.com/muscle/downloads.htm)

## Citation

[A Bayesian Taxonomic Classification Method for 16S rRNA Gene Sequences with Improved Species-level Accuracy](https://www.ncbi.nlm.nih.gov/pubmed/28486927). *Xiang Gao; Huaiying Lin; Kashi Revanna; Qunfeng Dong* BMC Bioinformatics 2017 May 10;18(1):247.

## Install
To check out the source code, go to https://github.com/qunfengdong/BLCA. To obtain the scripts and example fasta files, do the following:

```
$ git clone https://github.com/qunfengdong/BLCA.git
```

After the github repository is cloned, you will find a folder named BLCA. All the scripts and example data files will be included in it. **It is highly recommended to run your own analysis inside this directory (BLCA), meaning you should have your fasta files moved to here, so you don't have to change the default database directory.**

## Quick start

We do not include a pre-compiled database with this release, so the first step is to build a taxonomy database from the NCBI 16S microbial database. We achieve this by using script _1.subset_db_acc.py_ (or 1.subset_db_gg.py). After the database is built and stored on your local machine, you will supply the location of the taxonomy output file (16SMicrobial.taxID.taxonomy) from the last step along with your input fasta file (test.fasta) to _2.blca_main.py_, then you will get a blca output as test.fasta.blca.out.

## Getting started

### Step 1
* To compile, subset the 16S Microbial database. Please run:
```
$ python 1.subset_db_acc.py
```
More options available:
```
$ python 1.subset_db_acc.py -h

usage: 1.subset_db_acc.py [--dir DIR] [-d DATABASE] [--taxdmp TAXDMP]
                          [--taxdb TAXDB] [-h]

 << Bayesian-based LCA taxonomic classification method >>

   Please make sure the following softwares are in your PATH:
		 1.muscle (http://www.drive5.com/muscle/downloads.htm), muscle should be the program's name.
		 2.ncbi-blast suite (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
		 3.clustalo (http://www.clustal.org/omega/), clustalo should be the program's name.
		 4.Biopython should be installed locally.

optional arguments:
  --dir DIR             The local directory name where you want to store the formatted database. Default: db
  -d DATABASE, --database DATABASE
                        The database link that you want to download from and format. Default: ftp://ftp.ncbi.nlm.nih.gov/blast/db/16SMicrobial.tar.gz
  --taxdmp TAXDMP       The taxonomy database dmp link from NCBI. Default: ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip
  --taxdb TAXDB         The taxonomy database db link from NCBI. Default: ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
  -h, --help            show this help message and exit

No warrenty comes with this script. Author: hlin2@luc.edu. 
Any suggestions or bugs report are welcomed.
```
During the process of setting up the database, NCBI's 16SMicrobial.tar.gz file, and taxdmp.zip will be downloaded into a default folder: ./db/, and uncompressed. 16SMicrobial.ACC.taxonomy under the ./db directory is the taxonomy file should be supplied to the 2.blca_main.py as the database. 

### Alternative Step 1
* To format GreenGenes database, first you have to download the Greengenes fasta and taxonomy files from https://greengenes.secondgenome.com/?prefix=downloads/greengenes_database/gg_13_5/. The files you need are gg_13_5.fasta.gz and gg_13_5_taxonomy.txt.gz. After you make sure you download the targeted two files under BLCA folder, please run:
```
$ python 1.subset_db_gg.py
```
This script will unzip the downloaded files and create a new folder called "gg" to store all needed information.

More options available:
```
$ python 1.subset_db_gg.py -h
usage: 1.subset_db_gg.py [--dir DIR] [--ggfasta GGFASTA] [--ggtax GGTAX] [-t]
                         [-h]

 << Bayesian-based LCA taxonomic classification method >>

   Please make sure the following softwares are in your PATH:
      1.muscle (http://www.drive5.com/muscle/downloads.htm), muscle should be the program's name.
      2.ncbi-blast suite (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
      3.clustalo (http://www.clustal.org/omega/), clustalo should be the program's name.
      4.Biopython should be installed locally.

      This is the utility script to format Greengene Database before running the BLCA taxonomy profiling.
      >> Please first download the Greengenes fasta and taxonomy files from https://greengenes.secondgenome.com/?prefix=downloads/greengenes_database/gg_13_5/.

optional arguments:
  --dir DIR          The local directory name where you want to store the formatted database. Default: gg
  --ggfasta GGFASTA  The GreenGene database fasta file. Default: gg_13_5.fasta.gz
  --ggtax GGTAX      The GreenGene database taxonomy file. Default: gg_13_5_taxonomy.txt.gz
  -t, --fulltax      Extract a subset of GreenGene with only reads with full taxonomy information. This could take a while.
  -h, --help         show this help message and exit

No warrenty comes with this script. Author: hlin2@luc.edu. 
Any suggestions or bugs report are welcomed.
```
### SILVA LSU database (Credit to Dr. Daniel Swan)

Thanks to Dr. Swan's personal effort to build a BLCA-compatible blastn-database and taxonomy file for SILVA LSU database. 

* Download the pre-compiled database at this [link](https://drive.google.com/drive/folders/1t0TzC08y7_LyglsdihaXu27oWr7PiKLe).

* After you've downloaded the SILVA LSU fasta and taxonomy file, you will need to format the fasta file as the following:

``` 
$ makeblastdb -in SILVA_132_LSURef_tax_silva_BLCAparsed.fasta -dbtype nucl -parse_seqids -out SILVA_132_LSURef_tax_silva_BLCAparsed.fasta
```
* Then you can follow the instructions in the [Training your own database](#training-your-own-database) section.

### Split input fasta (Optional)
* If you have a big fasta file, and you want to run BLCA in "parallel", you can use [this python package](https://pypi.python.org/pypi/pyfasta/#command-line-interface) to split fasta sequences into multiple parts, then run BLCA on each individual part.

### Step 2 
* Run your analysis with the compiled database. Please run:
``` 
$ python 2.blca_main.py -i test.fasta
```
If you are running your analysis somewhere else other than in the BLCA directory, please do the following:
```
$ python /location/to/2.blca_main.py -i test.fasta -r /location/to/your/database/16SMicrobial.ACC.taxonomy -q /location/to/your/database/16SMicrobial
```
If you are using the Greengene database as your reference, please do the following:
```
$ python /location/to/2.blca_main.py -i test.fasta -r gg/gg_13_5_taxonomy.taxonomy -q gg/gg_13_5
```

More options are the following:
```
$ python 2.blca_main.py -h

usage: 2.blca_main.py -i FSA [-x] [-n NPER] [-j NSUB] [-d TOPPER] [-e ESET]
                      [-b BSET] [-c CVRSET] [--iset ISET] [-a ALIGN]
                      [-m MATCH] [-f MISMATCH] [-g NGAP] [-r TAX] [-q DB]
                      [-t GAP] [-o OUTFILE] [-p PROC] [-h]

 << Bayesian-based LCA taxonomic classification method >>

   Please make sure the following softwares are in your PATH:
    1.muscle (http://www.drive5.com/muscle/downloads.htm), muscle should be the program's name.
    2.ncbi-blast suite (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
    3.clustalo (http://www.clustal.org/omega/), clustalo should be the program's name.
    4.Biopython should be installed locally.

required arguments:
  -i FSA, --fsa FSA     Input fasta file

taxonomy profiling options [filtering of hits]:
  -x, --skipblast       skip blastn. Default: blastn is not skipped
  -n NPER, --nper NPER  number of times to bootstrap. Default: 100
  -j NSUB, --nsub NSUB  maximum number of subjects to include for each query reads. Default: 50
  -d TOPPER, --topper TOPPER
                        proportion of hits to include from top hit. Default: 0.1 [0-1]
  -e ESET, --eset ESET  minimum evalue to include for blastn. Default: 0.1
  -b BSET, --bset BSET  minimum bitscore to include for blastn hits. Default: 100
  -c CVRSET, --cvrset CVRSET
                        minimum coverage to include. Default: 0.85 [0-1]
  --iset ISET           minimum identity score to include. Default: 90 [0-100]

alignment control arguments:
  -a ALIGN, --align ALIGN
                        alignment tool: clustal omega or muscle. Default: clustalo
  -m MATCH, --match MATCH
                        alignment match score. Default: 1
  -f MISMATCH, --mismatch MISMATCH
                        alignment mismatch penalty. Default: -2.5
  -g NGAP, --ngap NGAP  alignment gap penalty. Default: -2

other arguments:
  -r TAX, --tax TAX     reference taxonomy file for the Database. Default: db/16SMicrobial.ACC.taxonomy
  -q DB, --db DB        refernece blast database. Default: db/16SMicrobial
  -t GAP, --gap GAP     extra number of nucleotides to include at the beginning and end of the hits. Default: 10
  -o OUTFILE, --outfile OUTFILE
                        output file name. Default: <fasta>.blca.out
  -p PROC, --proc PROC  how many processors are used in blastn step. Default: 2 processors
  -h, --help            show this help message and exit

No warrenty comes with this script. Author: hlin2@luc.edu. 
Any suggestions or bugs report are welcomed.
```

## Output
* A text file with sequence id in the first column, and taxonomy annotation with confidence scores after each level of annotaion (superkingdom, phylum, class, order, family, genus, species).
* Any reads that do not have a classification will be recorded as "Unclassified".
* There could be cases having the confidence score showing while there is no taxonomy assignment at genus/species level. It is due to the lack of taxonomy information in the database.

### Example output file:
```
seq94	superkingdom:Bacteria;100.0;phylum:Firmicutes;100.0;class:Clostridia;100.0;order:Clostridiales;100.0;family:Lachnospiraceae;100.0;genus:Lachnoclostridium;100.0;species:[Clostridium] symbiosum;100.0;
seq89   superkingdom:Bacteria;100.0;phylum:Proteobacteria;100.0;class:Gammaproteobacteria;100.0;order:Aeromonadales;57.4166666667;family:Aeromonadaceae;57.4166666667;genus:;57.4166666667;species:;100.0;
seq87   superkingdom:Bacteria;100.0;phylum:Firmicutes;100.0;class:Clostridia;100.0;order:Clostridiales;100.0;family:Ruminococcaceae;100.0;genus:;69.0019047619;species:;100.0;
seq93	superkingdom:Bacteria;100.0;phylum:Actinobacteria;100.0;class:Actinobacteria;100.0;order:Corynebacteriales;100.0;family:Nocardiaceae;100.0;genus:Rhodococcus;100.0;species:Rhodococcus zopfii;99.5;
seq96	Unclassified
```

## Training your own database

* BLCA main script 2.blca_main.py needs 
1. A BLAST formatted library from a fasta file containing sequences of your interest, using makeblastdb, as the following:
```
>NR_117221.1
AGTCGATCGATCGATCATCGCTCTCTAGAGAGAAAACCCGATCGATCGA...
>NR_144700.1
CGCGCGACGAGCAAGCGCAAACGGCAACGCGCGAAACCCGCGAGCGAGA...

$ makeblastdb -in YourDatabase.fasta -dbtype nucl -parse_seqids -out YourDatabase
```

2. A taxonomy file with two columns, sequence ID in fasta file, and its taxonomy from superkingdom to species in the following format (The deliminator between the sequence ID and taxonomy information should be a **tab [\t]**):
```
NR_117221.1     species:Mycobacterium arosiense;genus:Mycobacterium;family:Mycobacteriaceae;order:Corynebacteriales;class:Actinobacteria;phylum:Actinobacteria;superkingdom:Bacteria;
NR_144700.1     species:Virgibacillus massiliensis;genus:Virgibacillus;family:Bacillaceae;order:Bacillales;class:Bacilli;phylum:Firmicutes;superkingdom:Bacteria;
NR_108831.1     species:Bacillus endoradicis;genus:Bacillus;family:Bacillaceae;order:Bacillales;class:Bacilli;phylum:Firmicutes;superkingdom:Bacteria;
NR_113104.1     species:Prevotella enoeca;genus:Prevotella;family:Prevotellaceae;order:Bacteroidales;class:Bacteroidia;phylum:Bacteroidetes;superkingdom:Bacteria;
NR_027573.1     species:Intestinibacter bartlettii;genus:Intestinibacter;family:Peptostreptococcaceae;order:Clostridiales;class:Clostridia;phylum:Firmicutes;superkingdom:Bacteria;
```

3. Run 2.blca_main.py with the formatted database and taxonomy file.
```bash
$ python 2.blca_main.py -i test.fasta -r /location/to/your/database/YourDatabase.taxonomy -q /location/to/your/database/YourDatabase
```

## Version
* Version 2.2 An alternative public release

## Authors
* Dr. Xiang Gao, theoretical conception and algorithm development
* Dr. Qunfeng Dong, algorithm development
* Huaiying Lin, program coding and testing
* Kashi Revanna, program coding and package development

## Error report

Please report any errors or bugs to hlin2@luc.edu.

## License
GNU

## Acknowledgements
* BLAST program: Camacho C., Coulouris G., Avagyan V., Ma N., Papadopoulos J., Bealer K., & Madden T.L. (2008) "BLAST+: architecture and applications." BMC Bioinformatics 10:421.
* Clustal Omega: Sievers, Fabian, Andreas Wilm, David Dineen, Toby J. Gibson, Kevin Karplus, Weizhong Li, Rodrigo Lopez et al. "Fast, scalable generation of high‐quality protein multiple sequence alignments using Clustal Omega." Molecular systems biology 7, no. 1 (2011): 539.
* MUSCLE: Edgar, R.C. (2004) MUSCLE: multiple sequence alignment with high accuracy and high throughput.Nucleic Acids Res. 32(5):1792-1797. doi:10.1093/nar/gkh340
* Biopython: Cock PA, Antao T, Chang JT, Bradman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B and de Hoon MJL (2009) Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25, 1422-1423
