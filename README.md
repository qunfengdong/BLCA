Bayesian-based LCA taxonomic classification method
--------------------------------------------------
Bayesian LCA-based Taxonomic Classification Method (BLCA) is a Bayesian-based method that provides a solid probabilistic basis for evaluating the taxonomic assignments for the query sequences with bootstrap confidence scores, which is based on Bayesian posterior probability that quantitatively weigh each database hit sequence according to its similarity to the query sequence - the more similar database hit sequence to the query, the more its contribution to the taxonomic assignment of the query. 

We implemented the above algorithm as a simple python script here.

## Update
* **May 11 2017 update** to be compatiable for the latest blastn v2.5 and added a new parameter -j to limit the accepted hits number to 50. After another round of testing, we've decided to change the default value of coverage and identify filter to 0.80 and 90 respectively.

## Important Note -- Please do read
* BLCA currently is only compatiable with **blast 2.5.0+**, please make sure you have blast 2.5.0. Any error messages raised as the following could be a blastn version issue.
```shell
outout.write(le+":"+max(lexsum,key=lexsum.get)+";"+str(max(lexsum.values()))+";")
ValueError: max() arg is an empty sequence
```
* Read ID length should be no more than **28** due to limited ID space in muscle output. 

## Prerequisities
* Python 2.7
* Linux
* Biopython

### **The following programs should be in your PATH:**

* BLAST binary (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.5.0/)
* MUSCLE (http://www.drive5.com/muscle/downloads.htm)

## Citation

[A Bayesian Taxonomic Classification Method for 16S rRNA Gene Sequences with Improved Species-level Accuracy](https://www.ncbi.nlm.nih.gov/pubmed/28486927). *Xiang Gao; Huaiying Lin; Kashi Revanna; Qunfeng Dong* BMC Bioinformatics 2017 May 10;18(1):247.

## Install
To check out the source code, go to https://github.com/qunfengdong/BLCA. To obtain the scripts and example fasta files, do the following:

```shell
$ git clone https://github.com/qunfengdong/BLCA.git
```

After the github repository is cloned, you will find a folder named BLCA. All the scripts and example data files will be included in it. **It is highly recommended to run your own analysis inside this directory (BLCA), meaning you should have your fasta files moved to here, so you don't have to change the default database directory.**

## Quick start

We do not include a pre-compiled database with this release, so the first step is to build a taxonomy database from the NCBI 16S microbial database. We achieve this by using script _1.subset_db_acc.py_ (or 1.subset_db_gg.py). After the database is built and stored on your local machine, you will supply the loction of the taxonomy output file (16SMicrobial.taxID.taxonomy) from the last step along with your input fasta file (test.fasta) to _2.blca_main.py_, then you will get a blca output as test.fasta.blca.out.

## Getting started

### Step 1
* To compile, subset the 16S Microbial database. Please run:
```
$ python 1.subset_db_acc.py
```
More options available:
```
$ python 1.subset_db_acc.py -h

<< Bayesian-based LCA taxonomic classification method >>

   Please make sure the following softwares are in your PATH:
	1.muscle (http://www.drive5.com/muscle/downloads.htm), muscle should be the program's name.
	2.ncbi-blast suite (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
	3.Biopython is installed locally.

   This is the utility script to format 16S Microbial Database from NCBI before running the BLCA taxonomy profiling. This could be used for other subsets of NCBI formatted database for blast too.

Usage: python 1.subset_db_acc.py

Arguments:
 - Optional:
	-d		The database link that you want to download from and format. Default: ftp://ftp.ncbi.nlm.nih.gov/blast/db/16SMicrobial.tar.gz.
	-t		The taxonomy database link from NCBI. Default: ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip.
	-u		The taxdb from NCBI. Default: ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz 
 - Other:
	-h		Show program usage and quit
```
During the process of setting up the database, NCBI's 16SMicrobial.tar.gz file, and taxdmp.zip will be downloaded into a default folder: ./db/, and uncompressed. 16SMicrobial.ACC.taxonomy under the ./db directory is the taxonomy file should be supplied to the 2.blca_main.py as the database. 

### Alternative Step 1
* To format Greengenes database, first you have to download the Greengenes fasta and taxonomy files from http://greengenes.secondgenome.com/downloads/database/13_5. The files you need are gg_13_5.fasta.gz and gg_13_5_taxonomy.txt.gz. After you make sure you download the targeted two files under BLCA folder, please run:
```
$ python 1.subset_db_gg.py
```
This script will unzip the downloaded files and create a new folder called "gg" to store all needed information.

More options available:
```
$ python 1.subset_db_gg.py -h

<< Bayesian-based LCA taxonomic classification method >>

   Please make sure the following softwares are in your PATH:
	1.muscle (http://www.drive5.com/muscle/downloads.htm), muscle should be the program's name.
	2.ncbi-blast suite (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
	3.Biopython is installed locally.

   This is the utility script to format Greengene Database before running the BLCA taxonomy profiling.

Usage: python 1.subset_db_gg.py

Arguments:
 - Optional:
	-d		The database file that you want to format. Default: gg_13_5.fasta.gz.
	-t		The taxonomy database link from Greengenes. Default: gg_13_5_taxonomy.txt.gz. 
 - Other:
	-h		Show program usage and quit
```
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

<< Bayesian-based LCA taxonomic classification method >>

   Please make sure the following softwares are in your PATH:
	1.muscle (http://www.drive5.com/muscle/downloads.htm), muscle should be the program's name.
	2.ncbi-blast suite (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
	3.Biopython should be installed locally.

Usage: python 2.blca_main.py -i <fasta file> [option]

 
Arguments:
 - Required:
	-i		Input fasta file.
 - Taxonomy Profiling Options [filtering of hits]:
	-n		Number of times to bootstrap. Default: 100
	-j		Maximum number of subjects to include for each query reads. Default: 50
	-d		Proportion of hits to include from top hit. Default: 0.1 [0-1]
	-e		Minimum evalue to include for blastn. Default: 0.1
	-a		Minimum bitscore to include for blastn hits. Default: 100
	-c		Minimum coverage to include. Default: 0.85 [0-1]
	-b		Minimum identity score to include. Default: 90 [0-100]
	-r		Reference Taxonomy file for the Database. Default: db/16SMicrobial.ACC.taxonomy
	-q		Refernece blast database. Default: db/16SMicrobial
	-o		Output file name. Default: <fasta>.blca.out
 - Alignment Options:
	-m		Alignment match score. Default: 1
	-f		Alignment mismatch penalty. Default: -2.5
	-g		Alignment gap penalty. Default: -2
 - Other:
	-t		Extra number of nucleotides to include at the beginning and end of the hits. Default: 10
	-h		Show program usage and quit
```

## Output
* A text file with sequence id in the first column, and taxonomy annotation with confidence scores after each level of annotaion (superkingdom, phylum, class, order, family, genus, species).
* Any reads that do not have a classification will be recorded as "unclassified".

### Example output file:
```
seq94	superkingdom:Bacteria;100.0;phylum:Firmicutes;100.0;class:Clostridia;100.0;order:Clostridiales;100.0;family:Lachnospiraceae;100.0;genus:Lachnoclostridium;100.0;species:[Clostridium] symbiosum;100.0;
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

makeblastdb -in YourDatabase.fasta -dbtype nucl -parse_seqids -out YourDatabase
```

2. A taxonomy file with two columns, sequence ID in fasta file, and its taxonomy from superkingdom to species in the following format:
```
NR_117221.1     species:Mycobacterium arosiense;genus:Mycobacterium;family:Mycobacteriaceae;order:Corynebacteriales;class:Actinobacteria;phylum:Actinobacteria;superkingdom:Bacteria;
NR_144700.1     species:Virgibacillus massiliensis;genus:Virgibacillus;family:Bacillaceae;order:Bacillales;class:Bacilli;phylum:Firmicutes;superkingdom:Bacteria;
NR_108831.1     species:Bacillus endoradicis;genus:Bacillus;family:Bacillaceae;order:Bacillales;class:Bacilli;phylum:Firmicutes;superkingdom:Bacteria;
NR_113104.1     species:Prevotella enoeca;genus:Prevotella;family:Prevotellaceae;order:Bacteroidales;class:Bacteroidia;phylum:Bacteroidetes;superkingdom:Bacteria;
NR_027573.1     species:Intestinibacter bartlettii;genus:Intestinibacter;family:Peptostreptococcaceae;order:Clostridiales;class:Clostridia;phylum:Firmicutes;superkingdom:Bacteria;
```

3. Run 2.blca_main.py with the formatted database and taxonomy file.
```python
python 2.blca_main.py -i test.fasta -r /location/to/your/database/YourDatabase.taxonomy -q /location/to/your/database/YourDatabase
```

## Version
* Version 2.1 An alternative public release

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
* MUSCLE: Edgar, R.C. (2004) MUSCLE: multiple sequence alignment with high accuracy and high throughput.Nucleic Acids Res. 32(5):1792-1797. doi:10.1093/nar/gkh340
* Biopython: Cock PA, Antao T, Chang JT, Bradman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B and de Hoon MJL (2009) Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25, 1422-1423