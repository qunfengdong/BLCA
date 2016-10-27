Bayesian-based LCA taxonomic classification method
--------------------------------------------------
Bayesian LCA-based Taxonomic Classification Method (BLCA) is a Bayesian-based method that provides a solid probabilistic basis for evaluating the taxonomic assignments for the query sequences with bootstrap confidence scores, which is based on Bayesian posterior probability that quantitatively weigh each database hit sequence according to its similarity to the query sequence - the more similar database hit sequence to the query, the more its contribution to the taxonomic assignment of the query. 

We implemented the above algorithm as both a python3 package at https://github.com/qunfengdong/BLCA and a simple python script at https://github.com/qunfengdong/BLCA_commandline.

Depends on your preference, you could either download the python3 package, or the python scripts at current page.

**The following instructions are only for the command line version of BLCA.**

## Prerequisities
* Python 2.7
* Linux
* Biopython
* BLAST binary (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.3.0/)
* MUSCLE (http://www.drive5.com/muscle/downloads.htm)

## Install
Checkout the source code: https://github.com/qunfengdong/BLCA_commandline. To obtain the scripts and example fasta files, do the following:

```shell
$ git clone https://github.com/qunfengdong/BLCA_commandline.git

```

After the github repository is cloned, you will find a folder named BLCA_commandline. All the scripts and example data files will be included in it.

## Quick start

We do not include a pre-compiled database with this release, so the first step is to build a taxonomy database from the NCBI 16s library. And we achieve this by using script _1.subset_db_tax.py_ After the database is built and store on your local machine, you will supply the loction of


## Getting started



## Output
* A text file with sequence id in the first column, taxonomy annotation with confidence scores in the parentheses.
* Each taxonomy level is separated by a semicolon.
* If no taxonomy annotation is available, it is listed as 'NA'


### Example output file:
```
sequence_id  superkingdom	phylum	class	order	family	genus	species subspecies  organism
seq1  Bacteria (100)  Firmicutes (100)  Bacilli (100) Bacillales (100)  Bacillaceae (100) Bacillus (100)  Bacillus subtilis (60)  Bacillus subtilis subsp. spizizenii (43)  Bacillus subtilis subsp. spizizenii (24) 
seq2  Bacteria (100)  Proteobacteria (100)  Gammaproteobacteria (100) Pseudomonadales (100) Pseudomonadaceae (100)  Pseudomonas (100) Pseudomonas mucidolens (21) Pseudomonas cedrina subsp. fulgida (6)  Pseudomonas mucidolens (21)
```

### Description of each column
If you put them in Excel, they should be in the following order.
```
Column A: sequence_id
Column B: superkingdom
Column C: phylum
Column D: class
Column E: order
Column F: family
Column G: genus
Column H: species
Column I: subspecies
Column J: organism
```


## Version
* Version 0.6
first public release

## Authors
* Dr. Xiang Gao, theoretical conception and algorithm development
* Dr. Qunfeng Dong, algorithm development
* Kashi Revanna, program coding and package development
* Huaiying Lin, program coding and testing

## License
GNU

## Acknowledgements
* BLAST program: Camacho C., Coulouris G., Avagyan V., Ma N., Papadopoulos J., Bealer K., & Madden T.L. (2008) "BLAST+: architecture and applications." BMC Bioinformatics 10:421.
* MUSCLE: Edgar, R.C. (2004) MUSCLE: multiple sequence alignment with high accuracy and high throughput.Nucleic Acids Res. 32(5):1792-1797. doi:10.1093/nar/gkh340
* Biopython: Cock PA, Antao T, Chang JT, Bradman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B and de Hoon MJL (2009) Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25, 1422-1423


