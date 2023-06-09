# beastXMLgenerator

Scripts for generating the BEAST XML input files

## XML generator v1.1

Generating XML file for BEAST 1.10.4 
 
usage: 

<code>beast1XMLgenerator [-h] [-o OUTPUT] [-t TEMPLATE] [-f FASTA] [-g GFF] [-p PRIORS] [-m MCMC] [-l LOG] [-i ITERATIVE]</code>


optional arguments: 

<code>-h, --help                             [show this help message and exit.]
-o OUTPUT, --output OUTPUT             [The output XML file.]
-t TEMPLATE, --template TEMPLATE       [The XML template file.]
-f FASTA, --fasta FASTA                [The aligned fasta file.]
-g GFF, --gff GFF                      [The gff annotation file.]
-p PRIORS, --priors PRIORS             [The priors table in csv format.]
-m MCMC, --mcmc MCMC                   [The MCMC chain length.]
-l LOG, --log LOG                      [Write the log file.]
-i ITERATIVE, --iterative ITERATIVE    [If you want to generate single-dating file for each undated sample. Potential option are True or False.]</code> 


## XML generator v2.0 (ON TEST)
This script is still on test. 
 
Available options are:
- The input fasta (multiple sequence alignment);
- The annotation file; 
- The tip priors table;
- Partitions (what partitions you want to include or exclude)
- The taxon set;
- Population model;
- Subsitution model; 
- Tree root;
- MCMC chain length and log length;


