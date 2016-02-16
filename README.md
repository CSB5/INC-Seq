INC-Seq: Accurate single molecule reads using nanopore sequencing
======
Description:
------
This repository contains the code for analyzing INC-Seq data (http://biorxiv.org/content/early/2016/01/27/038042). The full datasets have been deposited into ENA (http://www.ebi.ac.uk/ena/data/view/PRJEB12294).

Requirements:
--------------
 - Python 2.7
 - Biopython 1.65
 - BLAST 2.2.28+
 - PBDAGCON (https://github.com/PacificBiosciences/pbdagcon)

Examples:
--------------
* Basic usage
```
./inc-seq.py -i data/inc_seq_test_read.fa -o consensus.fa
```
* Use graphmap as segment aligner
```
./inc-seq.py -i data/inc_seq_test_read.fa -o consensus.fa -a graphmap
```

