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

Usage:
--------------
```
usage: inc-seq.py [-h] -i INFASTA [-o OUTFILE] [-a ALIGNER] [-m MINRL]
                  [--anchor_seg_step ANCHOR_SEG_STEP]
                  [--anchor_length ANCHOR_LEN] [--anchor_cov ANCHOR_COV]
                  [--anchor_seq ANCHOR_SEQ] [--iterative] [--seg_cov SEG_COV]
                  [--copy_num_thre COPY_NUM_THRE]
                  [--length_difference_threshold LEN_DIFF_THRE]

The INC-Seq pipeline

optional arguments:
  -h, --help            show this help message and exit
  -i INFASTA, --input INFASTA
                        Input file in fasta format
  -o OUTFILE, --outfile OUTFILE
                        Output file
  -a ALIGNER, --aligner ALIGNER
                        The aligner used (blastn, graphmap, poa) [Default:
                        blastn]
  -m MINRL, --minReadLength MINRL
                        The reads shorter than this will be discarded
                        [Default:2000]
  --anchor_seg_step ANCHOR_SEG_STEP
                        Step of sliding window used as anchors [Default: 500]
                        (eg. -s 500 : start at 0, 500, 1000, ...)
  --anchor_length ANCHOR_LEN
                        The length of the anchor, should be smaller than the
                        unit length [Default: 500]
  --anchor_cov ANCHOR_COV
                        Anchor coverage required [Default: 0.8]
  --anchor_seq ANCHOR_SEQ
                        A single file containing the sequences used as the
                        anchor [Default: Use subsequences as anchors]
  --iterative           Iteratively run pbdagcon on consensus [Default: False]
  --seg_cov SEG_COV     Segment coverage required [Default: 0.8]
  --copy_num_thre COPY_NUM_THRE
                        Minimal copy number required [Default: 6]
  --length_difference_threshold LEN_DIFF_THRE
                        Segment length deviation from the median to be
                        considered as concordant [Default: 0.05]
```
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
* Use bpipe pipeline for pseudo-parallel computing
 * Split the reads into multiple files (300 reads per file) and run INC-Seq (4 instances) in parallel.
```
bpipe run -p READ_NUM=300 -n 4 pipeline.bpipe a_lot_of_incseq_reads.fa
```
