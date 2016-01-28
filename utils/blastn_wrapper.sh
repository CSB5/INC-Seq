#!/bin/bash
Q=$1
DB=$2
OUT=$3
NTHREDS=$4
NALN=$5

blastn -perc_identity 90  -query $Q -task blastn -evalue 0.1  -db $DB -outfmt "6 sseqid sstart send slen qstart qend qlen evalue length nident mismatch gaps sseq qseq qseqid" -num_alignments $NALN  -num_threads $NTHREDS > $OUT

