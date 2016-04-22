#!/usr/bin/env python

"""The INC-Seq pipeline
"""
import os
import sys
import argparse
import subprocess

from datetime import datetime
from Bio import SeqIO
from utils import findUnit, buildConsensus

def get_tmp(program):
    program_time = datetime.now().strftime('%Y-%m-%d_%H-%M-%S.%f')
    tmp_folder = "/dev/shm/" if os.path.exists("/dev/shm") else "/tmp/"
    tmp_folder += program + program_time + "/"
    return tmp_folder

def callBuildConsensus(aligner, record, aln, copy_num_thre, len_diff_thre, tmp_folder, seg_cov, iterative):
    if aligner == "blastn":
        consensus = buildConsensus.consensus_blastn(record, aln, copy_num_thre,
                                                    len_diff_thre, tmp_folder,
                                                    seg_cov, iterative)
    elif aligner == "graphmap":
        consensus = buildConsensus.consensus_graphmap(record, aln, copy_num_thre,
                                                      len_diff_thre, tmp_folder,
                                                      seg_cov, iterative)
    elif aligner == "poa":
        consensus = buildConsensus.consensus_poa(record, aln, copy_num_thre,
                                                 len_diff_thre, tmp_folder)
    elif aligner == "marginAlign":
        consensus = buildConsensus.consensus_marginAlign(record, aln, copy_num_thre,
                                                         len_diff_thre, tmp_folder,
                                                         seg_cov, iterative)
    return consensus

def main(arguments):
    #### parsing arguments
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', "--input",
                        required="True",
                        help="Input file in fasta format",
                        dest="inFasta")
    parser.add_argument('-o', '--outfile',
                        help="Output file",
                        dest = "outFile",
                        default=sys.stdout, type=argparse.FileType('w'))
    parser.add_argument("-a", "--aligner",
                        default='blastn',
                        dest="aligner",
                        help="The aligner used (blastn, graphmap, poa) [Default: blastn]")
    parser.add_argument("-m", "--minReadLength",
                        default=2000,
                        dest="minRL",
                        type=int,
                        help="The reads shorter than this will be discarded [Default:2000]")
    ##find unit specific
    parser.add_argument("--anchor_seg_step",
                        default=500,
                        dest="anchor_seg_step",
                        type=int,
                        help="Step of sliding window used as anchors [Default: 500] (eg. -s 500 : start at 0, 500, 1000, ...)")
    parser.add_argument("--anchor_length",
                        default=500,
                        dest="anchor_len",
                        type=int,
                        help="The length of the anchor, should be smaller than the unit length [Default: 500]")
    parser.add_argument("--anchor_cov",
                        default=0.8,
                        dest="anchor_cov",
                        type=float,
                        help="Anchor coverage required [Default: 0.8]")
    parser.add_argument("--anchor_seq",
                        dest="anchor_seq",
                        type=str,
                        help="A single file containing the sequences used as the anchor [Default: Use subsequences as anchors]")
    ##consensus building specific
    parser.add_argument("--iterative",
                        action = "store_true",
                        dest="iterative",
                        help="Iteratively run pbdagcon on consensus [Default: False]")
    parser.add_argument("--seg_cov",
                        default=0.8,
                        dest="seg_cov",
                        type=float,
                        help="Segment coverage required [Default: 0.8]")
    parser.add_argument("--copy_num_thre",
                        dest="copy_num_thre",
                        default = 6,
                        type = int,
                        help="Minimal copy number required [Default: 6]")
    parser.add_argument("--length_difference_threshold",
                        dest="len_diff_thre",
                        default = 0.05,
                        type = float,
                        help="Segment length deviation from the median to be considered as concordant [Default: 0.05]")

    args = parser.parse_args(arguments)
    
    #### parse input reads    
    seqs = SeqIO.parse(args.inFasta, "fasta")

    #### create temp folder
    tmp_folder = get_tmp('incseq_' + args.inFasta.split("/")[-1] + '_')
    os.makedirs(tmp_folder)
    
    counter = 0

    for record in seqs:
        seqlen = len(record.seq)
        sys.stderr.write("---------- Processing read %i ----------\n" % (counter + 1))
        counter += 1
        if seqlen < args.minRL:
            #### length filter
            sys.stderr.write("Failed to pass length filter!\n")
        else:
            #### find units
            if args.aligner == "blastn" or args.aligner == "graphmap" or args.aligner =="poa" or args.aligner == "marginAlign": ## FIXME graphmap implementation
                if args.anchor_seq:
                    ## anchor sequence provided, run with INC-Seq2 mode
                    aln = findUnit.find_unit_blastn(record, args.anchor_seq, tmp_folder, seqlen,
                                                    args.anchor_seg_step,
                                                    args.anchor_len,
                                                    args.anchor_cov)
                else:
                    ## use subsequences as anchors (INC-Seq mode)
                    aln = findUnit.find_unit_blastn(record, None, tmp_folder, seqlen,
                                                    args.anchor_seg_step,
                                                    args.anchor_len,
                                                    args.anchor_cov)

            #### build consensus
            consensus = callBuildConsensus(args.aligner, record, aln, args.copy_num_thre,
                                           args.len_diff_thre, tmp_folder,
                                           args.seg_cov, args.iterative)
 
            if consensus:                    
                    sys.stderr.write("Consensus called\t%s\tNumber of segments\t%d\n" %(record.id, consensus[1]))
                    args.outFile.write(consensus[0])
            else:
                sys.stderr.write("Consensus construction failed!\n")
    os.rmdir(tmp_folder)
    
if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
