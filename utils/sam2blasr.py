#!/usr/bin/env python

"""Convert the sam file from graphmap to blasr m5

"""

import os
import sys
import argparse
from re import split as rs
from Bio import SeqIO
from collections import deque

def split_cigar(cigar):
    ## split the cigar into a list of tuples
    cigar_splitted = rs('(S|I|D|M|X|=)', cigar)
    return zip(cigar_splitted[0::2], cigar_splitted[1::2])

def build_aln(ref, seq, cigar, tstart):
    ## construct the alignment
    qstart = 1
    qend = len(seq)
    if cigar[0][1] == 'S':
        qstart = int(cigar[0][0]) + 1 ## this is one based
        tstart = tstart - int(cigar[0][0])
    if cigar[-1][1] == 'S':
        qend = qend - int(cigar[-1][0])

    tmp_q = deque(seq[qstart-1:])
    tmp_t = deque(ref[tstart + qstart-1:])
    qseq = []
    tseq = []
    mathc = []
    for c in cigar:
        length = int(c[0])
        if c[1] in ['M', 'X', '=']:
            qseq.extend([tmp_q.popleft() for _i in xrange(length)])
            tseq.extend([tmp_t.popleft() for _i in xrange(length)])

        elif c[1] == 'I':
            qseq.extend([tmp_q.popleft() for _i in xrange(length)])
            tseq.extend(['-'] * length)
        elif c[1] == 'D':
            tseq.extend([tmp_t.popleft() for _i in xrange(length)])
            qseq.extend(['-'] * length)
    match = []
    for i,j in zip(qseq,tseq):
        match.append('|' if i==j else '*')

    qseq = ''.join(qseq)
    tseq = ''.join(tseq)
    match = ''.join(match)
    return [str(qstart), str(qend), qseq, match, tseq]

def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--in_sam",
                        required="True",
                        dest="sam",
                        help="The input sam file.")
    parser.add_argument("-r", "--reference",
                        required="True",
                        dest="ref",
                        help="Reference sequence(s).")
    parser.add_argument("-e", "--evalue",
                        default=None,
                        type=float,
                        dest="e_cutoff",
                        help="E-value cutoff.")
    parser.add_argument("--debug",
                        action = "store_true",
                        dest="debug",
                        help="Only print the alignments")
    parser.add_argument('-o', '--outfile',
                        help="Output file",
                        dest='outFile',
                        default=sys.stdout, 
                        type=argparse.FileType('w'))

    args = parser.parse_args(arguments)
    # read the reference sequence
    ref = {}
    refH = open(args.ref, "rU")
    for record in SeqIO.parse(refH, "fasta"):
        ref[record.id] = list(record.seq)
    refH.close()

    if args.sam == '-':
        sam = sys.stdin
    else:
        sam = open(args.sam, 'rU')

    counter = 0
    for l in sam:
        ## skip the headers
        if l[0] != '@':
            fields = l.strip().split("\t")
            #            if counter % 100 == 0:
            #                sys.stderr.write('=')
            # only for graphmap output sam files
            if args.e_cutoff != None:
                e_pass = False
                ZE = fields[-4].split(':')
                if ZE[0]!='ZE':
                    sys.exit("Wrong sam specification (ZE)!")
                evalue = float(ZE[-1])
                if evalue < args.e_cutoff:
                    e_pass = True
            else:
                e_pass = True
                
            if fields[2] != '*' and e_pass:
                counter += 1
                cigar = fields[5]
                read_seq = fields[9]
                tStart = fields[3]
                tName = fields[2]
                qStart, qEnd, qseq, match, tseq =  build_aln(ref[tName], list(read_seq), split_cigar(cigar), int(tStart)-1)
                if args.debug:
                    print qseq
                    print match
                    print tseq
                else:
                    qName = fields[0]
                    qLength = str(len(read_seq))
                    tLength = str(len(tseq))
                    tEnd = str(int(tStart) + len(tseq) - tseq.count('-') - 1)
                    score = '-3000'
                    numMatch = str(match.count('|'))
                    numIns = str(qseq.count('-'))
                    numDel = str(tseq.count('-'))
                    numMismatch = str(match.count('*') - int(numIns) - int(numDel))
                    mapQV = '254'
                    output = ' '.join([qName, qLength, qStart, qEnd, '+', tName, tLength, tStart, tEnd, '+', score, numMatch, numMismatch, numIns, numDel, mapQV, qseq, match, tseq])

                    args.outFile.write(output + '\n')
    ## qName qLength qStart qEnd qStrand tName tLength tStart tEnd tStrand score numMatch numMismatch numIns numDel mapQV qAlignedSeq matchPattern tAlignedSeq ##
    ##sys.stderr.write('\nNumber of aligned reads: ' + str(counter) + '\n')

    
if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
