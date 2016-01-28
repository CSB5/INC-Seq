#!/usr/bin/env python

"""Convert the sam file from graphmap to blasr m5

"""

import os
import sys
import argparse
import pysam
from Bio import SeqIO

def build_aln(r, refsq):
    rseq = []
    qseq = []
    qstart = r.query_alignment_start
    qend = r.query_alignment_end
    for (qapos, rpos) in r.get_aligned_pairs():
        ## qapos is the aligned index, i.e. this ignores clipping. add that
        if qapos is None:
            qpos = None
        else:
            qpos = qapos + qstart
        ## qpos and rpos now safe to use, but might be None (indel)
        rbase = qbase = "-"
        if rpos is not None:
            rbase = refsq[rpos].upper()
        if qpos is not None:
            qbase = r.seq[qpos].upper()
        rseq.append(rbase)
        qseq.append(qbase)
    match = []
    for i,j in zip(qseq,rseq):
        match.append('|' if i==j else '*')
    return (str(qstart+1), str(qend), ''.join(qseq), ''.join(match) , ''.join(rseq))


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
    refH = open(args.ref, "rU")
    ref = SeqIO.to_dict(SeqIO.parse(refH, "fasta"))
    refH.close()

    if args.sam == '-':
        sam = sys.stdin
    else:
        sam = args.sam

    samfile = pysam.AlignmentFile(sam, "r")

    counter = 0
    for r in samfile.fetch():
        if counter % 100 == 0:
            sys.stderr.write('=')
        if not r.is_unmapped:
            counter += 1
            tName = samfile.getrname(r.reference_id)
            refseq = str(ref[tName].seq)
            qStart, qEnd, qseq, match, tseq =  build_aln(r, refseq)
            if args.debug:
                print qseq
                print match
                print tseq
            else:
                qName = r.query_name
                qLength = str(r.query_length)
                tLength = str(len(refseq))
                tStart = str(r.reference_start + 1)
                tEnd = str(r.reference_end)
                score = '-3000'
                numMatch = str(match.count('|'))
                numIns = str(qseq.count('-'))
                numDel = str(tseq.count('-'))
                numMismatch = str(match.count('*') - int(numIns) - int(numDel))
                mapQV = '254'
                output = ' '.join([qName, qLength, qStart, qEnd, '+', tName, tLength, tStart, tEnd, '+', score, numMatch, numMismatch, numIns, numDel, mapQV, qseq, match, tseq])

                args.outFile.write(output + '\n')
    ## qName qLength qStart qEnd qStrand tName tLength tStart tEnd tStrand score numMatch numMismatch numIns numDel mapQV qAlignedSeq matchPattern tAlignedSeq ##
    sys.stderr.write('\nNumber of aligned reads: ' + str(counter) + '\n')

    
if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
