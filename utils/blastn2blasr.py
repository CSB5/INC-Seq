#!/usr/bin/env python

"""Convert blastn output to blasr m5 format.

"""

import os
import sys
import argparse


def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile',
                        required = "True",
                        dest="infile",
                        help="Blastn output generated from filter_best_match.py")
    parser.add_argument('-o', '--outfile',
                        help="blasr m5 file",
                        dest='outFile',
                        default=sys.stdout,
                        type=argparse.FileType('w'))

    args = parser.parse_args(arguments)

    if args.infile == '-':
        h = sys.stdin
    else:
        h = open(args.infile, 'rU')

    blastOutFMT = 'sseqid sstart send slen qstart qend qlen evalue score length nident mismatch gaps sseq qseq qseqid'.split()
    blasrFMT = 'qName qLength qStart qEnd qStrand tName tLength tStart tEnd tStrand score numMatch numMismatch numIns numDel mapQV qAlignedSeq matchPattern tAlignedSeq'
    
    for line in h:
        fields = line.strip().split()
        record = dict(zip(blastOutFMT, fields))

        output = [record['sseqid'], record['slen'], record['sstart'], record['send'], ]
        output += [ '+' if int(record['sstart']) < int(record['send']) else '-']
        output += [record['qseqid'], record['qlen'], record['qstart'], record['qend'], '+']
        output += ['-3000']  ## a fake score
        output += [record['nident'], record['mismatch']]
        output += [str(record['qseq'].count('-'))]
        output += [str(record['sseq'].count('-'))]
        output += ['254'] ## fake mapQV
        output += [record['sseq']]
        aln = ''
        for i,j in zip(record['qseq'],record['sseq']):
            aln += '|' if i==j else '*'
        output += [aln]
        output += [record['qseq']]
        
        args.outFile.write(' '.join(output)+'\n')
    h.close()
    
if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
