#!/usr/bin/env python

"""Convert blastn files to bed format

"""

import os
import sys
import argparse


def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile',
                        required = "True",
                        dest="infile",
                        help="Blastn output generated from customized format '6 sseqid sstart send ...'")
    parser.add_argument('-o', '--outfile',
                        help="bed file name",
                        dest='outFile',
                        default=sys.stdout,
                        type=argparse.FileType('w'))

    args = parser.parse_args(arguments)

    with open(args.infile) as f:
        for line in f:
            fields = line.strip().split()[0:3]
            name = fields[0]
            pos = [int (i) for i in fields[1:]]
            start = str(min(pos) - 1)
            end = str(max(pos))
            args.outFile.write('\t'.join([name, start, end])+'\n')


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
