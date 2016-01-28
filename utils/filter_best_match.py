#!/usr/bin/env python

"""Filter the blast output:
This script takes in a blastn results from find_unit.py, 
scans through the alignment, if there are overlapped ones, 
keeps the one with lower e value. If e value is the same, 
break ties using identity.

Critical column arrangement of the input:
-outfmt "6 sseqid sstart send slen qstart qend qlen ..."

"""

import os
import sys
import argparse

# MIN_OVER_LAP = 1 ## minimal overlap length to be considered as the same location (in this application, the results should not be overlapping)

def getOverlap(a1, a2, b1, b2):
    (a_min, a_max) = sorted([a1,a2])
    (b_min, b_max) = sorted([b1,b2])
    return max(0, min(a_max, b_max) - max(a_min, b_min))

def get_cor(aln):
    cors = aln.strip().split()[0:3]
    return (cors[0], min(int(cors[1]), int(cors[2])))
    

def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('-i', '--infile',
                        required = "True",
                        dest="infile",
                        help="Input file sorted based on positions of reads")
    parser.add_argument('-o', '--outfile', help="Output file with collapsed alignments [Default: stdout]", dest='outfile',
                        default=sys.stdout, type=argparse.FileType('w'))
    parser.add_argument("-c", "--coverage",
                        default=0.9,
                        type = float,
                        dest="coverage",
                        help="The query coverage threshold [default: 0.9]")
    parser.add_argument("-m", "--min_overlap",
                        default=1,
                        type = int,
                        dest="min_overlap",
                        help="Minimal overlap length to be considered for collapse [default: 1]")


    args = parser.parse_args(arguments)
    
    print_flag = 0 ## 0: do not print ;1: safe to print the previous record;
    update_flag = 0 ## 0: do not update; 1: update
    if args.infile == '-':
        inFile = sys.stdin
    else:
        inFile = open(args.infile, 'rU')
        
        ##    previous_fields = 'sseqid sstart send slen qstart qend qlen evalue bitscore length pident mismatch gaps gapopen'.split(' ')

    alns = sorted(inFile.readlines(), key=get_cor)

    previous_fields = alns[0].strip().split()
    if previous_fields != []:
        for record in alns[1:]:
            fields = record.strip().split()

            if (fields[0] == previous_fields[0]) and getOverlap(int(fields[1]),int(fields[2]),int(previous_fields[1]),int(previous_fields[2])) >= args.min_overlap:
                ## The two records overlap, print the one with lower e value
                if float(fields[7])  < float(previous_fields[7]):
                    print_flag = 0 
                    update_flag = 1
                elif float(fields[7])  > float(previous_fields[7]):
                    print_flag = 0
                    update_flag = 0
                else:
                    ## The two e-values equal
                    ## compare the identity
                    if int(fields[9]) > int(previous_fields[9]):
                        print_flag = 0
                        update_flag = 1
                    elif int(fields[9]) <= int(previous_fields[9]):
                        print_flag = 0
                        update_flag = 0
            else:
                ## The two records do not overlap (different reads or non-overlapping regions in the same reads)
                print_flag = 1
                update_flag = 1
            if (print_flag == 1):
                if abs(int(previous_fields[4])-int(previous_fields[5]))*1.0/int(previous_fields[6]) >= args.coverage:
                    args.outfile.write('\t'.join(previous_fields) + "\n")
            if (update_flag == 1):
                previous_fields = fields
            ## print the last record
        if abs(int(previous_fields[4])-int(previous_fields[5]))*1.0/int(previous_fields[6]) >= args.coverage:
            args.outfile.write('\t'.join(previous_fields) + "\n")
if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
