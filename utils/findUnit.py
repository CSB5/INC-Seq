#!/usr/bin/env python

import os
import sys
import hashlib

import subprocess
from aligners import *

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

                
################################################################################

###############################find unit########################################
def best_aln(aln):
    ## find the best aln according to e-value
    best_alignment = ''
    best_e = 1000
    for l in aln.split('\n'):
        if l != '':
            evalue = float(l.split()[7])
            if evalue <= best_e:
                best_alignment = l
                best_e = evalue
    return best_alignment

def find_unit_blastn(record, ref_anchor, tmp_folder, seqlen, query_seg_step, query_len, anchor_cov):
    tmpname = tmp_folder + hashlib.md5(record.id).hexdigest() + ".tmp"
    tmpRef = tmpname + ".ref.fasta"
    tmpQ = tmpname + ".q.fasta"
    blastOutFMT = '6 sseqid sstart send slen qstart qend qlen evalue score length nident mismatch gaps'

    alignments = {'alignments':'\n', 'number':0}

    ## write the ref seq (single seq)
    with open(tmpRef, "w") as ref_handle:
        SeqIO.write(record, ref_handle, "fasta")

    if ref_anchor:
        ## ref anchor is provided
        ## firstly map the ref anchor to the read (best mapping)
        ## then extract 100 bps from the ref anchor and use it as the new anchor
        stdout = blastn(ref_anchor, tmpRef, None, blastOutFMT,
                        1, anchor_cov, False)
        best_alignment=best_aln(stdout)
        if best_alignment != '':
            s_start = int((best_aln(stdout)).split()[1])
            with open(tmpQ, 'w') as q_handle:
                qrecord = SeqRecord(record.seq[s_start:s_start+400],
                                    record.id+ "RefAnchor",
                                    description= "")
                SeqIO.write(qrecord, q_handle, "fasta")
            stdout = blastn(tmpQ, tmpRef, None, blastOutFMT,
                            seqlen/query_len + 1, anchor_cov, False)
        alignments['number'] = max(stdout.count('\n') - 1,0)
        alignments['alignments'] = stdout
    else:
        ## try different anchors
        starts = xrange(0, seqlen/2, query_seg_step)
        
        for start in starts:
            ## write the query seq
            with open(tmpQ, 'w') as q_handle:
                qrecord = SeqRecord(record.seq[start:(start+query_len)],
                                    record.id+ str(start) + 'to' + str(query_len) + "bps",
                                    description= "")
                SeqIO.write(qrecord, q_handle, "fasta")
            stdout = blastn(tmpQ, tmpRef, None, blastOutFMT,
                            seqlen/query_len + 1, anchor_cov, False)
            num_alignments = stdout.count('\n') - 1 
            if num_alignments > alignments['number']:
                alignments['number'] = num_alignments
                alignments['alignments'] = stdout

    # finished one read, clean tmp files
    tmp = subprocess.check_output("rm  %s*" % (tmpname), shell = True)
    sys.stderr.write("Max number of segments found: %i \n" % (alignments['number']))
    if alignments['alignments'] != '\n' and alignments['alignments'] != '':
        return alignments['alignments']
    return None

################################################################################

