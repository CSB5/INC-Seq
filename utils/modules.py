#!/usr/bin/env python

import os
import sys
from datetime import datetime
import subprocess
import hashlib
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def get_tmp(program):
    program_time = datetime.now().strftime('%Y-%m-%d_%H-%M-%S.%f')
    tmp_folder = "/dev/shm/" if os.path.exists("/dev/shm") else "/tmp/"
    tmp_folder += program + program_time + "/"
    return tmp_folder

################################aligners########################################
def blastn(query, ref, wordSize, blastOutFMT, num_alignments, anchor_cov, toBlasr):
    script_dir = os.path.dirname(os.path.realpath(__file__))
    ## blastn
    cmd = ["blastn", "-query", query, "-task", "blastn"]# , "-evalue","0.000001"]
    if wordSize:
        cmd.extend(["-word_size", str(wordSize)])
    cmd.extend(["-subject", ref, "-num_alignments", str(num_alignments)])
    cmd.extend(["-outfmt", blastOutFMT])
    blastout = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    #cmd = ["sort","-k","1,1","-k","2,2n"]
    #sortedblastout = subprocess.Popen(cmd, stdin = blastout.stdout, stdout= subprocess.PIPE)
    cmd = [script_dir+"/filter_best_match.py","-i","-", "-c", str(anchor_cov)]
    filteredblastout = subprocess.Popen(cmd, stdin = blastout.stdout, stdout = subprocess.PIPE)
    if not toBlasr:
        stdout, stderr = filteredblastout.communicate()
        blastout.stdout.close()

    else:
        cmd = [script_dir+"/blastn2blasr.py", "-i", "-"]
        blast2blasr = subprocess.Popen(cmd, stdin = filteredblastout.stdout, stdout = subprocess.PIPE)            
        stdout, stderr = blast2blasr.communicate()
        blast2blasr.stdout.close()
        
    return stdout

def graphmap(query, ref):
    script_dir = os.path.dirname(os.path.realpath(__file__))
    ## run graphmap
    cmd = [script_dir+"/graphmap", "-t", "1", "-d", query, "-r", ref, "-v", "0", "-z", "0.000001"]
    graphmapout = subprocess.Popen(cmd, stdout= subprocess.PIPE)
    cmd = [script_dir+"/sam2blasr.py", "-i", "-", "-r", ref]
    sam2blasr = subprocess.Popen(cmd, stdin = graphmapout.stdout, stdout = subprocess.PIPE)
    stdout, stderr = sam2blasr.communicate()
    sam2blasr.stdout.close()

    return stdout
################################################################################

################################find primer location############################
def locate_primer(primer_fwd, primer_rev, consensus, tmp_folder, seqlen):
    tmpname = tmp_folder + hashlib.md5("primer").hexdigest() + ".tmp"
    tmpRef = tmpname + ".ref.fasta"
    tmpQ = tmpname + ".q.fasta"
    blastOutFMT = '6 sseqid sstart send slen qstart qend qlen evalue score length nident mismatch gaps qseqid'
    with open(tmpQ, 'w') as q_handle:
        for primer, name in zip((primer_fwd, primer_rev), ("fwd","rev")):
            qrecord = SeqRecord(Seq(primer),
                                name,
                                description= "")
            SeqIO.write(qrecord, q_handle, "fasta")

    with open(tmpRef, 'w') as q_handle:
        qrecord = SeqRecord(Seq(consensus),
                            "consensus",
                            description= "")
        SeqIO.write(qrecord, q_handle, "fasta")
    
    stdout = blastn(tmpQ, tmpRef, 4, blastOutFMT,
                    1, 0.3, False)
    alns = stdout.strip().split("\n")
    aln = ''
    lowest_e = 1
    for line in alns:
        fields = line.split("\t")
        e_value = float(fields[7])
        if e_value < lowest_e:
            lowest_e = e_value
            aln = line
    tmp = subprocess.check_output("rm  %s*" % (tmpname), shell = True)
    if aln == '':
        ## cannot find primer, discard
        return None
    else:
        fields = aln.split("\t")
        start, end, length, offset = [ int(x) for x in fields[1:5] ]
            
        if start < end:
            ## mapped in fwd orientation
            s_start = max(start - offset, 0)
        else:
            ## mapped in rev orientation
            s_start = min(start + offset, length)
        ## restore the correct position
        return (consensus[s_start:] + consensus[0:s_start])

            
    
################################################################################

###############################find unit########################################
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
        ## ref anchor is provided, use it
        with open(tmpQ, 'w') as q_handle:
            qrecord = SeqRecord(Seq(ref_anchor),
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

################################consensus building##############################
#--------------functions used for consensus building-----------
def median(mylist):
    sorts = sorted(mylist)
    length = len(sorts)
    if length == 0:
        return 0
    if not length % 2:
        return (sorts[length / 2] + sorts[length / 2 - 1]) / 2.0
    return sorts[length / 2]

def get_errors(alignments):
    errors = 0
    count = 0
    for l in alignments.split('\n'):
        if l != '':
            count += 1
            fields = l.split()
            errors += sum([ int(x) for x in fields[12:15] ])
    return errors*1.0/count

def post_processing(alignments):
    ## remove the self-alignments
    out = []
    for l in alignments.split('\n'):
        if l != '':
            fields = l.split()
            if fields[0] != fields[5]:
                out.append(l)
    return '\n'.join(out)

def segment_filter_orientation(aln):
    ## first filter by orientaion
    ## also return the lengths of the segments
    lengths = []
    coordinates = []
    for (i, j) in zip(aln[0:], aln[1:]):
        # check concordance
        if (int(i[0]) - int(i[1])) * (int(j[0]) - int(j[1])) > 0:
            # concordant prepare output
            # gather coordinates on the reads (be more relaxed at two ends)
            if int(i[0]) - int(i[1]) < 0: # in forward directions
                start = max(int(i[0]) - int(i[3]) + 1, 1) # avoid overshoot to negative coordinates
                ##start = int(i[0])
                end = int(j[0])
            else:                         # in reverse complement
                start = int(i[0])
                end = int(j[0]) + int(j[3]) - 1
                ##end = int(j[0])
            lengths.append(end - start + 1)
            coordinates += [(start,end)]
        else:
            # discordant
            coordinates.append('#')
            lengths.append("#")
    return (coordinates, lengths)

def segment_filter_lengths(coordinates, lengths, len_median, len_diff_thre):
    coordinates_filtered = []
    for l, c in zip(lengths, coordinates):
        if l == '#':
            # discordant
            coordinates_filtered.append(c)
        else:
            # concordant
            if abs(l-len_median) > len_median*len_diff_thre:
                # 0.066 for pacbio #expected length discrepency -- 0.22 (error) * (0.6-0.3) (insertion-deletion)
                # wrong length (indication of a chimera)
                coordinates_filtered.append("*")
            else:
                # correct length
                coordinates_filtered.append(c)
    coordinates_filtered.append("#") ## add a delimiter to the end
    return coordinates_filtered

def segment_filter_longest_strech(coordinates_filtered):
    candidate = []
    candidate_cur = []
    for cor in coordinates_filtered:
        if cor == "*" or cor == "#": ## segment boundary:
            if len(candidate_cur) > len(candidate):
                ## found a strech with more segments
                candidate = candidate_cur
            candidate_cur = []
        else:
            candidate_cur.append(cor)                                        
    return candidate

#---------------------------------------------------------------------
def consensus_blastn(record, alnFile, copy_num_thre, len_diff_thre, tmp_folder, seg_cov, iterative):
    if not alnFile:
        return None
    
    aln = []    
    for line in alnFile.strip().split('\n'):
        fields = line.split()
        aln += [ fields[1:] ]

    ## split into segments:
    ## the anchors flanking the segment must be concordant
    ##   1. in the same direction
    ##   2. the length must not be so different
    ##   3. the longest strech of concordant segments will be considered
    candidate_read_count = 0

    if len(aln) >= copy_num_thre:
        seg_coordinates = []
        ## filter for direction
        coordinates, lengths = segment_filter_orientation(aln)
        len_median = median([x for x in lengths if x != "#"])
        ## filter for length
        coordinates_filtered = segment_filter_lengths(coordinates, lengths, len_median, len_diff_thre)
        ## find the longest strech of segments in concordance
        candidate = segment_filter_longest_strech(coordinates_filtered)

        sys.stderr.write("Number of segments of the candidate strech: %d\n" %(len(candidate)))

        ## need some copies for correction
        if len(candidate) >= copy_num_thre:
            seg_coordinates = candidate
            candidate_read_count += 1
            sys.stderr.write("Candidate read found!\n")
        else:
            sys.stderr.write("Not enough alignmets!\n")
            return None
    else:
        sys.stderr.write("Not enough alignmets!\n")
        return None

    #### split into segments and call consensus
    ## split this read into a multiple fasta file
    tmpname = tmp_folder + hashlib.md5(record.id).hexdigest() + ".tmp"
    tmpRef = tmpname + ".ref.fasta"
    tmpQ = tmpname + ".q.fasta"

    seg_num = len(seg_coordinates)

    ## three fields added to facilate convertion to blasr m5 format
    blastOutFMT = '6 sseqid sstart send slen qstart qend qlen evalue score length nident mismatch gaps sseq qseq qseqid'
    ## try using each subread as the backbone
    alignments = {"alignments":'\n',"num":0, "errors":sys.maxint}
    subReads = []

    counter = 0
    # write refs (all the subreads)
    ref_handle = open(tmpRef, 'w')
    for s, e in seg_coordinates:
        counter += 1
        subRead = SeqRecord(record.seq[s-1:e], record.id+'_'+str(counter), description="")
        subReads.append(subRead)
        SeqIO.write(subRead, ref_handle, "fasta")
    ref_handle.close()
    # blast alignment
    for subRead in subReads:
        q_handle = open(tmpQ, 'w')
        SeqIO.write(subRead, q_handle, "fasta")
        q_handle.close()
        stdout = blastn(tmpQ, tmpRef, None, blastOutFMT, seg_num, seg_cov, True)
        num = stdout.count('\n') - 1
        errors = get_errors(stdout)
        if num > alignments["num"]:
            # prefer more alignments
            alignments["num"] = num
            alignments["error"] = errors
            alignments["alignments"] = stdout
        elif num == alignments["num"]:
            # for the same number of alignments, prefer the one with lower error rates
            if errors < alignments["errors"]:
                    alignments["num"] = num
                    alignments["error"] = errors
                    alignments["alignments"] = stdout              
        
    copy_num = alignments["alignments"].count("\n")
    if copy_num >= copy_num_thre:
        with open(tmpname + '.m5', 'w') as outH:
            outH.write(post_processing(alignments["alignments"]))
        consensus = subprocess.check_output(("pbdagcon -t2 -c1 -m1 %s" % (tmpname + '.m5')).split())

        ## run iteratively
        #----------------------------------------
        # write consensus seq 0
        if iterative:
            delta = 1
            consensus_p = consensus
            iteration = 0
            while (delta>0.001 and consensus and iteration<=10):
                iteration += 1
                with open(tmpname + '.con.iter.fa', 'w') as outH:
                    outH.write(consensus)
                stdout = blastn(tmpname+'.con.iter.fa', tmpRef, None, blastOutFMT, seg_num, seg_cov, True)
                copy_num = stdout.count('\n')
                ## write new m5
                with open(tmpname + '.con.iter.m5', 'w') as outH:
                    outH.write(stdout)
                consensus_p = consensus
                consensus = subprocess.check_output(("pbdagcon -t2 -c1 -m1 %s" % (tmpname + '.con.iter.m5')).split())
                with open(tmpname + '.con.iter.n.fa', 'w') as outH:
                    outH.write(consensus)
                # update delta
                tmp = blastn(tmpname+'.con.iter.fa', tmpname+'.con.iter.n.fa', None, blastOutFMT, 1, seg_cov, False)
                tmp_len, tmp_iden = tmp.strip().split()[9:11]
                delta = 1-float(tmp_iden)/int(tmp_len)
                sys.stderr.write("Delta: %f\n" %(delta))
            #----------------------------------------
            consensus = consensus_p
        tmp = subprocess.check_output("rm  %s*" % (tmpname), shell = True)
        return (consensus, copy_num)
    else:
        sys.stderr.write("Not enough aligned copy to correct!\n")
        tmp = subprocess.check_output("rm  %s*" % (tmpname), shell = True)
        return None



def consensus_graphmap(record, alnFile, copy_num_thre, len_diff_thre, tmp_folder, seg_cov, iterative):
    if not alnFile:
        return None
    
    aln = []    
    for line in alnFile.strip().split('\n'):
        fields = line.split()
        aln += [ fields[1:] ]

    ## split into segments:
    ## the anchors flanking the segment must be concordant
    ##   1. in the same direction
    ##   2. the length must not be so different
    ##   3. the longest strech of concordant segments will be considered
    candidate_read_count = 0

    if len(aln) >= copy_num_thre:
        seg_coordinates = []
        ## filter for direction
        coordinates, lengths = segment_filter_orientation(aln)
        len_median = median([x for x in lengths if x != "#"])
        ## filter for length
        coordinates_filtered = segment_filter_lengths(coordinates, lengths, len_median, len_diff_thre)
        ## find the longest strech of segments in concordance
        candidate = segment_filter_longest_strech(coordinates_filtered)

        sys.stderr.write("Number of segments of the candidate strech: %d\n" %(len(candidate)))

        ## need some copies for correction
        if len(candidate) >= copy_num_thre:
            seg_coordinates = candidate
            candidate_read_count += 1
            sys.stderr.write("Candidate read found!\n")
        else:
            sys.stderr.write("Not enough alignmets!\n")
            return None
    else:
        sys.stderr.write("Not enough alignmets!\n")
        return None

    #### split into segments and call consensus
    else: ## run graphmap
        ## split this read into a multiple fasta file
        tmpname = tmp_folder + hashlib.md5(record.id).hexdigest() + ".tmp"
        tmpRef = tmpname + ".ref.fasta"
        tmpQ = tmpname + ".q.fasta"

        ## try using each subread as the backbone, select the backbone with minimal errors
        alignments = {"alignments":'\n',"errors":sys.maxint}
        subReads = []

        counter = 0
        # write queries (all the subreads)
        q_handle = open(tmpQ, 'w')
        for s, e in seg_coordinates:
            counter += 1
            subRead = SeqRecord(record.seq[s-1:e], record.id+'_'+str(counter), description="")
            subReads.append(subRead)
            SeqIO.write(subRead, q_handle, "fasta")
        q_handle.close()
        # graphmap
        for subRead in subReads:
            ref_handle = open(tmpRef, 'w')
            SeqIO.write(subRead, ref_handle, "fasta")
            ref_handle.close()
            stdout = graphmap(tmpQ, tmpRef)
            errors = get_errors(stdout)
            if errors < alignments["errors"]:
                alignments["errors"] = errors
                alignments["alignments"] = stdout
            ## remove index
            tmp = subprocess.check_output("rm  %s*" % (tmpRef), shell = True)
            
    copy_num = alignments["alignments"].count("\n")
    if copy_num >= copy_num_thre:
        with open(tmpname + '.m5', 'w') as outH:
            outH.write(post_processing(alignments["alignments"]))
        consensus = subprocess.check_output(("pbdagcon -t2 -c1 -m1 %s" % (tmpname + '.m5')).split())

        ## run iteratively
        #----------------------------------------
        # write consensus seq 0
        if iterative:
            delta = 1
            consensus_p = consensus
            iteration = 0
            while (delta>0.001 and consensus and iteration<=10):
                iteration += 1
                with open(tmpname + '.con.iter.fa', 'w') as outH:
                    outH.write(consensus)
                stdout = blastn(tmpname+'.con.iter.fa', tmpRef, None, blastOutFMT, seg_num, seg_cov, True)
                copy_num = stdout.count('\n')
                ## write new m5
                with open(tmpname + '.con.iter.m5', 'w') as outH:
                    outH.write(stdout)
                consensus_p = consensus
                consensus = subprocess.check_output(("pbdagcon -t2 -c1 -m1 %s" % (tmpname + '.con.iter.m5')).split())
                with open(tmpname + '.con.iter.n.fa', 'w') as outH:
                    outH.write(consensus)
                # update delta
                tmp = blastn(tmpname+'.con.iter.fa', tmpname+'.con.iter.n.fa', None, blastOutFMT, 1, seg_cov, False)
                tmp_len, tmp_iden = tmp.strip().split()[9:11]
                delta = 1-float(tmp_iden)/int(tmp_len)
                sys.stderr.write("Delta: %f\n" %(delta))
            #----------------------------------------
            consensus = consensus_p
        tmp = subprocess.check_output("rm  %s*" % (tmpname), shell = True)
        return (consensus, copy_num)
    else:
        sys.stderr.write("Not enough aligned copy to correct!\n")
        tmp = subprocess.check_output("rm  %s*" % (tmpname), shell = True)
        return None

################################################################################