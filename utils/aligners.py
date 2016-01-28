#!/usr/bin/env python

import os
import sys
import subprocess


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

