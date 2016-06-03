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
        filteredblastout.stdout.close()

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


def marginAlign(query, ref, tmpName):
    script_dir = os.path.dirname(os.path.realpath(__file__))
    marginAlign_dir = "/mnt/projects/lich/dream_challenge/rollingcircle/finalized/ana_scripts/marginAlign"
    ## run
    cmd = [marginAlign_dir+"/marginAlign", query, ref, tmpName+".sam", "--jobTree", tmpName+".jobTree", "--em" ]
    with open(tmpName+".margin.log", 'w') as log:
        tmp = subprocess.check_output(" ".join(cmd), shell = True, stderr=log)
    cmd = [script_dir+"/sam2blasr.py", "-i", tmpName+".sam", "-r", ref]
    stdout = subprocess.check_output(" ".join(cmd), shell = True)
    tmp = subprocess.check_output("rm -rf %s %s" %(tmpName+".jobTree", tmpName+".sam"), shell = True)
    return stdout



def poa(fasta, tmpName, seqHeader):
    script_dir = os.path.dirname(os.path.realpath(__file__))
    cmd = [script_dir+"/poa", "-do_global", "-do_progressive",
           "-read_fasta", fasta,
           "-pir",  "pseudo",
           script_dir+"/blosum80.mat", "-hb"]
    with open(tmpName+".poa.log", 'w') as log:
        poaout = subprocess.check_output(" ".join(cmd), stderr=log, shell = True)
    
    ## post processing
    consensus = poaout.split(">")[-1]
    consensus = ''.join(consensus.split("\n")[1:])
    consensus = consensus.replace(".","")
    consensus = "\n".join([seqHeader, consensus]) + "\n"
    
    return consensus
