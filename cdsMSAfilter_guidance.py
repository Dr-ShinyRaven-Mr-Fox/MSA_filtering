#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 10:50:35 2023

@author: Dr-ShinyRaven-Mr-Fox (aka Bertrand Fouks)
"""

import os
from Bio import AlignIO
import argparse 

parser = argparse.ArgumentParser(prog="cdsMSAfilter_guidance",description="This program allows to filter a nucleotide alignment obtained by transposing protein MSAs and performed using GUIDANCE. There are two options: filter nucleotide MSAs using the GUIDANCE results and remove gaps in MSAs. WARNING: this program needs python3 and biopython to be installed (pip install biopython)",epilog="This program is freely available, please cite Fouks et al. 2023 iScience")
parser.add_argument("-PM", "--ProjectMode", choices =['filter+gap', 'filter', 'gap'], required=True, help="Three options are available: filter cds MSAs using guidance scores and remove gaps (filter+gap), filter cds MSAs using guidance scores only (filter), and only remove gaps (gap)") # mode: filter+gap, gap, filter
parser.add_argument("-GF", "--GuidanceFolder", required=True, help = "Please indicate the full path of the proteins MSA GUIDANCE folder")
parser.add_argument("-FS", "--FilteringScore", required=False, default = 0.93, type = float, help = "Please indicate the threshold you want from the GUIDANCE scores. Default is 0.93 as in GUIDANCE default setup")
parser.add_argument("-M", "--MSA", required=True, help = "Please indicate the full path and name of the nucleotide (cds) MSA obtained (e.g. Pal2nal) by transposing amino acid MSA into nucleotide MSA")
parser.add_argument("-FM", "--FilteredMSA", required=True, help = "Please indicate the full path and name of the output MSA")
args = parser.parse_args(args)
    
def rmgap(aln):
    lenaln = aln.get_alignment_length()
    newaln = aln[:,:1]
    i = 0
    while i < lenaln:
        if '-' in aln[:,i]:
            i += 3
        else:
            newaln = newaln + aln[:,i:i+3]
            i += 3
    gdaln = newaln[:,1:]
    return gdaln

def filtaln(caln,path,FilteringScore):
    rmvposS = []
    rmvposE = []
    for file in os.listdir(path):
        if file.endswith('csv'):
            with open('%s/%s' % (path,file)) as files:
                for line in files:
                    if line.startswith('#'):
                        continue
                    col = line.split(',')[0]
                    scor = float(line.strip().split(',')[-1])
                    if scor < FilteringScore:
                        remcol = int(col)
                        strtpos = (remcol*3)-3
                        endpos = (remcol*3)
                        rmvposS.append(strtpos)
                        rmvposE.append(endpos)
            alnF = AlignIO.read('%s' % caln,'fasta')
            ftaln = alnF[:,0:1]
            beg = 0
            i = 0
            for pos in rmvposS:
                if rmvposS[i] == 0:
                    beg = rmvposE[i]
                    i += 1
                    continue
                else:
                    ftaln = ftaln + alnF[:,beg:rmvposS[i]]
                    beg = rmvposE[i]
                    i += 1
                    if i == len(rmvposE):
                        ftaln = ftaln + alnF[:,beg:]
            gftaln = ftaln[:,1:]
    return gftaln

def cdsMSAfilter_guidance(ProjectMode,GuidanceFolder,FilteringScore,MSA,FilteredMSA):        
    if ProjectMode == 'filter+gap':
        fltaln = filtaln(args.MSA,args.GuidanceFolder,FilteringScore)
        gfltaln = rmgap(fltaln)
    elif ProjectMode == 'filter':
        gfltaln = filtaln(MSA,GuidanceFolder,FilteringScore)
    elif ProjectMode == 'gap':
        aln = AlignIO.read('%s' % MSA,'fasta')
        gfltaln = rmgap(aln)
    else:
        print("The mode specified '%s' does not exist" % ProjectMode)
        raise SystemExit(1)
    AlignIO.write(gfltaln,FilteredMSA,'fasta')  

if __name__ == "__main__":    
    import argparse 
    parser = argparse.ArgumentParser(prog="cdsMSAfilter_guidance",description="This program allows to filter a nucleotide alignment obtained by transposing protein MSAs and performed using GUIDANCE. There are two options: filter nucleotide MSAs using the GUIDANCE results and remove gaps in MSAs. WARNING: this program needs python3 and biopython to be installed (pip install biopython)",epilog="This program is freely available, please cite Fouks et al. 2023 iScience")
    parser.add_argument("-PM", "--ProjectMode", choices =['filter+gap', 'filter', 'gap'], required=True, help="Three options are available: filter cds MSAs using guidance scores and remove gaps (filter+gap), filter cds MSAs using guidance scores only (filter), and only remove gaps (gap)") # mode: filter+gap, gap, filter
    parser.add_argument("-GF", "--GuidanceFolder", required=True, help = "Please indicate the full path of the proteins MSA GUIDANCE folder")
    parser.add_argument("-FS", "--FilteringScore", required=False, default = 0.93, type = float, help = "Please indicate the threshold you want from the GUIDANCE scores. Default is 0.93 as in GUIDANCE default setup")
    parser.add_argument("-M", "--MSA", required=True, help = "Please indicate the full path and name of the nucleotide (cds) MSA obtained (e.g. Pal2nal) by transposing amino acid MSA into nucleotide MSA")
    parser.add_argument("-FM", "--FilteredMSA", required=True, help = "Please indicate the full path and name of the output MSA")
    args = parser.parse_args(args)
    cdsMSAfilter_guidance(args.ProjectMode,args.GuidanceFolder,args.FilteringScore,args.MSA,args.FilteredMSA)
    
