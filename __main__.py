#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 10:50:35 2023

@author: Dr-ShinyRaven-Mr-Fox (aka Bertrand Fouks)
"""

import os
from Bio import AlignIO
import argparse 
from cdsMSAfilter_guidance import cdsMSAfilter

parser = argparse.ArgumentParser(prog="cdsMSAfilter_guidance",description="This program allows to filter a nucleotide alignment obtained by transposing protein MSAs and performed using GUIDANCE. There are two options: filter nucleotide MSAs using the GUIDANCE results and remove gaps in MSAs. WARNING: this program needs python3 and biopython to be installed (pip install biopython)",epilog="This program is freely available, please cite Fouks et al. 2023 iScience")

parser.add_argument("-PM", "--ProjectMode", choices =['filter+gap', 'filter', 'gap'], required=True, help="Three options are available: filter cds MSAs using guidance scores and remove gaps (filter+gap), filter cds MSAs using guidance scores only (filter), and only remove gaps (gap)") # mode: filter+gap, gap, filter

parser.add_argument("-GF", "--GuidanceFolder", required=True, help = "Please indicate the full path of the proteins MSA GUIDANCE folder")

parser.add_argument("-FS", "--FilteringScore", required=False, default = 0.93, type = float, help = "Please indicate the threshold you want from the GUIDANCE scores. Default is 0.93 as in GUIDANCE default setup")

parser.add_argument("-M", "--MSA", required=True, help = "Please indicate the full path and name of the nucleotide (cds) MSA obtained (e.g. Pal2nal) by transposing amino acid MSA into nucleotide MSA")

parser.add_argument("-FM", "--FilteredMSA", required=True, help = "Please indicate the full path and name of the output MSA")

args = parser.parse_args()

cdsMSAfilter(args.ProjectMode,args.GuidanceFolder,args.FilteringScore,args.MSA,args.FilteredMSA)

