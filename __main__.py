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

cdsMSAfilter(args.ProjectMode,args.GuidanceFolder,args.FilteringScore,args.MSA,args.FilteredMSA)
