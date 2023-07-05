#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 10:50:35 2023

@author: Dr-ShinyRaven-Mr-Fox (aka Bertrand Fouks)
"""

import os
from Bio import AlignIO
import argparse


def rmgap(aln):  # add type hinting and docstring
    """
    add description
    :param aln: add description ...
    :return:
    """
    lenaln = aln.get_alignment_length()
    newaln = aln[:, :1]
    i = 0
    while i < lenaln:
        if '-' in aln[:, i]:
            i += 3
        else:
            newaln = newaln + aln[:, i:i + 3]  # wouldn't "newaln += aln[...]" be better here, but still correct?
            i += 3
    gdaln = newaln[:, 1:]
    return gdaln


def filtaln(caln, path, filtering_score):  # you should also here add type hinting and docstring
    """
    add description
    :param caln: add description ...
    :param path:
    :param filtering_score:
    :return:
    """
    rmvpos_s = []
    rmvpos_e = []
    for file in os.listdir(path):
        if file.endswith('csv'):
            # as stated below I would handle them as pathlib.Path() objects, not strings to avoid tons of problems
            # that can occur here. Use in any case "open(os.path.join(path, file))..." and don't merge file and path as
            # strings. That will likely cause massive problems on different systems or if user specifies path
            # with/without trailing "/"
            with open('%s/%s' % (path, file)) as files:  # path is not a good variable name, because name space clash with os.path
                for line in files:
                    if line.startswith('#'):
                        continue
                    # you could do the line.strip().split() just once and then directly assign int(x[0]) to remcol and
                    # scor accordingly. Now you do two times the line split and have an extra line for remcol +additional unused var col
                    col = line.split(',')[0]
                    scor = float(line.strip().split(',')[-1])
                    if scor < filtering_score:
                        remcol = int(col)
                        strtpos = (remcol * 3) - 3
                        endpos = (remcol * 3)
                        rmvpos_s.append(strtpos)
                        rmvpos_e.append(endpos)
            aln_f = AlignIO.read('%s' % caln, 'fasta')  # why not using directly caln?
            ftaln = aln_f[:, 0:1]
            beg = 0
            i = 0
            for _ in rmvpos_s:  # variable pos was never used, so if it's unnecessary use just underscore as placeholder
                if rmvpos_s[i] == 0:
                    beg = rmvpos_e[i]
                    i += 1
                    continue
                else:
                    ftaln = ftaln + aln_f[:, beg:rmvpos_s[i]]
                    beg = rmvpos_e[i]
                    i += 1
                    if i == len(rmvpos_e):
                        ftaln = ftaln + aln_f[:, beg:]
            gftaln = ftaln[:, 1:]
    return gftaln  # if no file ends with .csv gftaln will be tried to return before it was assigned and throw an error/break your program


# function names should only be lower case letters & argument names start with lower case letter
# & added type hinting and docstrings & all str type containing paths I would consider to turn into Path() objects
# => from pathlib import Path   and then e.g.  my_path = Path(guidance_folder)
def cds_msa_filter(project_mode: str, guidance_folder: str, filtering_score: float, msa: str, filtered_msa: str) \
        -> None:
    """
    add description what the function does
    :param project_mode: add param description ...
    :param guidance_folder:
    :param filtering_score:
    :param msa:
    :param filtered_msa:
    :return: None
    """
    if project_mode == 'filter+gap':
        fltaln = filtaln(msa, guidance_folder, filtering_score)
        gfltaln = rmgap(fltaln)
    elif project_mode == 'filter':
        gfltaln = filtaln(msa, guidance_folder, filtering_score)
    elif project_mode == 'gap':
        aln = AlignIO.read('%s' % msa, 'fasta')  # why not directly just msa variable?
        gfltaln = rmgap(aln)
    else:
        print(f"The specified mode {project_mode} does not exist!")  # I would work with f-string formatting everywhere
        raise SystemExit(1)  # I would consider here rather raising an Error (you can even specify your own "ModeError" or whatever, that will directly print the message above without print statement and automatically exit with corresponding exit code
    AlignIO.write(gfltaln, filtered_msa, 'fasta')


def main() -> None:
    """
    The main function containing argument parsing and running cdsMSAfilter with given arguments
    :return: None
    """
    parser = argparse.ArgumentParser(prog="cdsMSAfilter_guidance",
                                     description="This program allows to filter a nucleotide alignment obtained by "
                                                 "transposing protein MSAs and performed using GUIDANCE. There are two "
                                                 "options: filter nucleotide MSAs using the GUIDANCE results and remove"
                                                 " gaps in MSAs. WARNING: this program needs python3 and biopython to "
                                                 "be installed (pip install biopython)",
                                     epilog="This program is freely available, please cite Fouks et al. 2023 iScience")
    parser.add_argument("-PM", "--ProjectMode", choices=['filter+gap', 'filter', 'gap'], required=True,
                        help="Three options are available: filter cds MSAs using guidance scores and remove gaps "
                             "(filter+gap), filter cds MSAs using guidance scores only (filter), and only remove gaps "
                             "(gap)")  # mode: filter+gap, gap, filter
    parser.add_argument("-GF", "--GuidanceFolder", required=True,
                        help="Please indicate the full path of the proteins MSA GUIDANCE folder")
    parser.add_argument("-FS", "--FilteringScore", required=False, default=0.93, type=float,
                        help="Please indicate the threshold you want from the GUIDANCE scores. Default is 0.93 as in "
                             "GUIDANCE default setup")
    parser.add_argument("-M", "--MSA", required=True,
                        help="Please indicate the full path and name of the nucleotide (cds) MSA obtained (e.g. "
                             "Pal2nal) by transposing amino acid MSA into nucleotide MSA")
    parser.add_argument("-FM", "--FilteredMSA", required=True,
                        help="Please indicate the full path and name of the output MSA")
    args = parser.parse_args()

    cds_msa_filter(args.ProjectMode, args.GuidanceFolder, args.FilteringScore, args.MSA, args.FilteredMSA)


if __name__ == "__main__":
    main()
