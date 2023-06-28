# MSA_filtering

This program allows to filter a nucleotide multiple sequences alignment (MSA) based on GUIDANCE (http://guidance.tau.ac.il/) results of an amino-acid MSA. This nucleotide MSA need to be obtained by transposing the non-filtered amino-acid MSA of GUIDANCE (for example through the help of pal2nal: http://www.bork.embl.de/pal2nal/). 

There are two options: 

(1) filter nucleotide MSAs based on the GUIDANCE scores of the amino-acid MSA

(2) remove gaps in MSAs 

Both options can be combined in one run.

WARNING: this program needs python3 and biopython to be installed (pip install biopython), see requirements.txt file.

To install the program either copy the folder in .zip and extract the python script 'cdsMSAfilter_guidance.py' and then run it using on the command line: 'python cdsMSAfilter_guidance.py -h' to display the information on the options needed to be given to run the code, or install with pip, as such: 'pip install "git+https://github.com/Dr-ShinyRaven-Mr-Fox/MSA_filtering.git" and then run it as such: 'MSA_filtering -h'.

To uninstall the program if you installed it with pip, please run on the command line: 'pip uninstall MSA_filtering'

This program has been tested on ubuntu 22.04 and might not work for other OS.
