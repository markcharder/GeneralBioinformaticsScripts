#!/usr/bin/env bash

## This script was used to run the REPET TE annotation pipeline.

## Set environment variables.
export REPET_PATH=/home/mark/bioinformatics/src/REPET_linux-x64-2.2/
export PATH=/home/mark/bioinformatics/src/REPET_linux-x64-2.2/bin:$PATH
export PATH=/home/mark/bioinformatics/src/REPET_linux-x64-2.2/denovo_pipe:$PATH
export PATH=/home/mark/bioinformatics/src/gt-1.5.8-Linux_x86_64-64bit-complete/bin:$PATH
export PATH=/home/mark/bioinformatics/src/RECON-1.07/scripts/:$PATH
export PATH=/home/mark/bioinformatics/src/hmmer-3.1b2-linux-intel-x86_64/binaries:$PATH
export PATH=/home/mark/bioinformatics/src/mreps:$PATH
export PATH=/home/mark/bioinformatics/src/RepeatMasker:$PATH

## Run each step individually with no hangups so that outputs/error messages can be inspected upon completion.
nohup TEannot.py -P repet_annot -C TEannot.cfg -S 1
nohup TEannot.py -P repet_annot -C TEannot.cfg -S 2 -a BLR
nohup TEannot.py -P repet_annot -C TEannot.cfg -S 2 -a RM
nohup TEannot.py -P repet_annot -C TEannot.cfg -S 2 -a BLR -r
nohup TEannot.py -P repet_annot -C TEannot.cfg -S 2 -a RM -r
nohup TEannot.py -P repet_annot -C TEannot.cfg -S 3 -c BLR+RM
nohup TEannot.py -P repet_annot -C TEannot.cfg -S 4 -s TRF
nohup TEannot.py -P repet_annot -C TEannot.cfg -S 4 -s Mreps
nohup TEannot.py -P repet_annot -C TEannot.cfg -S 4 -s RMSSR
nohup TEannot.py -P repet_annot -C TEannot.cfg -S 5
nohup TEannot.py -P repet_annot -C TEannot.cfg -S 6 -b tblastx
nohup TEannot.py -P repet_annot -C TEannot.cfg -S 7
nohup TEannot.py -P repet_annot -C TEannot.cfg -S 8 -o gameXML
nohup TEannot.py -P repet_annot -C TEannot.cfg -S 8 -o GFF3
