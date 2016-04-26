#!/usr/bin/env bash

## This script was used to run the REPET TE de novo pipeline.

## Set environment variables.
export REPET_PATH=/home/mark/bioinformatics/src/REPET_linux-x64-2.2/
export PATH=/home/mark/bioinformatics/src/REPET_linux-x64-2.2/bin:$PATH
export PATH=/home/mark/bioinformatics/src/REPET_linux-x64-2.2/denovo_pipe:$PATH
export PATH=/home/mark/bioinformatics/src/gt-1.5.8-Linux_x86_64-64bit-complete/bin:$PATH
export PATH=/home/mark/bioinformatics/src/RECON-1.07/scripts/:$PATH
export PATH=/home/mark/bioinformatics/src/hmmer-3.1b2-linux-intel-x86_64/binaries:$PATH

## Run each step individually with no hangups so that outputs/error messages can be inspected upon completion.
nohup TEdenovo.py -P repet -C TEdenovo.cfg -S 1
nohup TEdenovo.py -P repet -C TEdenovo.cfg -S 2 -s Blaster
nohup TEdenovo.py -P repet -C TEdenovo.cfg -S 2 --struct
nohup TEdenovo.py -P repet -C TEdenovo.cfg -S 3 -s Blaster -c Grouper
nohup TEdenovo.py -P repet -C TEdenovo.cfg -S 3 --struct
nohup TEdenovo.py -P repet -C TEdenovo.cfg -S 4 -s Blaster -c Grouper -m Map
nohup TEdenovo.py -P repet -C TEdenovo.cfg -S 4 --struct -m Map
nohup TEdenovo.py -P repet -C TEdenovo.cfg -S 5 -s Blaster -c Grouper -m Map --struct
nohup TEdenovo.py -P repet -C TEdenovo.cfg -S 6 -s Blaster -c Grouper -m Map --struct
nohup TEdenovo.py -P repet -C TEdenovo.cfg -S 7 -s Blaster -c Grouper -m Map --struct
nohup TEdenovo.py -P repet -C TEdenovo.cfg -S 8 -s Blaster -c Grouper -m Map -f Blastclust --struct

