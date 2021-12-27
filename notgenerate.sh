#!/bin/bash
# HLMA product generation script for next-gen HDWX
# Created 27 December 2021 by Sam Gardner <stgardner4@tamu.edu>

if [ ! -d output/ ]
then
    mkdir output/
fi

if [ -f status.txt ]
then
  echo "lockfile found, exiting"
  exit
fi

if [ -f ~/mambaforge/envs/HDWX/bin/python3 ]
then
    ~/mambaforge/envs/HDWX/bin/python3 hlmaPlot.py
    ~/mambaforge/envs/HDWX/bin/python3 hlmaGR2A.py
fi
if [ -f ~/miniconda3/envs/HDWX/bin/python3 ]
then
    ~/mambaforge/envs/HDWX/bin/python3 hlmaPlot.py
    ~/mambaforge/envs/HDWX/bin/python3 hlmaGR2A.py
fi
rm status.txt