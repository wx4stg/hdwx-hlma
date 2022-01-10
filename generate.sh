#!/bin/bash
# HLMA product generation script for next-gen HDWX
# Created 27 December 2021 by Sam Gardner <stgardner4@tamu.edu>

if [ ! -d output/ ]
then
    mkdir output/
fi
if [ ! -d lightningin/ ]
then
    mkdir lightningin/
fi

if [ -f status.txt ]
then
  echo "lockfile found, exiting"
  exit
fi
touch status.txt

$CONDA_PREFIX/bin/python3 hlmaGR2A.py
$CONDA_PREFIX/bin/python3 hlmaPlot.py

rm status.txt
