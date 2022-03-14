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
if [ -f ../config.txt ]
then
    source ../config.txt
else
    condaEnvName="HDWX"
fi

if [ -f $condaRootPath/envs/$condaEnvName/bin/python3 ]
then
    $condaRootPath/envs/$condaEnvName/bin/python3 hlmaGR2A.py
    $condaRootPath/envs/$condaEnvName/bin/python3 hlmaPlot.py
    $condaRootPath/envs/$condaEnvName/bin/python3 cleanup.py
fi
