#!/bin/bash
# HLMA product generation script for python-based HDWX
# Created 27 December 2021 by Sam Gardner <stgardner4@tamu.edu>

if [ ! -d output/ ]
then
    mkdir output/
fi
if [ ! -d lightningin/ ]
then
    mkdir lightningin/
fi
if [ -f ../config.txt ]
then
    source ../config.txt
else
    condaEnvName="HDWX"
fi
if [ -f ../HDWX_helpers.py ]
then
    if [ -f ./HDWX_helpers.py ]
    then
        rm ./HDWX_helpers.py
    fi
    cp ../HDWX_helpers.py ./
fi
if [ -f $condaRootPath/envs/$condaEnvName/bin/python3 ]
then
    $condaRootPath/envs/$condaEnvName/bin/python3 hlmaGR2A.py
    bash continuous.sh
    $condaRootPath/envs/$condaEnvName/bin/python3 cleanup.py
fi

