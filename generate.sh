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
if [ -f ../config.txt ]
then
    source ../config.txt
else
    condaEnvName="HDWX"
fi

if [ -f $condaRootPath/envs/$condaEnvName/bin/python3 ]
then
    $condaRootPath/envs/$condaEnvName/bin/python3 hlmaGR2A.py
    if [ -f flash1min-lock.txt ]
    then
        pidToCheck=`cat flash1min-lock.txt`
        if ! kill -0 $pidToCheck
        then
            echo "starting 1 min flash"
            $condaRootPath/envs/$condaEnvName/bin/python3 hlmaPlot.py flash 1 &
            echo -n $! > flash1min-lock.txt
        else
            echo "1 minute flash locked"
        fi
    else
            echo "starting 1 min flash"
            $condaRootPath/envs/$condaEnvName/bin/python3 hlmaPlot.py flash 1 &
            echo -n $! > flash1min-lock.txt
    fi
    
    if [ -f src1min-lock.txt ]
    then
        pidToCheck=`cat src1min-lock.txt`
        if ! kill -0 $pidToCheck
        then
            $condaRootPath/envs/$condaEnvName/bin/python3 hlmaPlot.py src 1 &
            echo -n $! > src1min-lock.txt
        else
            echo "1 minute source locked"
        fi
    else
            echo "starting 1 min source"
            $condaRootPath/envs/$condaEnvName/bin/python3 hlmaPlot.py src 1 &
            echo -n $! > src1min-lock.txt
    fi


    if [ -f flash10min-lock.txt ]
    then
        pidToCheck=`cat flash10min-lock.txt`
        if ! kill -0 $pidToCheck
        then
            $condaRootPath/envs/$condaEnvName/bin/python3 hlmaPlot.py flash 10 &
            echo -n $! > flash10min-lock.txt
        else
            echo "10 minute flash locked"
        fi
    else
            echo "starting 10 min flash"
            $condaRootPath/envs/$condaEnvName/bin/python3 hlmaPlot.py flash 10 &
            echo -n $! > flash10min-lock.txt
    fi


    if [ -f src10min-lock.txt ]
    then
        pidToCheck=`cat src10min-lock.txt`
        if ! kill -0 $pidToCheck
        then
            $condaRootPath/envs/$condaEnvName/bin/python3 hlmaPlot.py src 10 &
            echo -n $! > src10min-lock.txt
        else
            echo "10 minute source locked"
        fi
    else
            echo "starting 10 min source"
            $condaRootPath/envs/$condaEnvName/bin/python3 hlmaPlot.py src 10 &
            echo -n $! > src10min-lock.txt
    fi
fi
