#!/usr/bin/env python3
# Python-based placefile creation of Houston Lightning Mapping Array data for python-based HDWX
# Created 21 December 2021 by Sam Gardner <stgardner4@tamu.edu>

import hlmaFetch
from os import path, listdir
from pyxlma.lmalib.io import read as lma_read
import pandas as pd
import numpy as np
from datetime import datetime as dt, timedelta
from pathlib import Path
import json
from shutil import copyfile
import warnings
basePath = path.abspath(path.dirname(__file__))
hasHelpers = False
if path.exists(path.join(basePath, "HDWX_helpers.py")):
    import HDWX_helpers
    hasHelpers = True


def makeSrcPlacefile(lmaFilePaths):
    # Silence error_bad_lines warning when reading in LMA data
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        # Read in LMA data
        lmaData, startTimeOfPlot = lma_read.dataset(lmaFilePaths)
    # Get time of the plot for metadata
    lastFileName = path.basename(lmaFilePaths[-1])
    timeOfPlot = dt.strptime(lastFileName.split("_")[1]+lastFileName.split("_")[2], "%y%m%d%H%M%S")+timedelta(seconds=int(lastFileName.split("_")[3].replace(".dat", "").replace(".gz", "")))
    elapsedTimeOfPlot = timeOfPlot - startTimeOfPlot
    # Create a GR2Analyst compatible placefile
    if elapsedTimeOfPlot.total_seconds() <= 90:
        outFileName = "1min-src"
        placeFileString = ";One-minute HLMA Data\n;Generated by pyxlma on python-based HDWX\n;Installation instructions: GR2Analyst -> Window -> Placefile manager -> + -> Add the URL to this file\n;Valid: "+dt.strftime(timeOfPlot, "%Y-%m-%d %H:%M")+"\n;Code by Sam Gardner <stgardner4@tamu.edu>\nThreshold: 999\nRefresh: 1\nIconFile: 1, 32, 32, 16, 16, https://hdwx.tamu.edu/products/wxgen3/gr2a/lightningicons.png\n"
        productID = 142
    else:
        outFileName = "10min-src"
        placeFileString = ";Ten-minute HLMA Data\n;Generated by pyxlma on python-based HDWX\n;Installation instructions: GR2Analyst -> Window -> Placefile manager -> + -> Add the URL to this file\n;Valid: "+dt.strftime(timeOfPlot, "%Y-%m-%d %H:%M")+"\n;Code by Sam Gardner <stgardner4@tamu.edu>\nThreshold: 999\nRefresh: 1\nIconFile: 1, 32, 32, 16, 16, https://hdwx.tamu.edu/products/wxgen3/gr2a/lightningicons.png\n"
        productID = 145
    # Add all of our lightning points
    point_mask = np.where((lmaData.event_chi2.data <= 1) & (lmaData.event_altitude.data <= 20000))[0]
    numRows = len(point_mask)
    if numRows > 0:
        lmaData = lmaData.isel(number_of_events=point_mask)
        if numRows > 1:
            scaleOfPoint = 1+((lmaData.event_time.data - lmaData.event_time.data[0]).astype(float)/float(lmaData.event_time.data[-1] - lmaData.event_time.data[0])*1024).astype(int)
        else:
            scaleOfPoint = np.array([1])
        outArr0 = np.char.add(np.full((numRows, 1), "Icon: ", dtype="|S6").astype(str), lmaData.event_latitude.data.astype(str).reshape(numRows, 1))
        outArr0 = np.char.add(outArr0, np.full((numRows, 1), ", ", dtype="|S2").astype(str))
        outArr0 = np.char.add(outArr0, lmaData.event_longitude.data.astype(str).reshape(numRows, 1))
        outArr0 = np.char.add(outArr0, np.full((numRows, 1), ", 000, 1, ", dtype="|S10").astype(str))
        outArr0 = np.char.add(outArr0, scaleOfPoint.astype(str).reshape(numRows, 1))
        outArr0 = np.char.add(outArr0, np.full((numRows, 1), ", ", dtype="|S2").astype(str))
        outArr0 = np.char.add(outArr0, np.char.replace(lmaData.event_time.data.astype(str), dt.utcnow().strftime("%Y-%m-%dT"), "").reshape(numRows, 1))
        outArr0 = np.char.add(outArr0, np.full((numRows, 1), "\\nAltitude (m):", dtype="|S16").astype(str))
        outArr0 = np.char.add(outArr0, lmaData.event_altitude.data.astype(str).reshape(numRows, 1))
        outArr0 = np.char.add(outArr0, np.full((numRows, 1), "\\nReduced chi^2:", dtype="|S17").astype(str))
        outArr0 = np.char.add(outArr0, lmaData.event_chi2.data.astype(str).reshape(numRows, 1))
        outArr0 = np.char.add(outArr0, np.full((numRows, 1), "\\nPower (dBW):", dtype="|S15").astype(str))
        outArr0 = np.char.add(outArr0, lmaData.event_power.data.astype(str).reshape(numRows, 1))
        outArr0 = np.char.add(outArr0, np.full((numRows, 1), "\\nStation Count:", dtype="|S17").astype(str))
        outArr0 = np.char.add(outArr0, lmaData.event_stations.data.astype(str).reshape(numRows, 1)).flatten()
        placeFileString = placeFileString+"\n".join(outArr0)
    print(placeFileString)
    # Create a path object for GR2A placefile's productPath
    gr2aProductPath = path.join("gr2a", "")
    # Target path for gr2a placefiles is output/gr2a/1min-src.txt
    gr2aSavePath = path.join(basePath, "output", gr2aProductPath, outFileName+".txt")
    # Create parent dir
    Path(path.dirname(gr2aSavePath)).mkdir(parents=True, exist_ok=True)
    # write out the placefile
    with open(gr2aSavePath, "w") as textWrite:
        textWrite.write(placeFileString)
    # Also write out php file for the website ( this is necessary because pointing browser at a txt will result in caching)
    phpText = "<?php\n    header('Content-type: text/plain');\n    echo file_get_contents(\"./"+outFileName+".txt\")\n?>"
    with open(path.join(basePath, "output", gr2aProductPath, outFileName+".php"), "w") as phpWrite:
        phpWrite.write(phpText)
    if not path.exists(path.join(basePath, "output", gr2aProductPath, "lightningicons.png")):
        if path.exists(path.join(basePath, "assets", "icons.png")):
            copyfile(path.join(basePath, "assets", "icons.png"), path.join(basePath, "output", gr2aProductPath, "lightningicons.png"))
    # write metadata for gr2a file
    if hasHelpers:
        HDWX_helpers.writeJson(basePath, productID, timeOfPlot, outFileName+".txt", timeOfPlot, ["0,0", "0,0"], 60)

if __name__ == "__main__":
    # get path to starting dir
    basePath = path.dirname(path.abspath(__file__))
    # get path to input files
    inputPath = path.join(basePath, "lightningin")
    # List the files in the input directory
    filesToPlot = sorted(listdir(inputPath))
    # If there are no files in the input dir, exit immediately
    if len(filesToPlot) == 0:
        exit()
    # Get the time of the latest complete minute of data
    for file in reversed(filesToPlot):
        lastMinFileDt = dt.strptime(file.split("_")[1][-4:]+file.split("_")[2], "%m%d%H%M%S").replace(year=dt.utcnow().year)
        if lastMinFileDt.second == 0:
            break
    # Get paths to those data
    oneMinuteFiles = hlmaFetch.getLmaFilesBetweenTimes(lastMinFileDt-timedelta(minutes=1), lastMinFileDt, True)
    print(oneMinuteFiles)
    tenMinuteFiles = hlmaFetch.getLmaFilesBetweenTimes(lastMinFileDt-timedelta(minutes=10), lastMinFileDt, True)
    # Check to see if the one-minute placefile has the latest input file generated
    # Read in metadata for one minute placefile
    oneMinMetadataPath = path.join(basePath, "output", "metadata", "products", "142", lastMinFileDt.strftime("%Y%m%d%H00")+".json")
    if path.exists(oneMinMetadataPath):
        with open(oneMinMetadataPath, "r") as jsonRead:
            lastOneMinMetadata = json.load(jsonRead)
        # If the latest lma input file is newer than the last generated placefile, generate a new one
        if int(lastOneMinMetadata["productFrames"][0]["valid"]) < int(lastMinFileDt.strftime("%Y%m%d%H%M")):
            makeSrcPlacefile(oneMinuteFiles)
    else:
        # If the json doesn't exist, then we definitely need to plot the latest file.
        makeSrcPlacefile(oneMinuteFiles)
    # Now let's make a placefile for the last ten minutes of data
    tenMinMetadataPath = path.join(basePath, "output", "metadata", "products", "145", lastMinFileDt.strftime("%Y%m%d%H00")+".json")
    if path.exists(tenMinMetadataPath):
        with open(tenMinMetadataPath, "r") as jsonRead:
            lastTenMinMetadata = json.load(jsonRead)
        if int(lastTenMinMetadata["productFrames"][0]["valid"]) < int(lastMinFileDt.strftime("%Y%m%d%H%M")):
            makeSrcPlacefile(tenMinuteFiles)
    else:
        makeSrcPlacefile(tenMinuteFiles)
