#!/usr/bin/env python3
# Python-based placefile creation of Houston Lightning Mapping Array data for python-based HDWX
# Created 21 December 2021 by Sam Gardner <stgardner4@tamu.edu>

from os import path, listdir, chmod
from pyxlma.lmalib.io import read as lma_read
import pandas as pd
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
    timeOfPlot = startTimeOfPlot + timedelta(minutes=len(lmaFilePaths))
    # Create a GR2Analyst compatible placefile
    if len(lmaFilePaths) == 1:
        placeFileString = ";One-minute HLMA Data\n;Generated by pyxlma on python-based HDWX\n;Installation instructions: GR2Analyst -> Window -> Placefile manager -> + -> Add the URL to this file\n;Valid: "+dt.strftime(timeOfPlot, "%Y-%m-%d %H:%M")+"\n;Code by Sam Gardner <stgardner4@tamu.edu>\nThreshold: 999\nRefresh: 1\nIconFile: 1, 32, 32, 16, 16, https://hdwx.tamu.edu/products/wxgen3/gr2a/lightningicons.png\n"
        productID = 142
    else:
        placeFileString = ";Ten-minute HLMA Data\n;Generated by pyxlma on python-based HDWX\n;Installation instructions: GR2Analyst -> Window -> Placefile manager -> + -> Add the URL to this file\n;Valid: "+dt.strftime(timeOfPlot, "%Y-%m-%d %H:%M")+"\n;Code by Sam Gardner <stgardner4@tamu.edu>\nThreshold: 999\nRefresh: 1\nIconFile: 1, 32, 32, 16, 16, https://hdwx.tamu.edu/products/wxgen3/gr2a/lightningicons.png\n"
        productID = 145
    # Add all of our lightning points
    for lat, lon, time, alt, chi2, power, stations in zip(lmaData.event_latitude.data, lmaData.event_longitude.data, lmaData.event_time.data, lmaData.event_altitude.data, lmaData.event_chi2.data, lmaData.event_power.data, lmaData.event_stations.data):
        # We want to filter out any source with a chi^2 over 50.
        if chi2 <= 1:
            if alt <= 20000:
                scaleOfPoint = 1+int(float(time - lmaData.event_time.data[0])/float(lmaData.event_time.data[-1] - lmaData.event_time.data[0])*1024)
                # Add an icon for every vhf source
                placeFileString = placeFileString+f"Icon: {lat}, {lon}, 000, 1, {scaleOfPoint}, {str(time)[11:]}\\nAltitude (m): {alt}\\nReduced chi^2: {chi2}\\nPower (dBW): {power}\\nStation Count: {stations}\n"
    # Create a path object for GR2A placefile's productPath
    gr2aProductPath = path.join("gr2a", "")
    # Target path for gr2a placefiles is output/gr2a/1min-src.txt
    gr2aSavePath = path.join(basePath, "output", gr2aProductPath, str(len(lmaFilePaths))+"min-src.txt")
    # Create parent dir
    Path(path.dirname(gr2aSavePath)).mkdir(parents=True, exist_ok=True)
    # write out the placefile
    with open(gr2aSavePath, "w") as textWrite:
        textWrite.write(placeFileString)
    # Also write out php file for the website ( this is necessary because pointing browser at a txt will result in caching)
    phpText = "<?php\n    header('Content-type: text/plain');\n    echo file_get_contents(\"./"+str(len(lmaFilePaths))+"min-src.txt\")\n?>"
    with open(path.join(basePath, "output", gr2aProductPath, str(len(lmaFilePaths))+"min-src.php"), "w") as phpWrite:
        phpWrite.write(phpText)
    if not path.exists(path.join(basePath, "output", gr2aProductPath, "lightningicons.png")):
        if path.exists(path.join(basePath, "assets", "icons.png")):
            copyfile(path.join(basePath, "assets", "icons.png"), path.join(basePath, "output", gr2aProductPath, "lightningicons.png"))
    # write metadata for gr2a file
    if hasHelpers:
        HDWX_helpers.writeJson(basePath, productID, timeOfPlot, str(len(lmaFilePaths))+"min-src.txt", timeOfPlot, ["0,0", "0,0"], 60)

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
    # Get the latest filename
    lastOneMinFile = filesToPlot[-1]
    # Get the time the last one-minute lma datafile as a datetime object
    lastOneMinFileDt = dt.strptime(lastOneMinFile.split("_")[1][-4:]+lastOneMinFile.split("_")[2], "%m%d%H%M%S").replace(year=dt.utcnow().year)
    # Check to see if the one-minute placefile has the latest input file generated
    # Read in metadata for one minute placefile
    oneMinMetadataPath = path.join(basePath, "output", "metadata", "products", "142", lastOneMinFileDt.strftime("%Y%m%d%H00")+".json")
    if path.exists(oneMinMetadataPath):
        with open(oneMinMetadataPath, "r") as jsonRead:
            lastOneMinMetadata = json.load(jsonRead)
        # If the latest lma input file is newer than the last generated placefile, generate a new one
        if int(lastOneMinMetadata["productFrames"][0]["valid"]) < int(lastOneMinFileDt.strftime("%Y%m%d%H%M")):
            makeSrcPlacefile([path.join(inputPath, lastOneMinFile)])
    else:
        # If the json doesn't exist, then we definitely need to plot the latest file.
        makeSrcPlacefile([path.join(inputPath, lastOneMinFile)])
    # Now let's make a placefile for the last ten minutes of data
    lastTenMinsFiles = filesToPlot[-10:]
    tenMinMetadataPath = path.join(basePath, "output", "metadata", "products", "145", lastOneMinFileDt.strftime("%Y%m%d%H00")+".json")
    if path.exists(tenMinMetadataPath):
        with open(tenMinMetadataPath, "r") as jsonRead:
            lastTenMinMetadata = json.load(jsonRead)
        if int(lastTenMinMetadata["productFrames"][0]["valid"]) < int(lastOneMinFileDt.strftime("%Y%m%d%H%M")):
            makeSrcPlacefile([path.join(inputPath, tenMinFile) for tenMinFile in lastTenMinsFiles])
    else:
        makeSrcPlacefile([path.join(inputPath, tenMinFile) for tenMinFile in lastTenMinsFiles])
