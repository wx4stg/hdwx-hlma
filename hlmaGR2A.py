#!/usr/bin/env python3
# Python-based placefile creation of Houston Lightning Mapping Array data for next-gen HDWX
# Created 21 December 2021 by Sam Gardner <stgardner4@tamu.edu>

from os import path, listdir
from pyxlma.lmalib.io import read as lma_read
import pandas as pd
from datetime import datetime as dt, timedelta
from pathlib import Path
import json
from shutil import copyfile
import warnings

def writeToStatus(stringToWrite):
    print(stringToWrite)
    stringToWrite = stringToWrite+"\n"
    if path.exists(path.join(basePath, "status.txt")):
        currentStatusFile = open(path.join(basePath, "status.txt"), "r")
        currentStr = open(path.join(basePath, "status.txt"), "r").read()
        currentStatusFile.close()
    else:
        currentStr = ""
    if stringToWrite not in currentStr:
        with open(path.join(basePath, "status.txt"), "a") as statw:
            statw.write(stringToWrite)
            statw.close()

def writeJson(productID, productPath, validTime):
    # If you have no idea what's going on or why I'm doing all this json stuff, 
    # check out http://weather-dev.geos.tamu.edu/wx4stg/api/ for documentation
    # Get description and GIS based on productID
    if productID == 142:
        productDesc = "GR2Analyst HLMA VHF Sources (1 minute)"
        filename = "1min-src.php"
    elif productID == 145:
        productDesc = "GR2Analyst HLMA VHF Sources (10 minutes)"
        filename = "10min-src.php"
    # For prettyness' sake, make all the publishTimes the same
    publishTime = dt.utcnow()
    # Create dictionary for the product. 
    productDict = {
        "productID" : productID,
        "productDescription" : productDesc,
        "productPath" : productPath,
        "productReloadTime" : 60,
        "lastReloadTime" : publishTime.strftime("%Y%m%d%H%M"),
        "isForecast" : True,
        "isGIS" : False,
        "fileExtension" : "php"
    }
    # Target path for the product json is just output/metadata/<productID>.json
    productDictJsonPath = path.join(basePath, "output", "metadata", str(productID)+".json")
    # Create output/metadata/ if it doesn't already exist
    Path(path.dirname(productDictJsonPath)).mkdir(parents=True, exist_ok=True)
    with open(productDictJsonPath, "w") as jsonWrite:
        # Write the json. indent=4 gives pretty/human-readable format
        json.dump(productDict, jsonWrite, indent=4)
    # Now we need to write a json for the product run in output/metadata/products/<productID>/<runTime>.json
    productRunDictPath = path.join(basePath, "output", "metadata", "products", str(productID), validTime.strftime("%Y%m%d%H00")+".json")
    # Create parent directory if it doesn't already exist.
    Path(path.dirname(productRunDictPath)).mkdir(parents=True, exist_ok=True)
    # If the json file already exists, read it in to to discover which frames have already been generated
    if path.exists(productRunDictPath):
        with open(productRunDictPath, "r") as jsonRead:
            oldData = json.load(jsonRead)
        # Add previously generated frames to a list, framesArray
        framesArray = oldData["productFrames"]
    else:
        # If that file didn't exist, then create an empty list instead
        framesArray = list()
    # Create a frame array to add to the runDict
    framesArray = [{
        "fhour" : 0,
        "filename" : filename,
        "gisInfo" : ["0,0", "0,0"],
        "valid" : int(validTime.strftime("%Y%m%d%H%M"))
    }]
    # Create a dictionary for the run
    productRunDict = {
        "publishTime" : publishTime.strftime("%Y%m%d%H%M"),
        "pathExtension" : "",
        "runName" : validTime.strftime("%d %b %Y %HZ"),
        "availableFrameCount" : 1,
        "totalFrameCount" : 1,
        "productFrames" : framesArray
    }
    # Write productRun dictionary to json
    with open(productRunDictPath, "w") as jsonWrite:
        json.dump(productRunDict, jsonWrite, indent=4)
    # Now we need to create a dictionary for the product type (TAMU)
    productTypeID = 1
    # Output for this json is output/metadata/productTypes/1.json
    productTypeDictPath = path.join(basePath, "output/metadata/productTypes/"+str(productTypeID)+".json")
    # Create output directory if it doesn't already exist
    Path(path.dirname(productTypeDictPath)).mkdir(parents=True, exist_ok=True)
    # Create empty list that will soon hold a dict for each of the products generated by this script
    productsInType = list()
    # If the productType json file already exists, read it in to discover which products it contains
    if path.exists(productTypeDictPath):
        with open(productTypeDictPath, "r") as jsonRead:
            oldProductTypeDict = json.load(jsonRead)
        # Add all of the products from the json file into the productsInType list...
        for productInOldDict in oldProductTypeDict["products"]:
            # ...except for the one that's currently being generated (prevents duplicating it)
            if productInOldDict["productID"] != productID:
                productsInType.append(productInOldDict)
    # Add the productDict for the product we just generated
    productsInType.append(productDict)
    # Create productType Dict
    productTypeDict = {
        "productTypeID" : productTypeID,
        "productTypeDescription" : "TAMU",
        "products" : sorted(productsInType, key=lambda dict: dict["productID"]) # productsInType, sorted by productID
    }
    # Write productType dict to json
    with open(productTypeDictPath, "w") as jsonWrite:
        json.dump(productTypeDict, jsonWrite, indent=4)

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
        placeFileString = ";One-minute HLMA Data\n;Generated by pyxlma on next-gen HDWX\n;Installation instructions: GR2Analyst -> Window -> Placefile manager -> + -> Add the URL to this file\n;Valid: "+dt.strftime(timeOfPlot, "%Y-%m-%d %H:%M")+"\n;Code by Sam Gardner <stgardner4@tamu.edu>\nThreshold: 999\nRefresh: 1\nIconFile: 1, 4, 4, 2, 2, http://weather-dev.geos.tamu.edu/wx4stg/gr2a/lightningicons.png\n"
        productID = 142
        writeToStatus("Generating 1-minute placefile for "+timeOfPlot.strftime("%H:%M"))
    else:
        placeFileString = ";Ten-minute HLMA Data\n;Generated by pyxlma on next-gen HDWX\n;Installation instructions: GR2Analyst -> Window -> Placefile manager -> + -> Add the URL to this file\n;Valid: "+dt.strftime(timeOfPlot, "%Y-%m-%d %H:%M")+"\n;Code by Sam Gardner <stgardner4@tamu.edu>\nThreshold: 999\nRefresh: 1\nIconFile: 1, 4, 4, 2, 2, http://weather-dev.geos.tamu.edu/wx4stg/gr2a/lightningicons.png\n"
        productID = 145
        writeToStatus("Generating 10-minute placefile for "+timeOfPlot.strftime("%H:%M"))
    # Add all of our lightning points
    for lat, lon, time, alt, chi2, power, stations in zip(lmaData.event_latitude.data, lmaData.event_longitude.data, lmaData.event_time.data, lmaData.event_altitude.data, lmaData.event_chi2.data, lmaData.event_power.data, lmaData.event_stations.data):
        # We want to filter out any source with a chi^2 over 50.
        if chi2 <= 1:
            if alt <= 20000:
                time = pd.to_datetime(time)
                # Get time elapsed since start of file range until the point
                timeOfPoint = time - startTimeOfPlot
                # Convert that time to a float
                timeOfPoint = timeOfPoint.seconds + timeOfPoint.microseconds*0.000001
                # There are 1020 icons in lightningicons.png, so based on the time scale, pick which one to use
                scaleOfPoint = 1020 - int(1020*timeOfPoint/(60*len(lmaFilePaths))) + 1
                # Add an icon for every vhf source
                placeFileString = placeFileString+"Icon: "+str(lat)+", "+str(lon)+", 000, 1, "+str(scaleOfPoint)+", "+time.strftime("%H:%M:%S.%f")+r"\nAltitude (m): "+str(alt)+r"\nReduced chi^2: "+str(chi2)+r"\nPower (dBW): "+str(power)+r"\nStation Count: "+str(stations)+"\n"
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
    writeJson(productID, gr2aProductPath, timeOfPlot)

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
        if lastOneMinMetadata["productFrames"][0]["valid"] < int(lastOneMinFileDt.strftime("%Y%m%d%H%M")):
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
        if lastTenMinMetadata["productFrames"][0]["valid"] < int(lastOneMinFileDt.strftime("%Y%m%d%H%M")):
            makeSrcPlacefile([path.join(inputPath, tenMinFile) for tenMinFile in lastTenMinsFiles])
    else:
        makeSrcPlacefile([path.join(inputPath, tenMinFile) for tenMinFile in lastTenMinsFiles])
