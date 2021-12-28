#!/usr/bin/env python3
# Python-based plotting of Houston Lightning Mapping Array data for next-gen HDWX
# Created 21 December 2021 by Sam Gardner <stgardner4@tamu.edu>

from os import path, listdir, remove
from pyxlma.lmalib.io import read as lma_read
from matplotlib import pyplot as plt
from matplotlib import image as mpimage
from matplotlib import colors
from cartopy import crs as ccrs
from cartopy import feature as cfeat
import numpy as np
from metpy.plots import USCOUNTIES
import pandas as pd
from datetime import datetime as dt, timedelta
from pathlib import Path
import json
import warnings
import logging

axExtent = [-99.5, -91, 26, 33.5]
def set_size(w, h, ax=None):
    if not ax: ax=plt.gca()
    l = ax.figure.subplotpars.left
    r = ax.figure.subplotpars.right
    t = ax.figure.subplotpars.top
    b = ax.figure.subplotpars.bottom
    figw = float(w)/(r-l)
    figh = float(h)/(t-b)
    ax.figure.set_size_inches(figw, figh)

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

def writeJson(productID, productPath, runPathExtension, validTime):
    # If you have no idea what's going on or why I'm doing all this json stuff, 
    # check out http://weather-dev.geos.tamu.edu/wx4stg/api/ for documentation
    # Get description and GIS based on productID
    if productID == 140:
        productDesc = "HLMA VHF 1-minute Sources"
        isGIS = True
        gisInfo = [str(axExtent[2])+","+str(axExtent[0]), str(axExtent[3])+","+str(axExtent[1])]
        fileExtension = "png"
        productFrameCount = 60
    elif productID == 141:
        productDesc = "HLMA VHF 1-minute Sources"
        isGIS = False
        gisInfo = ["0,0", "0,0"] # gisInfo is ["0,0", "0,0"] for non-GIS products
        fileExtension = "png"
        productFrameCount = 60
    elif productID == 143:
        productDesc = "HLMA VHF 10-minute Sources"
        isGIS = True
        gisInfo = [str(axExtent[2])+","+str(axExtent[0]), str(axExtent[3])+","+str(axExtent[1])]
        fileExtension = "png"
        productFrameCount = 60
    elif productID == 144:
        productDesc = "HLMA VHF 10-minute Sources"
        isGIS = False
        gisInfo = ["0,0", "0,0"]
        fileExtension = "png"
        productFrameCount = 60
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
        "isGIS" : isGIS,
        "fileExtension" : fileExtension
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
    # Now we need to add the frame we just wrote, as well as any that exist in the output directory that don't have metadata yet. 
    # To do this, we first check if the output directory is not empty.
    productRunPath = path.join(basePath, "output", productPath, runPathExtension)
    if len(listdir(productRunPath)) > 0:
        # If there are files inside, list them all
        frameNames = listdir(productRunPath)
        # get an array of integers representing the minutes past the hour of frames that have already been generated
        frameMinutes = [int(framename.replace(".png", "")) for framename in frameNames]
        # Loop through the previously-generated minutes and generate metadata for each
        for frameMin in frameMinutes:
            frmDict = {
                "fhour" : 0, # forecast hour is 0 for non-forecasts
                "filename" : str(frameMin)+".png",
                "gisInfo" : gisInfo,
                "valid" : int(validTime.strftime("%Y%m%d%H00"))+frameMin
            }
            # If this dictionary isn't already in the framesArray, add it
            if frmDict not in framesArray:
                framesArray.append(frmDict)
    productRunDict = {
        "publishTime" : publishTime.strftime("%Y%m%d%H%M"),
        "pathExtension" : runPathExtension,
        "runName" : validTime.strftime("%d %b %Y %HZ"),
        "availableFrameCount" : len(framesArray),
        "totalFrameCount" : productFrameCount,
        "productFrames" : sorted(framesArray, key=lambda dict: dict["valid"]) # productFramesArray, sorted by increasing valid Time
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

def makeSourcePlots(lmaFilePaths):
    # Silence error_bad_lines warning when reading in LMA data
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        # Read in LMA data
        lmaData, startTimeOfPlot = lma_read.dataset(lmaFilePaths)
    # Set parameters for 1 or 10 minute data
    if len(lmaFilePaths) == 1:
        gisProductID = 140
        staticProductID = 141
        writeToStatus("Plotting 1-minute data for "+startTimeOfPlot.strftime("%H:%M"))
    else:
        gisProductID = 143
        staticProductID = 144
        writeToStatus("Plotting 10-minute data for "+startTimeOfPlot.strftime("%H:%M"))
    # Create fig/ax
    fig = plt.figure()
    ax = plt.axes(projection=ccrs.epsg(3857))
    # We want to color our points by the time that they occurred, relative to the rest of the dataset
    times = [(pd.to_datetime(npt) - startTimeOfPlot).seconds + (pd.to_datetime(npt) - startTimeOfPlot).microseconds*0.000001 for npt in lmaData.event_time.data]
    # Normalize the colormap for the number of minutes we're plotting
    norm = colors.Normalize(0, 60*len(lmaFilePaths))
    # We want to "mask out" points where the event chi^2 is greater than 2
    chi2Mask = np.where(lmaData.event_chi2.data > 2.0, 1, 0)
    # Plot data
    vhfSct = ax.scatter(np.ma.masked_array(lmaData.event_longitude.data, mask=chi2Mask), np.ma.masked_array(lmaData.event_latitude.data, mask=chi2Mask), 3, np.ma.masked_array(times, mask=chi2Mask), transform=ccrs.PlateCarree(), zorder=4, cmap="rainbow", norm=norm)
    # Get a handle to a pixel
    px = 1/plt.rcParams["figure.dpi"]
    # Make the GIS plot have a decent resolution. This wont end up being exactly 1920 by 1080 
    # as it's overridden by the next line "set_extent", but should improve the resolution enough
    set_size(1920*px, 1080*px, ax=ax)
    # Now we set the extent to the coordinates we want
    ax.set_extent(axExtent, crs=ccrs.PlateCarree())
    # Add the number of one minute files ingested to the start time to get the end time
    timeOfPlot = startTimeOfPlot + timedelta(minutes=len(lmaFilePaths))
    # Create a path object to 'productPath' (as defined by the HDWX API), in this case gisproducts/hlma/vhf/ 
    # Use os.path to preserve Windows compatibility
    gisProductPath = path.join("gisproducts", "hlma", "vhf-"+str(len(lmaFilePaths))+"min")
    # Create a path oobject to 'runPathExtension', <year>/<month>/<day>/<hour>00/
    runPathExt = path.join(dt.strftime(timeOfPlot, "%Y"), dt.strftime(timeOfPlot, "%m"), dt.strftime(timeOfPlot, "%d"), dt.strftime(timeOfPlot, "%H")+"00")
    # Target path for the "GIS"/transparent image is output/<gisProductPath>/<runPathExtension>/<minute>.png
    gisSavePath = path.join(basePath, "output", gisProductPath, runPathExt, dt.strftime(timeOfPlot, "%M")+".png")
    # Create target directory if it doesn't already exist
    Path(path.dirname(gisSavePath)).mkdir(parents=True, exist_ok=True)
    # Get the exact extent of just the axes without the matplotlib auto-generated whitespace
    extent = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
    # save the figure, but trim the whitespace
    # we do this because including the whitespace would make the data not align to the GIS information in the metadata
    fig.savefig(gisSavePath, transparent=True, bbox_inches=extent)
    # Write metadata for the product
    writeJson(gisProductID, gisProductPath, runPathExt, timeOfPlot)
    # For the "static"/non-GIS/opaque image, add county/state/coastline borders
    logging.captureWarnings(True)
    ax.add_feature(USCOUNTIES.with_scale("5m"), edgecolor="gray", zorder=2)
    ax.add_feature(cfeat.STATES.with_scale("10m"), linewidth=0.5, zorder=3)
    ax.add_feature(cfeat.COASTLINE.with_scale("10m"), linewidth=0.5, zorder=3)
    # Move the data axes to maximize the amount of space available to it
    ax.set_position([0.05, 0.11, .9, .87])
    # Create a "colorbar axes". We want it to be beneath the plot, under the far left third.
    cbax = fig.add_axes([0, 0, (ax.get_position().width/3), .025])
    fig.colorbar(vhfSct, cax=cbax, orientation="horizontal", label="Seconds after "+startTimeOfPlot.strftime("%-d %b %Y %H%MZ"))
    cbax.set_position([0.05, ax.get_position().y0-.01-cbax.get_position().height, cbax.get_position().width, cbax.get_position().height])
    # Create a "title axes". We want it to be one third the width of the data axes, height will be handled by matplotlib automatically,
    # and we'll worry about positioning later
    tax = fig.add_axes([0,0,(ax.get_position().width/3),.05])
    # Add a descriptive title
    tax.text(0.5, 0.5, "Houston LMA "+str(len(lmaFilePaths))+"-minute VHF Sources\nValid "+timeOfPlot.strftime("%-d %b %Y %H%MZ"), horizontalalignment="center", verticalalignment="center", fontsize=16)
    # add credit
    tax.set_xlabel("Python HDWX -- Send bugs to stgardner4@tamu.edu")
    # Turn off the axis spines and ticks to give the appearance of floating text
    # we can't use 'tax.axis("off")' here as that would hide the xlabel
    plt.setp(tax.spines.values(), visible=False)
    tax.tick_params(left=False, labelleft=False)
    tax.tick_params(bottom=False, labelbottom=False)
    # Now we handle positioning, we want the plot to be centered below the data axes, as close to the bottom of the figure as possible
    # without cutting off the xlabel, and retain the original width/height
    tax.set_position([(ax.get_position().width/2)-(tax.get_position().width/2)+ax.get_position().x0,ax.get_position().y0-.01-tax.get_position().height,tax.get_position().width,tax.get_position().height], which="both")
    # Add a "logo axes" to display the TAMU ATMO logo. Again, one third the width of the data axes.
    lax = fig.add_axes([0,0,(ax.get_position().width/3),1])
    # The logo axes must have the same aspect ratio as the image we're trying to display or else the image will be stretched/compressed
    lax.set_aspect(2821/11071)
    # turn off the axis
    lax.axis("off")
    # turn off axis spines
    plt.setp(lax.spines.values(), visible=False)
    # read in the image
    atmoLogo = mpimage.imread("assets/atmoLogo.png")
    # show the image
    lax.imshow(atmoLogo)
    # We want the logo axes to be all the way to the right, and as low as possible without cutting anything off
    lax.set_position([.95-lax.get_position().width, ax.get_position().y0-.01-lax.get_position().height, lax.get_position().width, lax.get_position().height], which="both")
    # Make sure image is opaque
    fig.set_facecolor("white")
    # Set size to 1080p, resolution of the weather center monitors
    fig.set_size_inches(1920*px, 1080*px)
    # Create a path object to 'productPath' (as defined by the HDWX API), in this case gisproducts/hlma/vhf/ 
    staticProductPath = path.join("products", "hlma", "vhf-"+str(len(lmaFilePaths))+"min")
    # Target path for the "static"/non-GIS/transparent image is output/<gisProductPath>/<runPathExtension>/<minute>.png
    staticSavePath = path.join(basePath, "output", staticProductPath, runPathExt, dt.strftime(timeOfPlot, "%M")+".png")
    # Create save directory if it doesn't already exist
    Path(path.dirname(staticSavePath)).mkdir(parents=True, exist_ok=True)
    # Write the image
    fig.savefig(staticSavePath)
    # Write metadata for the product
    writeJson(staticProductID, staticProductPath, runPathExt, timeOfPlot)
    # Close figure when done (memory management)
    plt.close(fig)

if __name__ == "__main__":
    # get path to starting dir
    basePath = path.dirname(path.abspath(__file__))
    # get path to input files
    inputPath = path.join(basePath, "lightningin")
    # Get current time
    now = dt.utcnow()
    # Get time one hour ago
    oneHourAgo = now - timedelta(hours=1)
    # Create (empty, but add to it in a sec...) list representing already plotted frames
    alreadyPlottedOneMinFrames = list()
    # Get path to last hour's json metadata
    lastHourOneMinMetadataPath = path.join(basePath, "output", "metadata", "products", "140", dt.strftime(oneHourAgo, "%Y%m%d%H00")+".json")
    # Read in last hour's metadata
    if path.exists(lastHourOneMinMetadataPath):
        with open(lastHourOneMinMetadataPath, "r") as jsonRead:
            lastHourOneMinData = json.load(jsonRead)
        # Add already-generated frames to alreadyPlottedOneMinFrames list
        # The valid time gets converted first from an int to a string, then the string is trimmed to only include HHMM
        [alreadyPlottedOneMinFrames.append(str(frame["valid"])[-4:]+"00") for frame in lastHourOneMinData["productFrames"]]
    # Do the same thing for this hour's metadata
    thisHourMetadataPath = path.join(basePath, "output", "metadata", "products", "140", dt.strftime(now, "%Y%m%d%H00")+".json")
    if path.exists(thisHourMetadataPath):
        with open(thisHourMetadataPath, "r") as jsonRead:
            thisHourOneMinData = json.load(jsonRead)
        [alreadyPlottedOneMinFrames.append(str(frame["valid"])[-4:]+"00") for frame in thisHourOneMinData["productFrames"]]
    # Plot every file in the input directory
    inputDirContents = sorted(listdir(inputPath))
    for file in inputDirContents:
        timeOfFileArr = file.split("_")
        # The time in the filename is the *start*, but the time in the json is the end, so add one minute to the filename time
        timeOfFile = dt.strptime("20"+timeOfFileArr[1]+timeOfFileArr[2], "%Y%m%d%H%M%S") + timedelta(minutes=1)
        if timeOfFile.strftime("%H%M%S") in alreadyPlottedOneMinFrames:
            continue
        if timeOfFile < now - timedelta(hours=2):
            remove(path.join(inputPath, file))
            continue
        makeSourcePlots([path.join(inputPath, file)])
    # Now let's do the same thing, but for 10-minute intervals of data
    alreadyPlottedTenMinFrames = list()
    lastHourTenMinMetadataPath = path.join(basePath, "output", "metadata", "products", "143", dt.strftime(oneHourAgo, "%Y%m%d%H00")+".json")
    if path.exists(lastHourTenMinMetadataPath):
        with open(lastHourTenMinMetadataPath, "r") as jsonRead:
            lastHourTenMinData = json.load(jsonRead)
        [alreadyPlottedTenMinFrames.append(str(frame["valid"])[-4:]+"00") for frame in lastHourTenMinData["productFrames"]]
    thisHourTenMinMetadataPath = path.join(basePath, "output", "metadata", "products", "143", dt.strftime(now, "%Y%m%d%H00")+".json")
    if path.exists(thisHourTenMinMetadataPath):
        with open(thisHourTenMinMetadataPath, "r") as jsonRead:
            thisHourTenMinData = json.load(jsonRead)
        [alreadyPlottedTenMinFrames.append(str(frame["valid"])[-4:]+"00") for frame in thisHourTenMinData["productFrames"]]
    inputDirContents = sorted(listdir(inputPath))
    for i in range(10, len(inputDirContents)):
        lastFileInRange = inputDirContents[i]
        timeOfLastFileArr = lastFileInRange.split("_")
        timeOfLastFile = dt.strptime("20"+timeOfLastFileArr[1]+timeOfLastFileArr[2], "%Y%m%d%H%M%S") + timedelta(minutes=1)
        if timeOfLastFile.strftime("%H%M%S") in alreadyPlottedTenMinFrames:
            continue
        filesToPlot = [path.join(inputPath, fileToInclude) for fileToInclude in inputDirContents[(i-10):i]]
        makeSourcePlots(filesToPlot)
