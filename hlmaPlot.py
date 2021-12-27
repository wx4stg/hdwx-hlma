#!/usr/bin/env python3
# Python-based plotting of Houston Lightning Mapping Array data for next-gen HDWX
# Created 21 December 2021 by Sam Gardner <stgardner4@tamu.edu>

from os import path, listdir, remove
from pyxlma.lmalib.io import read as lma_read
from matplotlib import pyplot as plt
from matplotlib import image as mpimage
from cartopy import crs as ccrs
from cartopy import feature as cfeat
from metpy.plots import USCOUNTIES
import pandas as pd
from datetime import datetime as dt, timedelta
from pathlib import Path
import json


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

def writeJson(productID, productPath, runPathExtension, validTime):
    # If you have no idea what's going on or why I'm doing all this json stuff, 
    # check out http://weather-dev.geos.tamu.edu/wx4stg/api/ for documentation
    # Get description and GIS based on productID
    if productID == 140:
        productDesc = "HLMA VHF Sources"
        isGIS = True
        gisInfo = [str(axExtent[2])+","+str(axExtent[0]), str(axExtent[3])+","+str(axExtent[1])]
        fileExtension = "png"
        productFrameCount = 60
    elif productID == 141:
        productDesc = "HLMA VHF Sources"
        isGIS = False
        gisInfo = ["0,0", "0,0"] # gisInfo is ["0,0", "0,0"] for non-GIS products
        fileExtension = "png"
        productFrameCount = 60
    elif productID == 142:
        productDesc = "GR2Analyst HLMA VHF Sources"
        isGIS = True
        gisInfo = [str(axExtent[2])+","+str(axExtent[0]), str(axExtent[3])+","+str(axExtent[1])]
        fileExtension = "php"
        productFrameCount = 1
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

def makeOneMinPlots(lmaFilePath):
    # Read in LMA data
    lmaData = lma_read.lmafile(lmaFilePath).readfile()
    # Create fig/ax
    fig = plt.figure()
    ax = plt.axes(projection=ccrs.epsg(3857))
    # Plot Data
    ax.scatter(lmaData["lon"], lmaData["lat"], 1, "red", transform=ccrs.PlateCarree(), zorder=4)
    # Get a handle to a pixel
    px = 1/plt.rcParams["figure.dpi"]
    # Make the GIS plot have a decent resolution. This wont end up being exactly 1920 by 1080 
    # as it's overridden by the next line "set_extent", but should improve the resolution enough
    set_size(1920*px, 1080*px, ax=ax)
    # Now we set the extent to the coordinates we want
    ax.set_extent(axExtent, crs=ccrs.PlateCarree())
    # Extract the start time of the plot
    startTimeOfPlot = lmaData["Datetime"][0].to_pydatetime()
    # Since this is a one-minute file, the end time is just the start time plus one minute
    timeOfPlot = startTimeOfPlot + timedelta(minutes=1)
    # Create a path object to 'productPath' (as defined by the HDWX API), in this case gisproducts/hlma/vhf/ 
    # Use os.path to preserve Windows compatibility
    gisProductPath = path.join("gisproducts", "hlma", "vhf")
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
    writeJson(140, gisProductPath, runPathExt, timeOfPlot)
    # For the "static"/non-GIS/opaque image, add county/state/coastline borders
    ax.add_feature(USCOUNTIES.with_scale("5m"), edgecolor="gray", zorder=2)
    ax.add_feature(cfeat.STATES.with_scale("50m"), linewidth=0.5, zorder=3)
    ax.add_feature(cfeat.COASTLINE.with_scale("50m"), linewidth=0.5, zorder=3)
    # Move the data axes to maximize the amount of space available to it
    ax.set_position([0.05, 0.11, .9, .87])
    # Create a "title axes". We want it to be one third the width of the data axes, height will be handled by matplotlib automatically,
    # and we'll worry about positioning later
    tax = fig.add_axes([0,0,(ax.get_position().width/3),.05])
    # Add a descriptive title
    tax.text(0.5, 0.5, "Houston LMA 1-minute VHF Sources\nValid "+startTimeOfPlot.strftime("%-d %b %Y %H%MZ")+"through"+timeOfPlot.strftime("%-d %b %Y %H%MZ"), horizontalalignment="center", verticalalignment="center", fontsize=16)
    # add credit
    tax.set_xlabel("Python HDWX -- Send bugs to stgardner4@tamu.edu")
    # Turn off the axis spines and ticks to give the appearance of floating text
    # we can't use 'tax.axis("off")' here as that would hide the xlabel
    plt.setp(tax.spines.values(), visible=False)
    tax.tick_params(left=False, labelleft=False)
    tax.tick_params(bottom=False, labelbottom=False)
    # Now we handle positioning, we want the plot to be centered below the data axes, as close to the bottom of the figure as possible
    # without cutting off the xlabel, and retain the original width/height
    tax.set_position([(ax.get_position().width/2)-(tax.get_position().width/2)+ax.get_position().x0,0.05,tax.get_position().width,tax.get_position().height], which="both")
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
    lax.set_position([(ax.get_position().width)-(lax.get_position().width/2)+ax.get_position().x0, .03, lax.get_position().width, lax.get_position().height], which="both")
    # Make sure image is opaque
    fig.set_facecolor("white")
    # Set size to 1080p, resolution of the weather center monitors
    fig.set_size_inches(1920*px, 1080*px)
    # Create a path object to 'productPath' (as defined by the HDWX API), in this case gisproducts/hlma/vhf/ 
    staticProductPath = path.join("products", "hlma", "vhf")
    # Target path for the "static"/non-GIS/transparent image is output/<gisProductPath>/<runPathExtension>/<minute>.png
    staticSavePath = path.join(basePath, "output", staticProductPath, runPathExt, dt.strftime(timeOfPlot, "%M")+".png")
    # Create save directory if it doesn't already exist
    Path(path.dirname(staticSavePath)).mkdir(parents=True, exist_ok=True)
    # Write the image
    fig.savefig(staticSavePath)
    # Write metadata for the product
    writeJson(141, staticProductPath, runPathExt, timeOfPlot)
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
    alreadyPlottedFrames = list()
    # Get path to last hour's json metadata
    lastHourMetadataPath = path.join(basePath, "output", "metadata", "products", "140", dt.strftime(oneHourAgo, "%Y%m%d%H00")+".json")
    # Read in last hour's metadata
    if path.exists(lastHourMetadataPath):
        with open(lastHourMetadataPath, "r") as jsonRead:
            lastHourData = json.load(jsonRead)
        # Add already-generated frames to alreadyPlottedFrames list
        # The valid time gets converted first from an int to a string, then the string is trimmed to only include HHMM
        [alreadyPlottedFrames.append(str(frame["valid"])[-4:]+"00") for frame in lastHourData["productFrames"]]
    # Do the same thing for this hour's metadata
    thisHourMetadataPath = path.join(basePath, "output", "metadata", "products", "140", dt.strftime(now, "%Y%m%d%H00")+".json")
    if path.exists(thisHourMetadataPath):
        with open(thisHourMetadataPath, "r") as jsonRead:
            thisHourData = json.load(jsonRead)
        [alreadyPlottedFrames.append(str(frame["valid"])[-4:]+"00") for frame in thisHourData["productFrames"]]
    # Plot every file in the input directory
    for file in sorted(listdir(inputPath)):
        timeOfFileArr = file.split("_")
        if timeOfFileArr[2] in alreadyPlottedFrames:
            continue
        timeOfFile = dt.strptime("20"+timeOfFileArr[1]+timeOfFileArr[2], "%Y%m%d%H%M%S")
        if timeOfFile < now - timedelta(hours=2):
            remove(path.join(inputPath, file))
            continue
        makeOneMinPlots(path.join(inputPath, file))