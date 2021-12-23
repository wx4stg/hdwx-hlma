#!/usr/bin/env python3
# Python-based plotting of Houston Lightning Mapping Array data for next-gen HDWX
# Created 21 December 2021 by Sam Gardner <stgardner4@tamu.edu>

from os import path, listdir
from pyxlma.lmalib.io import read as lma_read
from matplotlib import pyplot as plt
from matplotlib import image as mpimage
from cartopy import crs as ccrs
from cartopy import feature as cfeat
from metpy.plots import USCOUNTIES
import pandas as pd
from datetime import datetime as dt
from pathlib import Path


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

def writeJson():
    print("yay")

def makeLmaPlot(lmaFilePath):
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
    # Extract the time of the plot
    timeOfPlot = lmaData["Datetime"][0].to_pydatetime()
    timeOfPlotStr = dt.strftime(timeOfPlot, "%Y%m%d%H%M")
    extent = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
    path.sep
    pathExt = path.join("hlma", "vhf", dt.strftime(timeOfPlot, "%Y"), dt.strftime(timeOfPlot, "%m"), dt.strftime(timeOfPlot, "%d"), dt.strftime(timeOfPlot, "%H")+"00", dt.strftime(timeOfPlot, "%M")+".png")
    gisSavePath = path.join(basePath, "output", "gisproducts", pathExt)
    Path(path.dirname(gisSavePath)).mkdir(parents=True, exist_ok=True)
    fig.savefig(gisSavePath, transparent=True, bbox_inches=extent)
    ax.add_feature(USCOUNTIES.with_scale("5m"), edgecolor="gray", zorder=2)
    ax.add_feature(cfeat.STATES.with_scale("50m"), linewidth=0.5, zorder=3)
    ax.add_feature(cfeat.COASTLINE.with_scale("50m"), linewidth=0.5, zorder=3)
    ax.set_position([0.05, 0.11, .9, .87])
    tax = fig.add_axes([0,0,(ax.get_position().width/3),.05])
    tax.text(0.5, 0.5, "Houston LMA VHF Sources\nValid "+timeOfPlot.strftime("%-d %b %Y %H%MZ"), horizontalalignment="center", verticalalignment="center", fontsize=16)
    tax.set_xlabel("Python HDWX -- Send bugs to stgardner4@tamu.edu")
    plt.setp(tax.spines.values(), visible=False)
    tax.tick_params(left=False, labelleft=False)
    tax.tick_params(bottom=False, labelbottom=False)
    tax.set_position([(ax.get_position().width/2)-(tax.get_position().width/2)+ax.get_position().x0,0.05,(ax.get_position().width/3),.05], which="both")
    lax = fig.add_axes([0,0,(ax.get_position().width/3),1])
    lax.set_aspect(2821/11071)
    lax.axis("off")
    plt.setp(lax.spines.values(), visible=False)
    atmoLogo = mpimage.imread("assets/atmoLogo.png")
    lax.imshow(atmoLogo)
    lax.set_position([(ax.get_position().width)-(lax.get_position().width/2)+ax.get_position().x0, .03, lax.get_position().width, lax.get_position().height], which="both")
    fig.set_facecolor("white")
    fig.set_size_inches(1920*px, 1080*px)
    staticSavePath = path.join(basePath, "output", "products", pathExt)
    Path(path.dirname(staticSavePath)).mkdir(parents=True, exist_ok=True)
    fig.savefig(staticSavePath)


if __name__ == "__main__":
    # get path to starting dir
    basePath = path.dirname(path.abspath(__file__))
    # get path to input files
    inputPath = path.join(basePath, "lightningin")
    # Plot every file in the input directory
    for file in listdir(inputPath):
        makeLmaPlot(path.join(inputPath, file))