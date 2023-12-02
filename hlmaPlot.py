#!/usr/bin/env python3
# Python-based plotting of Houston Lightning Mapping Array data for python-based HDWX
# Created 21 December 2021 by Sam Gardner <stgardner4@tamu.edu>

from os import path, listdir, remove
import sys
import pyart
from pyxlma.lmalib.io import read as lma_read
from pyxlma.lmalib.flash.cluster import cluster_flashes
from pyxlma.lmalib.grid import create_regular_grid, assign_regular_bins, events_to_grid
from pyxlma.plot.xlma_plot_feature import color_by_time, plot_points, subset
from pyxlma.plot.xlma_base_plot import subplot_labels, inset_view, BlankPlot
from matplotlib import pyplot as plt
from matplotlib import dates as pltdates
from matplotlib.patheffects import withStroke
from cartopy import crs as ccrs
from cartopy import feature as cfeat
import numpy as np
from metpy.plots import USCOUNTIES
import pandas as pd
from datetime import datetime as dt, timedelta
from pathlib import Path
import json
import warnings
import radarDataFetch
import xarray as xr


axExtent = [-99.5, -91, 26, 33.5]

# hlmaPlot.py <src/flash> <accumulation>
basePath = path.abspath(path.dirname(__file__))
hasHelpers = False
if path.exists(path.join(basePath, "HDWX_helpers.py")):
    import HDWX_helpers
    hasHelpers = True

def set_size(w, h, ax=None):
    if not ax: ax=plt.gca()
    l = ax.figure.subplotpars.left
    r = ax.figure.subplotpars.right
    t = ax.figure.subplotpars.top
    b = ax.figure.subplotpars.bottom
    figw = float(w)/(r-l)
    figh = float(h)/(t-b)
    ax.figure.set_size_inches(figw, figh)

def addMRMSToFig(fig, ax, cbax, taxtext, time, productID, data="ReflectivityAtLowestAltitude"):
    if time.minute % 2 != 0:
        return
    mrmsGribName = radarDataFetch.fetchRadarClosestToTime(time, data)
    if ".grib" in mrmsGribName:
        datasetFilePath = path.join(basePath, "radarInput", mrmsGribName)
        radarDS = xr.open_dataset(datasetFilePath)
        radarDS = radarDS.sel(latitude=slice(axExtent[3], axExtent[2]), longitude=slice(axExtent[0]+360, axExtent[1]+360))
        if data == "ReflectivityAtLowestAltitude":
            radarData = np.ma.masked_array(radarDS.unknown.data, mask=np.where(radarDS.unknown.data > 5, 0, 1))
            cmap = "pyart_ChaseSpectral"
            vmin=-10
            vmax=80
        elif data == "RadarOnly_QPE_01H":
            radarData = np.ma.masked_array(radarDS.unknown.data, mask=np.where(radarDS.unknown.data > 0, 0, 1))/25.4
            cmap = "viridis"
            vmin=0
            vmax=10
            labels = ax.contour(radarDS.longitude, radarDS.latitude, radarData, levels=range(1, 99, 1), cmap="viridis", vmin=0, vmax=10, linewidths=0.5, transform=ccrs.PlateCarree(), zorder=5)
        rdr = ax.pcolormesh(radarDS.longitude, radarDS.latitude, radarData, cmap=cmap, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree(), zorder=5, alpha=0.5)
        if productID == None:
            return rdr
        if cbax is not None:
            cbax.set_position([.01,0.05,(ax.get_position().width/3),.01])
            cbax.xaxis.set_ticks_position("top")
            cbax.tick_params(axis="x", labelsize=9)
            cbax.xaxis.set_label_position("top")
        if taxtext is not None:
            taxtext.set_text(taxtext.get_text().replace("Houston", "MRMS LL Reflectivity + Houston"))
        cbaxRdr = fig.add_axes([.01,0.035,(ax.get_position().width/3),.01])
        fig.colorbar(rdr, cax=cbaxRdr, orientation="horizontal", extend="max")
        cbaxRdr.set_xlabel("Reflectivity (dBZ)", fontsize=9, labelpad=1)
        cbaxRdr.tick_params(axis="x", labelsize=9)
        if productID == 151:
            lightType = "src"
            cbax.set_xlabel("Seconds after "+(time - timedelta(minutes=1)).strftime("%-d %b %Y %H%MZ"))
        elif productID == 153:
            lightType = "flash"
            cbax.set_xlabel("Flash Extent Density (Flashes/km^2/min)", fontsize=9)
        elif productID == 156:
            lightType = "src-analysis"
        productPath = path.join("products", "hlma", "mrms-"+lightType)
        runPathExtension = path.join(time.strftime("%Y"), time.strftime("%m"), time.strftime("%d"), time.strftime("%H")+"00")
        Path(path.join(basePath, "output", productPath, runPathExtension)).mkdir(parents=True, exist_ok=True)
        if hasHelpers:
            HDWX_helpers.saveImage(fig, path.join(basePath, "output", productPath, runPathExtension, time.strftime("%M.png")))
            HDWX_helpers.writeJson(basePath, productID, time, time.strftime("%M.png"), time, ["0,0", "0,0"], 60)
        else:
            fig.savefig(path.join(basePath, "output", productPath, runPathExtension, time.strftime("%M.png")))


def makeFlashPlots(lmaFilePaths):
    # Silence error_bad_lines warning when reading in LMA data
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        # Read in LMA data
        lmaData, startTimeOfPlot = lma_read.dataset(lmaFilePaths)
    numMins = 0
    if "0060.dat" in lmaFilePaths[0]:
        numMins = 1
    else:
        numMins = 10
    elapsedTimeOfData = numMins*len(lmaFilePaths)
    if elapsedTimeOfData == 1:
        gisProductID = 146
        staticProductID = 147
    elif elapsedTimeOfData == 10:
        gisProductID = 148
        staticProductID = 149
    # Get end time for grid creation
    timeOfPlot = startTimeOfPlot + numMins*timedelta(minutes=len(lmaFilePaths))
    dttuple = (np.datetime64(startTimeOfPlot), np.datetime64(timeOfPlot))
    grid_dt = np.asarray(60, dtype='m8[s]')
    grid_t0 = np.asarray(dttuple[0]).astype('datetime64[ns]')
    grid_t1 = np.asarray(dttuple[1]).astype('datetime64[ns]')
    time_range = (grid_t0, grid_t1+grid_dt, grid_dt)
    # We only want events with chi^2 less than 1
    lmaData = lmaData[{"number_of_events":(lmaData.event_chi2 <= 5.0)}]
    lmaDataClustered = cluster_flashes(lmaData)
    lat_range = (axExtent[2], axExtent[3], 0.025)
    lon_range = (axExtent[0], axExtent[1], 0.025)
    alt_range = (0, 18e3, 1.0e3)
    grid_edge_ranges ={
        'grid_latitude_edge':lat_range,
        'grid_longitude_edge':lon_range,
        'grid_altitude_edge':alt_range,
        'grid_time_edge':time_range,
    }
    grid_center_names ={
        'grid_latitude_edge':'grid_latitude',
        'grid_longitude_edge':'grid_longitude',
        'grid_altitude_edge':'grid_altitude',
        'grid_time_edge':'grid_time',
    }
    event_coord_names = {
        'event_latitude':'grid_latitude_edge',
        'event_longitude':'grid_longitude_edge',
        'event_altitude':'grid_altitude_edge',
        'event_time':'grid_time_edge',
    }
    grid_ds = create_regular_grid(grid_edge_ranges, grid_center_names)
    ds_ev = assign_regular_bins(grid_ds, lmaDataClustered, event_coord_names, pixel_id_var="event_pixel_id", append_indices=True)
    grid_spatial_coords=['grid_time', None, 'grid_latitude', 'grid_longitude']
    event_spatial_vars = ('event_altitude', 'event_latitude', 'event_longitude')
    griddedLmaData = events_to_grid(ds_ev, grid_ds, min_points_per_flash=3, pixel_id_var="event_pixel_id", event_spatial_vars=event_spatial_vars, grid_spatial_coords=grid_spatial_coords)
    fig = plt.figure()
    ax = plt.axes(projection=ccrs.epsg(3857))
    griddedLmaData = griddedLmaData.isel(grid_time=0)
    try:
        flashPcm = ax.pcolormesh(griddedLmaData.flash_extent_density.grid_longitude, griddedLmaData.flash_extent_density.grid_latitude, griddedLmaData.flash_extent_density.data, cmap="plasma", vmin=1, vmax=10, transform=ccrs.PlateCarree(), zorder=4)
    except Exception as e:
        if "GEOSContains" in str(e):
            return
        else:
            raise e
    # Plot station locations
    ax.scatter(lmaDataClustered["station_longitude"], lmaDataClustered["station_latitude"], 8, "white", "o", linewidths=.5, edgecolors="black", transform=ccrs.PlateCarree(), zorder=4)
    # Get a handle to a pixel
    px = 1/plt.rcParams["figure.dpi"]
    # Make the GIS plot have a decent resolution. This wont end up being exactly 1920 by 1080 
    # as it's overridden by the next line "set_extent", but should improve the resolution enough
    set_size(1920*px, 1080*px, ax=ax)
    # Now we set the extent to the coordinates we want
    ax.set_extent(axExtent, crs=ccrs.PlateCarree())
    # Create a path oobject to 'runPathExtension', <year>/<month>/<day>/<hour>00/
    runPathExt = path.join(dt.strftime(timeOfPlot, "%Y"), dt.strftime(timeOfPlot, "%m"), dt.strftime(timeOfPlot, "%d"), dt.strftime(timeOfPlot, "%H")+"00")
    if "--no-gis" not in sys.argv:
        # Create a path object to 'productPath' (as defined by the HDWX API), in this case gisproducts/hlma/vhf/ 
        # Use os.path to preserve Windows compatibility
        gisProductPath = path.join("gisproducts", "hlma", "flash-"+str(numMins*len(lmaFilePaths))+"min")
        # Target path for the "GIS"/transparent image is output/<gisProductPath>/<runPathExtension>/<minute>.png
        gisSavePath = path.join(basePath, "output", gisProductPath, runPathExt, dt.strftime(timeOfPlot, "%M")+".png")
        # Create target directory if it doesn't already exist
        Path(path.dirname(gisSavePath)).mkdir(parents=True, exist_ok=True)
        # Get the exact extent of just the axes without the matplotlib auto-generated whitespace
        extent = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
        # save the figure, but trim the whitespace
        # we do this because including the whitespace would make the data not align to the GIS information in the metadata
        if hasHelpers:
            HDWX_helpers.saveImage(fig, gisSavePath, transparent=True, bbox_inches=extent)
            # Write metadata for the product
            HDWX_helpers.writeJson(basePath, gisProductID, timeOfPlot, path.basename(gisSavePath), timeOfPlot, [str(axExtent[2])+","+str(axExtent[0]), str(axExtent[3])+","+str(axExtent[1])], 60)
        else:
            fig.savefig(gisSavePath, transparent=True, bbox_inches=extent)
    # For the "static"/non-GIS/opaque image, add county/state/coastline borders
    ax.add_feature(USCOUNTIES.with_scale("5m"), edgecolor="gray", zorder=2)
    ax.add_feature(cfeat.STATES.with_scale("10m"), linewidth=0.5, zorder=3)
    ax.add_feature(cfeat.COASTLINE.with_scale("10m"), linewidth=0.5, zorder=3)
    # Apply standard HDWX branding to the image
    HDWX_helpers.dressImage(fig, ax, "Houston LMA "+str(numMins*len(lmaFilePaths))+"-minute Flash Extent Density", timeOfPlot, plotHandle=flashPcm, cbticks=range(1, 11), cbextend="max", colorbarLabel="Flash Extent Density (Flashes/km^2/min)")
    title = fig.axes[2].get_children()[0]
    cbax = fig.axes[1]
    # Create a path object to 'productPath' (as defined by the HDWX API), in this case gisproducts/hlma/vhf/ 
    staticProductPath = path.join("products", "hlma", "flash-"+str(numMins*len(lmaFilePaths))+"min")
    # Target path for the "static"/non-GIS/transparent image is output/<gisProductPath>/<runPathExtension>/<minute>.png
    staticSavePath = path.join(basePath, "output", staticProductPath, runPathExt, dt.strftime(timeOfPlot, "%M")+".png")
    # Create save directory if it doesn't already exist
    Path(path.dirname(staticSavePath)).mkdir(parents=True, exist_ok=True)
    # Write the image
    if hasHelpers:
        HDWX_helpers.saveImage(fig, staticSavePath)
        # Write metadata for the product
        HDWX_helpers.writeJson(basePath, staticProductID, timeOfPlot, path.basename(staticSavePath), timeOfPlot, ["0,0", "0,0"], 60)
    else:
        fig.savefig(staticSavePath)
    if len(lmaFilePaths) == 1:
        # Now we can add MRMS to the figure
        addMRMSToFig(fig, ax, cbax, title, timeOfPlot, 153)
    # Close figure when done (memory management)
    plt.close(fig)
    # Make a two panel lightning/radar plot
    if len(lmaFilePaths) == 10:
        twoPanelFig = plt.figure()
        twoPanelFig.set_size_inches(1920*px, 1080*px)
        twoPanelAx1 = twoPanelFig.add_axes([0.025, 0.2, 0.35, 0.9], projection=ccrs.epsg(3857))
        try:
            flashPcm = twoPanelAx1.pcolormesh(griddedLmaData.flash_extent_density.grid_longitude, griddedLmaData.flash_extent_density.grid_latitude, griddedLmaData.flash_extent_density.data, cmap="plasma", vmin=1, vmax=10, transform=ccrs.PlateCarree(), zorder=4)
        except Exception as e:
            if "GEOSContains" in str(e):
                return
            else:
                raise e
        twoPanelAx1.set_extent(axExtent, crs=ccrs.PlateCarree())
        if timeOfPlot.minute % 2 == 0:
            timeForRadar = timeOfPlot
        else:
            timeForRadar = timeOfPlot - timedelta(minutes=1)
        rdr1 = addMRMSToFig(twoPanelFig, twoPanelAx1, None, None, timeForRadar, None, "ReflectivityAtLowestAltitude")
        cbax11 = twoPanelFig.add_axes([0.425, twoPanelAx1.get_position().y0, 0.01, twoPanelAx1.get_position().height])
        twoPanelFig.colorbar(rdr1, cax=cbax11, orientation="vertical", label="Reflectivity (dBZ)", extend="max")
        cbax11.yaxis.set_ticks_position("left")
        cbax11.yaxis.set_label_position("left")
        cbax12 = twoPanelFig.add_axes([0.45, twoPanelAx1.get_position().y0, 0.01, twoPanelAx1.get_position().height])
        twoPanelFig.colorbar(flashPcm, cax=cbax12, orientation="vertical", label="Flash Extent Density (flash/km^2/min)", extend="max")
        twoPanelAx2 = twoPanelFig.add_axes([0.5, .2, 0.35, 0.9], projection=ccrs.epsg(3857))
        twoPanelAx2.set_extent(axExtent, crs=ccrs.PlateCarree())
        latRange = [axExtent[2], axExtent[3]]
        lonRange = [axExtent[0], axExtent[1]]
        altRange = [0, 21]
        lonSet, latSet, altSet, timeSet, selectedData = subset(lmaData.event_longitude.values, lmaData.event_latitude.values, lmaData.event_altitude.values/1000, pd.Series(lmaData.event_time), lmaData.event_chi2.values, lmaData.event_stations.values, lonRange, latRange, altRange, [startTimeOfPlot, timeOfPlot], 5.0, 6.0)
        vmin, vmax, relcolors = color_by_time(timeSet, [startTimeOfPlot, timeOfPlot])
        rdr2 = addMRMSToFig(twoPanelFig, twoPanelAx2, None, None, timeForRadar, None, "RadarOnly_QPE_01H")
        cbax21 = twoPanelFig.add_axes([0.9, twoPanelAx1.get_position().y0, 0.01, twoPanelAx1.get_position().height])
        twoPanelFig.colorbar(rdr2, cax=cbax21, orientation="vertical", label="1-hr radar estimated rainfall (in.)", extend="max")
        cbax21.yaxis.set_ticks_position("left")
        cbax21.yaxis.set_label_position("left")
        cbax22 = twoPanelFig.add_axes([0.925, twoPanelAx1.get_position().y0, 0.01, twoPanelAx1.get_position().height])
        # if len(relcolors) > 0:
        vhfSct = twoPanelAx2.scatter(lonSet, latSet, 1, relcolors, ",", transform=ccrs.PlateCarree(), zorder=4, cmap="rainbow", vmin=vmin, vmax=vmax)
        twoPanelFig.colorbar(vhfSct, cax=cbax22, label=f"Seconds after {startTimeOfPlot.strftime('%-d %b %Y %H%MZ')}")

        twoPanelAx1.scatter(lmaDataClustered["station_longitude"], lmaDataClustered["station_latitude"], 8, "white", "o", linewidths=.5, edgecolors="black", transform=ccrs.PlateCarree(), zorder=4)
        twoPanelAx2.scatter(lmaDataClustered["station_longitude"], lmaDataClustered["station_latitude"], 8, "white", "o", linewidths=.5, edgecolors="black", transform=ccrs.PlateCarree(), zorder=4)

        twoPanelAx1.add_feature(USCOUNTIES.with_scale("5m"), edgecolor="gray", zorder=2)
        twoPanelAx1.add_feature(cfeat.STATES.with_scale("10m"), linewidth=0.5, zorder=3)
        twoPanelAx1.add_feature(cfeat.COASTLINE.with_scale("10m"), linewidth=0.5, zorder=3)

        twoPanelAx2.add_feature(USCOUNTIES.with_scale("5m"), edgecolor="gray", zorder=2)
        twoPanelAx2.add_feature(cfeat.STATES.with_scale("10m"), linewidth=0.5, zorder=3)
        twoPanelAx2.add_feature(cfeat.COASTLINE.with_scale("10m"), linewidth=0.5, zorder=3)

        twoPanelAx1.set_title("Reflectivity + Flash Extent Density")
        twoPanelAx2.set_title("1-hr QPE + VHF Sources")

        twoPanelAx3 = twoPanelFig.add_axes([twoPanelAx1.get_position().x0, 0.14, cbax22.get_position().x1 - twoPanelAx1.get_position().x0, 0.15])
        timeDensity = lmaData.event_time.groupby(lmaData.event_time.dt.minute).count()
        firstFFDate = pd.Timestamp(lmaData.event_time.data[-1]).to_pydatetime().replace(second=0, microsecond=0)
        datesForFFPlot = [firstFFDate + timedelta(minutes=i+1, seconds=30) for i in range(len(timeDensity.minute))]
        twoPanelAx3.plot(datesForFFPlot, timeDensity.data, color="black", linewidth=1)
        twoPanelAx3.scatter(datesForFFPlot, timeDensity.data, color="black", s=5)
        twoPanelAx3.text(datesForFFPlot[-1], 0.01, f"{int(timeDensity.data[-1])} src", ha="center", va="bottom", path_effects=[withStroke(linewidth=3, foreground="white")], transform=twoPanelAx3.get_xaxis_transform())
        twoPanelAx3.xaxis.set_major_formatter(pltdates.DateFormatter("%H:%MZ"))
        twoPanelAx3.set_xlabel("Time")
        twoPanelAx3.set_ylabel("Number of Sources")
        twoPanelAx3.yaxis.set_label_position("right")
        lax = twoPanelFig.add_axes([0.75, 0.01, 0.25, 0.15])
        if hasHelpers:
            HDWX_helpers.dressImage(twoPanelFig, None, "HLMA Flash Flood Analysis", timeOfPlot, lax=lax)
        lax.set_position([0.74, 0.01, lax.get_position().width, lax.get_position().height])

        lolItsThreePanelsNow = path.join(basePath, "output", "products", "hlma", "ffanalysis", runPathExt, dt.strftime(timeOfPlot, "%M")+".png")
        Path(path.dirname(lolItsThreePanelsNow)).mkdir(parents=True, exist_ok=True)
        if hasHelpers:
            HDWX_helpers.saveImage(twoPanelFig, lolItsThreePanelsNow)
            HDWX_helpers.writeJson(basePath, 157, timeOfPlot, dt.strftime(timeOfPlot, "%M.png"), timeOfPlot, ["0,0", "0,0"], 60)
        else:
            twoPanelFig.savefig(lolItsThreePanelsNow)


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
        lmaPlotID = 154
    else:
        gisProductID = 143
        staticProductID = 144
        lmaPlotID = 155
    # Add the number of one minute files ingested to the start time to get the end time
    timeOfPlot = startTimeOfPlot + timedelta(minutes=len(lmaFilePaths))
    # Create fig/ax
    fig = plt.figure()
    ax = plt.axes(projection=ccrs.epsg(3857))
    # We want to "mask out" points where the event chi^2 is greater than 1 and points outside our lat/lon/alt range
    latRange = [axExtent[2], axExtent[3]]
    lonRange = [axExtent[0], axExtent[1]]
    altRange = [0, 21]
    # Subset the data to our parameters (again... I was unsuccessful when attempting to get pyxlma to read in the masked arrays from earlier)
    lonSet, latSet, altSet, timeSet, selectedData = subset(lmaData.event_longitude.values, lmaData.event_latitude.values, lmaData.event_altitude.values/1000, pd.Series(lmaData.event_time), lmaData.event_chi2.values, lmaData.event_stations.values, lonRange, latRange, altRange, [startTimeOfPlot, timeOfPlot], 5.0, 6.0)
    vmin, vmax, relcolors = color_by_time(timeSet, [startTimeOfPlot, timeOfPlot])
    if vmin < 1:
        vmin = 0
    if vmax > (len(lmaFilePaths)*60 - 1):
        vmax = len(lmaFilePaths)*60
    # Plot data
    vhfSct = ax.scatter(lonSet, latSet, 1, relcolors, ",", transform=ccrs.PlateCarree(), zorder=4, cmap="rainbow", vmin=vmin, vmax=vmax)
    # Plot station locations
    ax.scatter(lmaData["station_longitude"], lmaData["station_latitude"], 8, "white", "o", linewidths=.5, edgecolors="black", transform=ccrs.PlateCarree(), zorder=4)
    # Get a handle to a pixel
    px = 1/plt.rcParams["figure.dpi"]
    # Make the GIS plot have a decent resolution. This wont end up being exactly 1920 by 1080 
    # as it's overridden by the next line "set_extent", but should improve the resolution enough
    set_size(1920*px, 1080*px, ax=ax)
    # Now we set the extent to the coordinates we want
    # Create a path oobject to 'runPathExtension', <year>/<month>/<day>/<hour>00/
    ax.set_extent(axExtent, crs=ccrs.PlateCarree())
    runPathExt = path.join(dt.strftime(timeOfPlot, "%Y"), dt.strftime(timeOfPlot, "%m"), dt.strftime(timeOfPlot, "%d"), dt.strftime(timeOfPlot, "%H")+"00")
    if "--no-gis" not in sys.argv:
        # Create a path object to 'productPath' (as defined by the HDWX API), in this case gisproducts/hlma/vhf/ 
        # Use os.path to preserve Windows compatibility
        gisProductPath = path.join("gisproducts", "hlma", "vhf-"+str(len(lmaFilePaths))+"min")
        # Target path for the "GIS"/transparent image is output/<gisProductPath>/<runPathExtension>/<minute>.png
        gisSavePath = path.join(basePath, "output", gisProductPath, runPathExt, dt.strftime(timeOfPlot, "%M")+".png")
        # Create target directory if it doesn't already exist
        Path(path.dirname(gisSavePath)).mkdir(parents=True, exist_ok=True)
        # Get the exact extent of just the axes without the matplotlib auto-generated whitespace
        extent = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
        # save the figure, but trim the whitespace
        # we do this because including the whitespace would make the data not align to the GIS information in the metadata
        # Write metadata for the product
        if hasHelpers:
            HDWX_helpers.saveImage(fig, gisSavePath, transparent=True, bbox_inches=extent)
            HDWX_helpers.writeJson(basePath, gisProductID, timeOfPlot, path.basename(gisSavePath), timeOfPlot, [str(axExtent[2])+","+str(axExtent[0]), str(axExtent[3])+","+str(axExtent[1])], 60)
        else:
            fig.savefig(gisSavePath, transparent=True, bbox_inches=extent)
    # For the "static"/non-GIS/opaque image, add county/state/coastline borders
    ax.add_feature(USCOUNTIES.with_scale("5m"), edgecolor="gray", zorder=2)
    ax.add_feature(cfeat.STATES.with_scale("10m"), linewidth=0.5, zorder=3)
    ax.add_feature(cfeat.COASTLINE.with_scale("10m"), linewidth=0.5, zorder=3)
    # Apply standard HDWX branding to the image
    if len(relcolors) > 0:
        HDWX_helpers.dressImage(fig, ax, "Houston LMA "+str(len(lmaFilePaths))+"-minute VHF Sources", timeOfPlot, plotHandle=vhfSct, cbextend="neither", colorbarLabel="Seconds after "+startTimeOfPlot.strftime("%-d %b %Y %H%MZ"))
    else:
        HDWX_helpers.dressImage(fig, ax, "Houston LMA "+str(len(lmaFilePaths))+"-minute VHF Sources", timeOfPlot, plotHandle=None, cbextend="neither", colorbarLabel="No VHF Sources")
    title = fig.axes[2].get_children()[0]
    cbax = fig.axes[1]
    # Create a path object to 'productPath' (as defined by the HDWX API), in this case gisproducts/hlma/vhf/ 
    staticProductPath = path.join("products", "hlma", "vhf-"+str(len(lmaFilePaths))+"min")
    # Target path for the "static"/non-GIS/transparent image is output/<gisProductPath>/<runPathExtension>/<minute>.png
    staticSavePath = path.join(basePath, "output", staticProductPath, runPathExt, dt.strftime(timeOfPlot, "%M")+".png")
    # Create save directory if it doesn't already exist
    Path(path.dirname(staticSavePath)).mkdir(parents=True, exist_ok=True)
    # Write the image
    if hasHelpers:
        HDWX_helpers.saveImage(fig, staticSavePath)
        # Write metadata for the product
        HDWX_helpers.writeJson(basePath, staticProductID, timeOfPlot, path.basename(staticSavePath), timeOfPlot, ["0,0", "0,0"], 60)
    else:
        fig.savefig(staticSavePath)
    if len(lmaFilePaths) == 1:
        # Now we can add MRMS to the figure
        addMRMSToFig(fig, ax, cbax, title, timeOfPlot, 151)
    # Close figure when done (memory management)
    plt.close(fig)
    # The LMA community has come up with a pretty cool plot of charge density vs lat, lon, time, and altitude, lets make one!
    lmaPlot = BlankPlot(startTimeOfPlot, bkgmap=True, xlim=lonRange, ylim=latRange, zlim=altRange, tlim=[startTimeOfPlot, timeOfPlot], title="Houston LMA "+str(len(lmaFilePaths))+"-minute VHF Sources\nValid "+startTimeOfPlot.strftime("%-d %b %Y %H%MZ")+" through "+timeOfPlot.strftime("%H%MZ"))
    # Plot station locations
    lmaPlot.ax_plan.scatter(lmaData["station_longitude"], lmaData["station_latitude"], 8, "white", "o", linewidths=.5, edgecolors="black", transform=ccrs.PlateCarree(), zorder=4)
    lmaPlotFig = plt.gcf()
    # Add our data
    plot_points(lmaPlot, lonSet, latSet, altSet, timeSet, "rainbow", 5, vmin, vmax, relcolors, edge_color="black", edge_width=0.25)
    # Create save directory if it doesn't already exist
    lmaProductPath = path.join("products", "hlma", "vhf-"+str(len(lmaFilePaths))+"min-analysis")
    lmaSavePath = path.join(basePath, "output", lmaProductPath, runPathExt, dt.strftime(timeOfPlot, "%M")+".png")
    Path(path.dirname(lmaSavePath)).mkdir(parents=True, exist_ok=True)
    # Write the image
    if hasHelpers:
        HDWX_helpers.saveImage(lmaPlotFig, lmaSavePath)
        # Write metadata for the product
        HDWX_helpers.writeJson(basePath, lmaPlotID, timeOfPlot, path.basename(lmaSavePath), timeOfPlot, ["0,0", "0,0"], 60)
    else:
        lmaPlotFig.savefig(lmaSavePath)
    if len(lmaFilePaths) == 1:
        addMRMSToFig(lmaPlotFig, lmaPlot.ax_plan, None, None, timeOfPlot, 156)

def getAlreadyPlottedFrames(productID):
    # Create (empty, but add to it in a sec...) list representing already plotted frames
    alreadyPlottedFrames = list()
    # Get path to last hour's json metadata
    lastHourMetadataPath = path.join(basePath, "output", "metadata", "products", str(productID), dt.strftime(oneHourAgo, "%Y%m%d%H00")+".json")
    # Read in last hour's metadata
    if path.exists(lastHourMetadataPath):
        with open(lastHourMetadataPath, "r") as jsonRead:
            lastHourData = json.load(jsonRead)
        # Add already-generated frames to alreadyPlottedFrames list
        # The valid time gets converted first from an int to a string, then the string is trimmed to only include HHMM
        [alreadyPlottedFrames.append(str(frame["valid"])[-4:]+"00") for frame in lastHourData["productFrames"]]
    # Do the same thing for this hour's metadata
    thisHourMetadataPath = path.join(basePath, "output", "metadata", "products", str(productID), dt.strftime(now, "%Y%m%d%H00")+".json")
    if path.exists(thisHourMetadataPath):
        with open(thisHourMetadataPath, "r") as jsonRead:
            thisHourOneMinData = json.load(jsonRead)
        [alreadyPlottedFrames.append(str(frame["valid"])[-4:]+"00") for frame in thisHourOneMinData["productFrames"]]
    return alreadyPlottedFrames

if __name__ == "__main__":
    # read "src" or "flash" from command line. If neither are provided, we'll assume both source and flash plots are desired
    shouldPlotSrc = True
    shouldPlotFlash = True
    shouldPlot1min = True
    shouldPlot10min = True
    if len(sys.argv) > 1:
        if sys.argv[1] == "src":
            shouldPlotFlash = False
        elif sys.argv[1] == "flash":
            shouldPlotSrc = False
    if len(sys.argv) > 2:
        if sys.argv[2] == "1":
            shouldPlot10min = False
        elif sys.argv[2] == "10":
            shouldPlot1min = False
    # get path to input files
    inputPath = path.join(basePath, "lightningin")
    # Get current time
    now = dt.utcnow()
    # Get time one hour ago
    oneHourAgo = now - timedelta(hours=1)
    if shouldPlot1min:
        if shouldPlotSrc:
            # We're plotting 1 minute VHF Sources!
            # Get the frames that have already been plotted
            alreadyPlottedFrames = getAlreadyPlottedFrames(141)
            # Plot every file in the input directory
            inputDirContents = sorted(listdir(inputPath), reverse=True)
            for file in inputDirContents:
                timeOfFileArr = file.split("_")
                # The time in the filename is the *start*, but the time in the json is the end, so add one minute to the filename time
                timeOfFile = dt.strptime("20"+timeOfFileArr[1]+timeOfFileArr[2], "%Y%m%d%H%M%S") + timedelta(minutes=1)
                if timeOfFile < now - timedelta(hours=2):
                    remove(path.join(inputPath, file))
                    continue
                if timeOfFile.strftime("%H%M%S") in alreadyPlottedFrames:
                    continue
                if timeOfFile < oneHourAgo:
                    continue
                makeSourcePlots([path.join(inputPath, file)])
        if shouldPlotFlash:
            # We're plotting 1 minute flash extent density!
            # Get the frames that have already been plotted
            alreadyPlottedFrames = getAlreadyPlottedFrames(147)
            # Plot every file in the input directory
            inputDirContents = sorted(listdir(inputPath), reverse=True)
            for file in inputDirContents:
                timeOfFileArr = file.split("_")
                # The time in the filename is the *start*, but the time in the json is the end, so add one minute to the filename time
                timeOfFile = dt.strptime("20"+timeOfFileArr[1]+timeOfFileArr[2], "%Y%m%d%H%M%S") + timedelta(minutes=1)
                if timeOfFile < now - timedelta(hours=2):
                    remove(path.join(inputPath, file))
                    continue
                if timeOfFile.strftime("%H%M%S") in alreadyPlottedFrames:
                    continue
                if timeOfFile < oneHourAgo:
                    continue
                makeFlashPlots([path.join(inputPath, file)])
    if shouldPlot10min:
        if shouldPlotSrc:
            # We're plotting 10 minute VHF Sources!
            # Get the frames that have already been plotted
            alreadyPlottedFrames = getAlreadyPlottedFrames(144)
            inputDirContents = sorted(listdir(inputPath), reverse=True)
            for i in range(0, len(inputDirContents)-10):
                lastFileInRange = inputDirContents[i]
                timeOfLastFileArr = lastFileInRange.split("_")
                timeOfLastFile = dt.strptime("20"+timeOfLastFileArr[1]+timeOfLastFileArr[2], "%Y%m%d%H%M%S") + timedelta(minutes=1)
                if timeOfLastFile.strftime("%H%M%S") in alreadyPlottedFrames:
                    continue
                if timeOfLastFile < oneHourAgo:
                    continue
                filesToPlot = [path.join(inputPath, fileToInclude) for fileToInclude in inputDirContents[i:i+10]]
                makeSourcePlots(filesToPlot)
        if shouldPlotFlash:
            # We're plotting 10 minute flash extent density!
            # Get the frames that have already been plotted
            alreadyPlottedFrames = getAlreadyPlottedFrames(149)
            inputDirContents = sorted(listdir(inputPath), reverse=True)
            for i in range(0, len(inputDirContents)-10):
                lastFileInRange = inputDirContents[i]
                timeOfLastFileArr = lastFileInRange.split("_")
                timeOfLastFile = dt.strptime("20"+timeOfLastFileArr[1]+timeOfLastFileArr[2], "%Y%m%d%H%M%S") + timedelta(minutes=1)
                if timeOfLastFile.strftime("%H%M%S") in alreadyPlottedFrames:
                    continue
                if timeOfLastFile < oneHourAgo:
                    continue
                filesToPlot = [path.join(inputPath, fileToInclude) for fileToInclude in inputDirContents[i:i+10]]
                makeFlashPlots(filesToPlot)
