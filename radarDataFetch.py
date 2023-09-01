#!/usr/bin/env python3
# Radar data fetching for HLMA "crossover" products
# Created 1 May 2022 by Sam Gardner <stgardner4@tamu.edu>

from os import path, remove
from pathlib import Path
import pandas as pd
from datetime import datetime as dt, timedelta
import requests
import gzip
import shutil
import urllib


basePath = path.realpath(path.dirname(__file__))

def downloadFile(fileName, data):
    output = path.join(basePath, "radarInput", fileName)
    if path.exists(path.join(basePath, "radarInput", fileName.replace(".gz", ""))):
        return output.replace(".gz", "")
    print("Downloading "+fileName)
    urlToFetch = f"https://mrms.ncep.noaa.gov/data/2D/{data}/{fileName}"
    mrmsData = requests.get(urlToFetch)
    if mrmsData.status_code == 200:
        with open(output, "wb") as fileWrite:
            fileWrite.write(mrmsData.content)
        with gzip.open(output, "rb") as f_in:
            with open(output.replace(".gz", ""), "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        remove(output)
        return output.replace(".gz", "")

def fetchRadarClosestToTime(time, data):
    Path(path.join(basePath, "radarInput")).mkdir(parents=True, exist_ok=True)
    try:
        gribList = pd.read_html(f"https://mrms.ncep.noaa.gov/data/2D/{data}/")[0].dropna(how="any")
    except urllib.error.URLError as e:
        import subprocess
        from io import BytesIO
        gribListCurlProc = subprocess.run(["curl", f"https://mrms.ncep.noaa.gov/data/2D/{data}/"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        gribList = pd.read_html(BytesIO(gribListCurlProc.stdout.encode("utf-8")))[0].dropna(how="any")
    gribList = gribList[~gribList.Name.str.contains("latest") == True].reset_index()
    if data == "ReflectivityAtLowestAltitude":
        gribList["pyDateTimes"] = [dt.strptime(filename, "MRMS_ReflectivityAtLowestAltitude_00.50_%Y%m%d-%H%M%S.grib2.gz") for filename in gribList["Name"]]
    elif data == "RadarOnly_QPE_01H":
        gribList["pyDateTimes"] = [dt.strptime(filename, "MRMS_RadarOnly_QPE_01H_00.00_%Y%m%d-%H%M%S.grib2.gz") for filename in gribList["Name"]]
    gribList = gribList.set_index(["pyDateTimes"])
    for moasicTime in reversed(gribList.index):
        if moasicTime < time:
            if moasicTime >= (time - timedelta(minutes=2)):
                return downloadFile(gribList[gribList.index == moasicTime]["Name"][0], data)
