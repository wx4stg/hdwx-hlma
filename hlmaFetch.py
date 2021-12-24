#!/usr/bin/env python3
# Script to pass HLMA data from thor to weather-dev for next-gen HDWX
# Created 23 December 2021 by Sam Gardner <stgardner4@tamu.edu>

from datetime import datetime as dt, timedelta
from os import path, listdir, terminal_size
import json
from shutil import copyfile
from pathlib import Path

if __name__ == '__main__':
    # Get base directory
    basePath = path.dirname(path.abspath(__file__))
    # Get the current time
    now = dt.utcnow()
    # Get the time one hour ago
    oneHourAgo = now - timedelta(hours=1)
    lastHourJsonPath = path.join(basePath, "output", "metadata", "products", "100", dt.strftime(oneHourAgo, "%Y%m%d%H00")+".json")
    with open(lastHourJsonPath, "r") as jsonRead:
        lastHourDict = json.load(jsonRead)
    generatedFrames = [frame["valid"] for frame in lastHourDict["productFrames"]]
    currentHourJsonPath = path.join(basePath, "output", "metadata", "products", "100", dt.strftime(now, "%Y%m%d%H00")+".json")
    with open(currentHourJsonPath, "r") as jsonRead:
        currentHourDict = json.load(jsonRead)
    [generatedFrames.append(frame["valid"]) for frame in currentHourDict["productFrames"]]
    # Generates a list of the hypothetical file paths to this hour and last hour of HLMA data
    lastHourInt = int(dt.strftime(oneHourAgo, "%Y%m%d%H00"))
    currentHourInt = int(dt.strftime(now, "%Y%m%d%H00"))
    currentTimeInt = int(dt.strftime(now, "%Y%m%d%H%M"))
    filesToCopy = list()
    lmaDataBasePath = path.join("home", "lma_admin", "lma", "realtime", "processed_data")
    dataInputPath = path.join(basePath, "lightningin")
    Path.mkdir(dataInputPath, parents=True, exist_ok=True)
    for min in range(0, 60):
        if (lastHourInt + min) not in generatedFrames:
            hypPathLastHour = path.join(lmaDataBasePath, dt.strftime(oneHourAgo, "%Y")[-2:]+dt.strftime(oneHourAgo, "%m%d"), dt.strftime(oneHourAgo ,"%H"), dt.strftime(oneHourAgo, "%m"))
            contentsOfLastHour = listdir(hypPathLastHour)
            for file in contentsOfLastHour:
                if file.endswith(".dat.gz"):
                    copyfile(path.join(hypPathLastHour, file), dataInputPath)
        if (currentHourInt + min) not in generatedFrames:
            if (currentHourInt + min) < currentTimeInt:
                hypPathCurrentHour = path.join(lmaDataBasePath, dt.strftime(now, "%Y")[-2:]+dt.strftime(now, "%m%d"), dt.strftime(now ,"%H"), dt.strftime(now, "%m"))
                contentsOfCurrentHour = listdir(hypPathCurrentHour)
                for file in contentsOfCurrentHour:
                    if file.endswith(".dat.gz"):
                        copyfile(path.join(hypPathCurrentHour, file), dataInputPath)
