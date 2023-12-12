#!/usr/bin/env python3
# Script to pass HLMA data from thor to the generation server for python-based HDWX
# Created 23 December 2021 by Sam Gardner <stgardner4@tamu.edu>

from datetime import datetime as dt, timedelta
from os import path, listdir, sep, remove
import json
from shutil import copyfile
from pathlib import Path


def getLmaFilesBetweenTimes(startTime, endTime, giveFullPath=False):
    # get base directory
    basePath = path.dirname(path.abspath(__file__))
    inputPath = path.join(basePath, "lightningin")
    filesAvailable = sorted(listdir(inputPath))
    filesToGrab = list()
    for file in reversed(filesAvailable):
        fileDt = dt.strptime(file.split("_")[1][-4:]+file.split("_")[2], "%m%d%H%M%S").replace(year=dt.utcnow().year)
        if fileDt < endTime:
                if giveFullPath:
                    filesToGrab.append(path.join(inputPath, file))
                else:
                    filesToGrab.append(file)
                if fileDt <= startTime:
                    return list(reversed(filesToGrab))


if __name__ == '__main__':
    # Get base directory
    basePath = path.dirname(path.abspath(__file__))
    # Get the current time
    now = dt.utcnow()
    # Get the time one hour ago
    oneHourAgo = now - timedelta(hours=1)
    lastHourJsonPath = path.join(basePath, "output", "metadata", "products", "155", dt.strftime(oneHourAgo, "%Y%m%d%H00")+".json")
    generatedFrames = list()
    if path.exists(lastHourJsonPath):
        with open(lastHourJsonPath, "r") as jsonRead:
            lastHourDict = json.load(jsonRead)
        [generatedFrames.append(frame["valid"]) for frame in lastHourDict["productFrames"]]
    currentHourJsonPath = path.join(basePath, "output", "metadata", "products", "155", dt.strftime(now, "%Y%m%d%H00")+".json")
    if path.exists(currentHourJsonPath):
        with open(currentHourJsonPath, "r") as jsonRead:
            currentHourDict = json.load(jsonRead)
        [generatedFrames.append(frame["valid"]) for frame in currentHourDict["productFrames"]]
    # Generates a list of the hypothetical file paths to this hour and last hour of HLMA data
    lastHourInt = int(dt.strftime(oneHourAgo, "%Y%m%d%H00"))
    currentHourInt = int(dt.strftime(now, "%Y%m%d%H00"))
    currentTimeInt = int(dt.strftime(now, "%Y%m%d%H%M"))
    filesToCopy = list()
    lmaDataBasePath = path.join(path.abspath(sep), "home", "lma_admin", "lma", "realtime", "processed_data", "persec")
    if not path.exists(lmaDataBasePath):
        exit()
    dataInputPath = path.join(basePath, "lightningin")
    Path(dataInputPath).mkdir(parents=True, exist_ok=True)
    for min in range(0, 60):
        if str(lastHourInt + min) not in generatedFrames:
            hypPathLastHour = path.join(lmaDataBasePath, dt.strftime(oneHourAgo, "%Y")[-2:]+dt.strftime(oneHourAgo, "%m%d"), dt.strftime(oneHourAgo ,"%H"), str(f'{min:02}'))
            if path.exists(hypPathLastHour):
                contentsOfLastHour = listdir(hypPathLastHour)
                for file in contentsOfLastHour:
                    if file.endswith(".dat.gz"):
                        copyfile(path.join(hypPathLastHour, file), path.join(dataInputPath, file))
        if str(currentHourInt + min) not in generatedFrames:
            if (currentHourInt + min) < currentTimeInt:
                hypPathCurrentHour = path.join(lmaDataBasePath, dt.strftime(now, "%Y")[-2:]+dt.strftime(now, "%m%d"), dt.strftime(now ,"%H"), str(f'{min:02}'))
                if path.exists(hypPathCurrentHour):
                    contentsOfCurrentHour = listdir(hypPathCurrentHour)
                    for file in contentsOfCurrentHour:
                        if file.endswith(".dat.gz"):
                            copyfile(path.join(hypPathCurrentHour, file), path.join(dataInputPath, file))
    for file in listdir(dataInputPath):
        timeOfFileArr = file.split("_")
        timeOfFile = dt.strptime("20"+timeOfFileArr[1]+timeOfFileArr[2], "%Y%m%d%H%M%S").replace(minute=0)
        if timeOfFile < now - timedelta(hours=2):
            remove(path.join(dataInputPath, file))