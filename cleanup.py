#!/usr/bin/env python3
# Purges no-longer-needed files from hlma plotting
# Created on 19 December 2021 by Sam Gardner <stgardner4@tamu.edu>

from datetime import datetime as dt, timedelta
from os import path, walk, remove
from shutil import rmtree

if __name__ == "__main__":
    now = dt.now()
    basePath = path.dirname(path.abspath(__file__))
    outputPath = path.join(basePath, "output")
    if path.exists(outputPath):
        for root, dirs, files in walk(outputPath):
            for name in files:
                filepath = path.join(basePath, root, name)
                if filepath.endswith(".json"):
                    deleteAfter = timedelta(days=2)
                else:
                    deleteAfter = timedelta(minutes=20)
                createTime = dt.fromtimestamp(path.getmtime(filepath))
                if createTime < now - deleteAfter:
                    remove(filepath)
    if path.exists(path.join(basePath, "radarInput")):
        rmtree(path.join(basePath, "radarInput"))