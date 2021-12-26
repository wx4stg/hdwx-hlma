#!/usr/bin/env python3
# Script for generating an icon of every color (used for generating the GR2A lightning icons)
# Created 25 December 2021 by Sam Gardner <stgardner4@tamu.edu>

from PIL import Image
import numpy as np
from os import path

def make_icons():
    # Replace the color RGB = [50, 100, 150] in template with target hue
    templateImage = Image.open("assets/icontemplate.png")
    templateImage = templateImage.convert("RGBA")
    origData = np.array(templateImage)
    r, g, b, a = origData.T
    replaceMask = (r == 50) & (g == 100) & (b == 150)
    sum = 0
    for a in range(0, 255):
        # R = 255, G = 0->255, B = 0
        sum += 1
        print(sum)
        targetR = 255
        targetG = a
        targetB = 0
        modData = origData.copy()
        modData[..., :-1][replaceMask.T] = (targetR, targetG, targetB)
        imgToAdd = Image.fromarray(modData)
        if path.exists("assets/icons.png"):
            currentImage = Image.open("assets/icons.png")
            newImage = Image.new("RGBA", (currentImage.size[0]+imgToAdd.size[0], currentImage.size[1]))
            newImage.paste(currentImage, (0, 0))
        else:
            newImage = Image.new("RGBA", imgToAdd.size)
        newImage.paste(imgToAdd, (newImage.size[0]-imgToAdd.size[0], 0))
        newImage.save("assets/icons.png")
    for b in range(0, 255):
        # R = 255->0, G = 255, B = 0
        sum += 1
        print(sum)
        targetR = 255 - b
        targetG = 255
        targetB = 0
        modData = origData.copy()
        modData[..., :-1][replaceMask.T] = (targetR, targetG, targetB)
        imgToAdd = Image.fromarray(modData)
        if path.exists("assets/icons.png"):
            currentImage = Image.open("assets/icons.png")
            newImage = Image.new("RGBA", (currentImage.size[0]+imgToAdd.size[0], currentImage.size[1]))
            newImage.paste(currentImage, (0, 0))
        else:
            newImage = Image.new("RGBA", imgToAdd.size)
        newImage.paste(imgToAdd, (newImage.size[0]-imgToAdd.size[0], 0))
        newImage.save("assets/icons.png")

    for c in range(0, 255):
        # R = 0, G = 255, B = 0->255
        sum += 1
        print(sum)
        targetR = 0
        targetG = 255
        targetB = c
        modData = origData.copy()
        modData[..., :-1][replaceMask.T] = (targetR, targetG, targetB)
        imgToAdd = Image.fromarray(modData)
        if path.exists("assets/icons.png"):
            currentImage = Image.open("assets/icons.png")
            newImage = Image.new("RGBA", (currentImage.size[0]+imgToAdd.size[0], currentImage.size[1]))
            newImage.paste(currentImage, (0, 0))
        else:
            newImage = Image.new("RGBA", imgToAdd.size)
        newImage.paste(imgToAdd, (newImage.size[0]-imgToAdd.size[0], 0))
        newImage.save("assets/icons.png")

    for d in range(0, 255):
        # R = 0, G = 255->0, B = 255
        sum += 1
        print(sum)
        targetR = 0
        targetG = 255 - d
        targetB = 255
        modData = origData.copy()
        modData[..., :-1][replaceMask.T] = (targetR, targetG, targetB)
        imgToAdd = Image.fromarray(modData)
        if path.exists("assets/icons.png"):
            currentImage = Image.open("assets/icons.png")
            newImage = Image.new("RGBA", (currentImage.size[0]+imgToAdd.size[0], currentImage.size[1]))
            newImage.paste(currentImage, (0, 0))
        else:
            newImage = Image.new("RGBA", imgToAdd.size)
        newImage.paste(imgToAdd, (newImage.size[0]-imgToAdd.size[0], 0))
        newImage.save("assets/icons.png")

    for e in range(0, 256):
        # R = 0, G = 255->0, B = 255
        sum += 1
        print(sum)
        targetR = e
        targetG = 0
        targetB = 255
        modData = origData.copy()
        modData[..., :-1][replaceMask.T] = (targetR, targetG, targetB)
        imgToAdd = Image.fromarray(modData)
        if path.exists("assets/icons.png"):
            currentImage = Image.open("assets/icons.png")
            newImage = Image.new("RGBA", (currentImage.size[0]+imgToAdd.size[0], currentImage.size[1]))
            newImage.paste(currentImage, (0, 0))
        else:
            newImage = Image.new("RGBA", imgToAdd.size)
        newImage.paste(imgToAdd, (newImage.size[0]-imgToAdd.size[0], 0))
        newImage.save("assets/icons.png")

if __name__ == '__main__':
    make_icons()