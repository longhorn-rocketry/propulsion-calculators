#!/usr/bin/python3
"""
FUEL GRAIN REGRESSION CALCULATOR v1.0
Written by Daniel Teal in 2019 for the LRA

This script can approximate the regression of the inside of a hybrid rocket engine fuel grain during combustion.

In general, the script takes a black-and-white image of the fuel grain cross-section and a regression rate, then, at a given number of timesteps afterward, calculates the regressed cross-section, area, and perimeter. Due to quantization effects from individual pixels, this is not perfectly accurate, but all numbers should be within 10% or so.

USAGE:

First, create a black and white PNG image of the cross-section of the fuel grain. The inside of the grain, where the gases will combust, should be white, and everything else in the image should be black (or else the program will try to regress whatever else it is that's not black, too). That is, the image should be, say, a white circle on a black square.

Next, create a FuelGrainState object as follows:

    my_fuelgrain_variable = FuelGrainState(filename='my_image.png', scale=0.01, rate=0.5)

where the filename is the location of the file (either an absolute pathname like 'C:/Users/Username/Documents/image.png' or a pathname relative to the location of this script when it is run, like 'image.png' when the image is in the same folder as this script (run in a terminal currently at that directory)), the scale is the number of inches per image pixel (something like 0.01 in/px (a 4"x4" fuel grain cross-section would have a 400x400px image) or smaller works well), and the rate is the fuel grain regression rate in in/s.

Now, run the calculation with the calculate_regression function as follows:

    calculate_regression(my_fuelgrain_variable, timestep=0.2, num_steps=20, plot=True)

This will calculate the regression in discrete timesteps into the future. The 'timestep' is the time difference between each point, and the 'num_steps' is the number of steps calculated, so, e.g., timestep=0.2s and num=steps=20 simulates the fuel rate regression over 0.2*20=4s. If 'plot' is true, then this will display GRAPHS and pretty colors in a new window. Finally, the function returns a list of lists: [times, areas, perimeters], where 'times' is a list of times in seconds at each timestep (e.g., [0, 0.2, 0.4, ...]), 'areas' is a list of areas in in^2 at each timestep, and perimeters is a list of the perimeter in inches of the fuel grain at each timestep. This data might be useful for plugging this script into another simulation or making a table or something.

WARNING: there is a sweet spot for the pixel resolution of the image and timestep size for a given regression rate. If the timestep is way too big, then any shape will appear to immediately wear away into a circle neglecting most features of the inner area shape. On the other hand, if the timestep is too slow and the image isn't high enough resolution, then each timestep may only regress the pattern by a single pixel (or none!). This is bad because the actual regression rate might have been, say, 1.5px, but the quantization made the algorithm round off with 50% error. In general, aim to have each timestep regress the fuel grain image by 3-10 px.
"""

import math
import sys
import numpy as np
import matplotlib.pyplot as plt
import skimage.filters as flt
import skimage.measure as ms
import subprocess
import time

plt.style.use('dark_background')

class FuelGrainState:
    def __init__(self, image=None, filename='', scale=1, rate=1, time=0):
        if not filename == '': # load image from file if given ...
            self.image = plt.imread(filename)
            if self.image.shape[2] == 4: # convert from RGBA to RGB
                self.image = self.image[:,:,:2]
            if len(self.image.shape) == 3: # convert to grayscale
                self.image = np.sum(self.image, axis=2) / 3
            # convert to binary
            removed_locations = self.image > 0.01
            self.image[removed_locations] = 1
            self.image[removed_locations == False] = 0
        else: # ... or get image directly
            self.image = image
        self.scale = scale
        self.rate = rate
        self.time = time
    def calculate_area(self):
        return np.count_nonzero(self.image) * (self.scale**2) # convert to in^2
    def calculate_perimeter(self):
        return ms.perimeter(self.image) * (self.scale) # convert to in
    def regress(self, timestep):
        blurred_image = flt.gaussian(self.image, sigma=self.rate*timestep/self.scale)
        removed_locations = blurred_image > 0.13 # magic number that depends on gaussian blur distribution to get correct thickness removed
        blurred_image[removed_locations] = 1
        blurred_image[removed_locations == False] = 0
        blurred_image = np.maximum(blurred_image, self.image) # make sure fuel isn't added (can happen for some blur/threshold combinations)
        return FuelGrainState(image=blurred_image, scale=self.scale, rate=self.rate, time=self.time+timestep)

def calculate_regression(initial_state, timestep=1, num_steps=5, plot=True):
    prev_state = initial_state
    times = [0]
    areas = [initial_state.calculate_area()]
    perimeters = [initial_state.calculate_perimeter()]
    stacked_image = np.array(initial_state.image)
    for _ in range(num_steps):
        new_state = prev_state.regress(timestep)
        times.append(times[-1] + timestep)
        areas.append(new_state.calculate_area())
        perimeters.append(new_state.calculate_perimeter())
        stacked_image = stacked_image + new_state.image
        prev_state = new_state
    if plot:
        graphics_axes = plt.subplot(1,2,1)
        width = stacked_image.shape[0] * initial_state.scale
        height = stacked_image.shape[1] * initial_state.scale
        plt.imshow(stacked_image, extent=(-width/2, width/2, -height/2, height/2))
        plt.xlabel('Distance (in)')
        plt.ylabel('Distance (in)')
        plt.title('Fuel Grain Regression')
        perimeter_axes = plt.subplot(2,2,2)
        perimeter_axes.plot(times, perimeters)
        xmin, xmax, ymin, ymax = plt.axis()
        plt.ylim(0, ymax) # make graph start at 0
        plt.xlabel('Time (s)')
        plt.ylabel('Perimeter (in)')
        plt.title('Perimeter Change')
        plt.subplot(2,2,4)
        area_axes = plt.plot(times, areas)
        plt.xlabel('Time (s)')
        plt.ylabel('Area (in$^2$)')
        plt.title('Interior Area Change')
        plt.show()
        print('done plotting')
        return [times, areas, perimeters]

def hashfile(filepath):
    """Generate SHA2-256 hash of file at given path; return lowercase string."""
    return subprocess.check_output(['sha256sum', filepath]).decode('utf-8').split()[0]

if __name__ == '__main__':
    filename = sys.argv[1]
    filehash = ''
    while not filehash == hashfile(filename): # repeat if file has changed (e.g., in external editor)
        filehash = hashfile(filename)
        grain = FuelGrainState(filename=filename, rate=0.5, scale=0.01)
        calculate_regression(grain, timestep=0.2, num_steps=20)
