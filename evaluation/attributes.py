import math
import random
import time
import scipy.spatial



import numpy as np



from synapseaware.utilities import dataIO
from synapseaware.utilities.constants import *



def EvaluateWidths(prefix, label):
    start_time = time.time()
    
    # get the resolution, surface voxels, and radii for this prefix label pair
    resolution = dataIO.Resolution(prefix)
    zres, yres, xres = dataIO.GridSize(prefix)
    
    # surface information
    surface_point_cloud = np.array(dataIO.ReadPoints(prefix, label, 'surfaces'), dtype=np.int64)    
    npoints = len(surface_point_cloud)

    np_point_cloud = np.zeros((npoints, 3), dtype=np.int32)
    for index, iv in enumerate(surface_point_cloud):
        iz = iv // (yres * xres)
        iy = (iv - iz * yres * xres) // xres
        ix = iv % xres

        np_point_cloud[index,:] = (resolution[OR_X] * ix, resolution[OR_Y] * iy, resolution[OR_Z] * iz)
        index += 1


    widths = dataIO.ReadWidths(prefix, label)
    skeletons = dataIO.ReadPoints(prefix, label, 'connectomes')
    
    # keep track of the error over time
    mean_absolute_error = 0.0

    count = 0
    for iv, index in enumerate(skeletons):
        if random.random() < 0.80: continue

        # some of the soma locations will not be in the widths 
        if not index in widths: continue

        iz = index // (yres * xres)
        iy = (index - iz * yres * xres) // xres
        ix = index % xres

        # create a 2D vector for this point
        vec = np.zeros((1, 3), dtype=np.int32)
        vec[0,:] = (resolution[OR_X] * ix, resolution[OR_Y] * iy, resolution[OR_Z] * iz)

        # get the radius at this index 
        radius = widths[index]

        minimum_distance = scipy.spatial.distance.cdist(np_point_cloud, vec).min()
        
        error = abs(radius - minimum_distance)

        mean_absolute_error += error
        count += 1

    print ('Mean Absolute Error: {:0.2f} nanometers'.format(mean_absolute_error / count))

    print (time.time() - start_time)

    output_filename = 'width-errors/{}-{:06d}.txt'.format(prefix, label)
    with open(output_filename, 'w') as fd:
        fd.write('{} {}\n'.format(mean_absolute_error, count))

