import os
import struct



from numba import jit
import numpy as np
import scipy.spatial, scipy.optimize


from synapseaware.utilities import dataIO
from synapseaware.utilities.constants import *



def ReadPredictions(prefix, method, label):
    filename = '{}/{}/{:06d}-endpoints.pts'.format(method, prefix, label)
    if not os.path.exists(filename):
        SavePredictions(prefix, method, label)

    with open(filename, 'rb') as fd:
        zres, yres, xres, npoints, = struct.unpack('qqqq', fd.read(32))
        endpoints = struct.unpack('%sq' % npoints, fd.read(8 * npoints))
    
    return endpoints



@jit(nopython=True)
def CreateEndpoints(points, zres, yres, xres):
    endpoints = []

    for point in points:
        iz = point // (yres * xres)
        iy = (point - iz * yres * xres) // xres
        ix = point % xres

        # check all 26 neighbors
        nneighbors = 0
        for iw in range(iz - 1, iz + 2):
            if iw < 0 or iw > zres - 1: continue
            for iv in range(iy - 1, iy + 2): 
                if iv < 0 or iv > yres - 1: continue
                for iu in range(ix - 1, ix + 2):
                    if iu < 0 or iu > xres - 1: continue
                
                    neighbor_index = iw * yres * xres + iv * xres + iu

                    if neighbor_index in points:
                        nneighbors += 1

        if nneighbors == 2:
            endpoints.append(point)

    return endpoints

    


def CorrectIsthmusEndpoints(prefix, method):
    # go through every label for this method and dataset
    directory = '{}/{}'.format(method, prefix)

    labels = []
    for filename in sorted(os.listdir(directory)):
        if 'endpoints' in filename: continue
        labels.append(int(filename[:-4]))
    
    for label in labels:
        filename = '{}/{}/{:06d}.pts'.format(method, prefix, label)
        print (filename)

        zres, yres, xres = dataIO.GridSize(prefix)

        with open(filename, 'rb') as fd:
            zres, yres, xres, npoints = struct.unpack('qqqq', fd.read(32))
            points = list(struct.unpack('%sq' % npoints, fd.read(8 * npoints)))
        
        for ip in range(npoints):
            if points[ip] < 0: 
                points[ip] = -1 * points[ip]

        endpoints = CreateEndpoints(set(points), zres, yres, xres)

        output_filename = '{}/{}/{:06d}-endpoints.pts'.format(method, prefix, label)

        with open(output_filename, 'wb') as fd:
            nendpoints = len(endpoints)
            fd.write(struct.pack('qqqq', zres, yres, xres, nendpoints))
            fd.write(struct.pack('%sq' % nendpoints, *endpoints))

        

def SavePredictions(prefix, method, label):
    filename = '{}/{}/{:06d}.pts'.format(method, prefix, label)

    with open(filename, 'rb') as fd:
        zres, yres, xres, npoints, = struct.unpack('qqqq', fd.read(32))
        points = struct.unpack('%sq' % npoints, fd.read(8 * npoints))

    endpoints = []

    for point in points:
        if point < 0: endpoints.append(-1 * point)

    output_filename = '{}/{}/{:06d}-endpoints.pts'.format(method, prefix, label)

    with open(output_filename, 'wb') as fd:
        nendpoints = len(endpoints)
        fd.write(struct.pack('qqqq', zres, yres, xres, nendpoints))
        fd.write(struct.pack('%sq' % nendpoints, *endpoints))



def SynapseEvaluate(prefix, method, label):
    # go through every label for this method and dataset
    synapse_filename = 'synapses/{}/{:06d}.pts'.format(prefix, label)
    if not os.path.exists(synapse_filename): return
    endpoint_filename = '{}/{}/{:06d}.pts'.format(method, prefix, label)
    if not os.path.exists(endpoint_filename): return
    nri_filename = 'nris/{}/{}-{:06d}.txt'.format(prefix, method.replace('/', '-'), label)
    if os.path.exists(nri_filename): return

    # get the grid size
    zres, yres, xres = dataIO.GridSize(prefix)
    resolution = dataIO.Resolution(prefix)

    max_distance = 800

    # read the true synapse locations
    synapses = dataIO.ReadPoints(prefix, label, 'synapses')

    # read the predicted locations
    predictions = ReadPredictions(prefix, method, label)

    ngt_pts = len(synapses)
    npr_pts = len(predictions)
    

    gt_pts = np.zeros((ngt_pts, 3), dtype=np.int64)
    pr_pts = np.zeros((npr_pts, 3), dtype=np.int64)
    npoints = ngt_pts * npr_pts

    for pt in range(ngt_pts):
        # get x, y, z locations
        index = synapses[pt]

        iz = index // (yres * xres)
        iy = (index - iz * yres * xres) // xres
        ix = index % xres 

        # coordinates are (x, y, z)
        gt_pts[pt,0] = resolution[OR_X] * ix
        gt_pts[pt,1] = resolution[OR_Y] * iy
        gt_pts[pt,2] = resolution[OR_Z] * iz
    
    for pt in range(npr_pts):
        # get x, y, z locations
        index = predictions[pt]

        iz = index // (yres * xres)
        iy = (index - iz * yres * xres) // xres
        ix = index % xres

        # coordinates are (x, y, z)
        pr_pts[pt,0] = resolution[OR_X] * ix
        pr_pts[pt,1] = resolution[OR_Y] * iy
        pr_pts[pt,2] = resolution[OR_Z] * iz

    cost_matrix = scipy.spatial.distance.cdist(gt_pts, pr_pts)
    matching = scipy.optimize.linear_sum_assignment(cost_matrix)

    valid_matches = set()
    for match in zip(matching[0], matching[1]):
        # valid pairs must be within max_distance in nanometers
        if cost_matrix[match[0], match[1]] > max_distance: continue

        valid_matches.add(match)

    ncorrect_synapses = len(valid_matches)
    nadded_synapses = npr_pts - len(valid_matches)
    nmissed_synapses = ngt_pts = len(valid_matches)

    # the number of true positives is the number of paths between the valid locations
    true_positives = ncorrect_synapses * (ncorrect_synapses - 1) // 2
    # the number of false positives is every pair of paths between true and added synapses
    false_positives = ncorrect_synapses * nadded_synapses
    # the number of false negatives is every synapse pair that is divided
    false_negatives = ncorrect_synapses * nmissed_synapses

    if true_positives == 0:
        nri = 0
    else:
        precision = true_positives / float(true_positives + false_positives)
        recall = true_positives / float(true_positives + false_negatives)

        nri = 2 * (precision * recall) / (precision + recall)
    
    with open(nri_filename, 'w') as fd:
        fd.write('{} {} {}\n'.format(true_positives, false_positives, false_negatives))
        fd.write('{}\n'.format(nri))

    return nri
