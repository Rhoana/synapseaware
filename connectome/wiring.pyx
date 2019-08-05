import os
import time


cimport cython
cimport numpy as np
import ctypes
import numpy as np



from synapseaware.utilities import dataIO



cdef extern from 'cpp-wiring.h':
    void CppUpdateResolution(float resolution[3])
    void CppSkeletonGeneration(const char *prefix, long label, const char *lookup_table_directory)
    void CppSkeletonRefinement(const char *prefix, long label, double resolution[3])




# extract the wiring diagram for this prefix and label
def GenerateSkeleton(prefix, label):
    # start running time statistics
    start_time = time.time()

    # generate the widths 
    cdef np.ndarray[float, ndim=1, mode='c'] cpp_resolution = np.ascontiguousarray(dataIO.Resolution(prefix)).astype(np.float32)
    CppUpdateResolution(&(cpp_resolution[0]))

    # the look up table is in the synapseaware/connectome folder
    lut_directory = os.path.dirname(__file__)

    # call the topological skeleton algorithm
    CppSkeletonGeneration(prefix.encode('utf-8'), label, lut_directory.encode('utf-8'))

    # print out statistics for wiring extraction
    print ('Generated skeletons in {:0.2f} seconds'.format(time.time() - start_time))



# post process the volume to correct segment errors
def RefineSkeleton(prefix, label):
    # start running time statistics
    start_time = time.time()

    # get the resolution for this data
    cdef np.ndarray[double, ndim=1, mode='c'] cpp_resolution = np.ascontiguousarray(dataIO.Resolution(prefix))

    # call the post processing algorithm
    CppSkeletonRefinement(prefix.encode('utf-8'), label, &(cpp_resolution[0]))

    # print out statistics 
    print ('Refined skeletons in {:0.2f} seconds'.format(time.time() - start_time))
