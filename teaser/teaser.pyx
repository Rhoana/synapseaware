import time


cimport cython
import ctypes
cimport numpy as np
import numpy as np


from synapseaware.utilities import dataIO



cdef extern from 'cpp-teaser.h':
    void CppTeaserSkeletonization(const char *prefix, long label, float resolution[3])



# use TEASER algorithm to generate skeletons
def TEASER(prefix, label):
    start_time = time.time()

    # call the teaser skeletonization algorithm
    cdef np.ndarray[float, ndim=1, mode='c'] cpp_resolution = np.ascontiguousarray(dataIO.Resolution(prefix)).astype(np.float32)
    CppTeaserSkeletonization(prefix.encode('utf-8'), label, &(cpp_resolution[0]))

    print ('TEASER skeletonization time in {:0.2f} seconds'.format(time.time() - start_time))
