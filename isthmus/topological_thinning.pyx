import os
import time


cimport cython
cimport numpy as np
import ctypes
import numpy as np


from synapseaware.utilities import dataIO


cdef extern from 'cpp-topological-thinning.h':
    void CppTopologicalThinning(const char *prefix, const char *lookup_table_directory, long label, float resolution[3])


# generate skeletons for this volume
def TopologicalThinning(prefix, label):
    start_time = time.time()
    
    # convert the numpy arrays to c++
    lut_directory = os.path.dirname(__file__)

    # call the topological skeleton algorithm
    cdef np.ndarray[float, ndim=1, mode='c'] cpp_resolution = np.ascontiguousarray(dataIO.Resolution(prefix)).astype(np.float32)
    CppTopologicalThinning(prefix.encode('utf-8'), lut_directory.encode('utf-8'), label, &(cpp_resolution[0]))
    
    print ('Topological thinning time in {:0.2f} seconds'.format(time.time() - start_time))
