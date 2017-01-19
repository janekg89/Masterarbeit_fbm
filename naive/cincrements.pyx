
import numpy as np
cimport numpy as np
np.import_array()

cdef extern from "_generatefracincrements.h":
    cdef cppclass Increments:
        Increments(int, int)
        void generateIncrements(double, double, double)
        double * returnIncrements()

cdef class pyIncrements:
    cdef Increments *thisptr
    cdef int shape0
    cdef int shape1
    def __cinit__(self, n, particles):
        self.thisptr = new Increments(n, particles)
        self.shape0 = particles
        self.shape1 = n
    def __dealloc__(self):
        del self.thisptr
    def generateIncrements(self, D, tau, alpha):
        self.thisptr.generateIncrements(D, tau, alpha)
    def returnIncrements(self):
        cdef np.npy_intp shape[3]
        shape[0] = <np.npy_intp> self.shape0
        shape[1] = <np.npy_intp> 3
        shape[2] = <np.npy_intp> self.shape1
        # Create a 1D array, of length 'size'
        data_ptr = <void*> self.thisptr.returnIncrements()
        ndarray = np.PyArray_SimpleNewFromData(3, shape, np.NPY_DOUBLE, data_ptr)
        return ndarray
