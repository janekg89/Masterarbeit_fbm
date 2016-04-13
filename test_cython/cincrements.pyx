cimport numpy
from libcpp.vector cimport vector
cdef extern from "_generatefracincrements.h":
     vector[double]_generateIncrements(int N, double* D,double* tau, double* alpha,int particles)

def generateIncrements(N, D, tau, alpha,particles):
    """
    Calculates the squared distance. x and p have to be numpy arrays of double with the same lenght that is given by dim for this to work
    :param dim:
        the dimension of x AND p
    :param x:
       first point
    :param p:
        second point
    :return: double distance
        square of the distance of both points
    """
    _D = <double*> numpy.PyArray_DATA(D)
    _tau = <double*> numpy.PyArray_DATA(tau)
    _alpha=<double*> numpy.PyArray_DATA(alpha)
    return _generateIncrements(N, _D, _tau, _alpha,particles)

