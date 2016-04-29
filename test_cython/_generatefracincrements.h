//double *_generateIncrements (int N, double* D,double* tau, double* alpha,int particles);
#ifndef INCREMENTS_INCLUDED
#define INCREMENTS_INCLUDED

#include <boost/multi_array.hpp>
#include <complex>
#include </home/mi/janekg89/Dokumente/Masterarbeit/test_cython/include/fftw3.h>
#include<stdio.h>
#include <iostream>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<stdlib.h>
#include <algorithm>


class Increments {
public:
    Increments(int N, int particles);
    int N;
    int particles;
    boost::multi_array<double,3> increments;
    void generateIncrements(double D,double tau, double alpha);
    double * returnIncrements() { return increments.data(); }
};

#endif // INCREMENTS_INCLUDED