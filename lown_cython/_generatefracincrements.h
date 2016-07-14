//double *_generateIncrements (int N, double* D,double* tau, double* alpha,int particles);
#ifndef INCREMENTS_INCLUDED
#define INCREMENTS_INCLUDED

#include <boost/multi_array.hpp>
#include <complex>
#include <fftw3.h>
#include <stdio.h>
#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdlib.h>
#include <algorithm>
//#include "Random.h"


class Increments1 {
public:
    Increments1(int N, int particles);

    int N;
    int particles;
    boost::multi_array<double, 3> increments;

    void generateIncrements1(double D, double tau, double alpha);

    double *returnIncrements() { return increments.data(); }
};

#endif // INCREMENTS_INCLUDED
