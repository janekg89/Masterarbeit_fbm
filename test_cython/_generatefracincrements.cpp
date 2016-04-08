#include<Python.h>
#include<stdio.h>
#include<complex.h>
#include<vector>
#include</home/mi/janekg89/Dokumente/Masterarbeit/test_cython/include/fftw3.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>



 double _generateIncrements (int N, double* D,double* tau, double* alpha)
{

        int i;


        fftw_complex *in, *out;
        fftw_plan plan;

        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);

        plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);


        // = D * tau;
        const gsl_rng_type * T;
        gsl_rng * r;

        double sigma = sqrt(*tau);
        std::vector <double> result(N);
        //complex <double> *z;

        gsl_rng_env_setup();


        T = gsl_rng_default;
        r = gsl_rng_alloc (T);

        for(i = 0; i < N; i++)
        {

            in[i] = gsl_ran_gaussian(r, sigma)+I*0.;
            //in[i] = (i+1.)+(3.*i-1.)*I;
        }

        fftw_execute(plan);
        fftw_destroy_plan(plan);

        for(i=0; i < N; i++)
        {
            //in [i]=out[i];
            result.at(i)=creal(out[i]);

        }

        gsl_rng_free(r);
        fftw_free(in);
        fftw_free(out);
        //return result.at(1);
        return creal(out[1]);

        //printf(out[1][0].type);


}
    