#include "_generatefracincrements.h"

Increments1::Increments1(int N, int particles) {
    this->N = N;
    this->particles = particles;
    increments.resize(boost::extents[this->particles][3][this->N]);
}
/*
namespace janek {
    boost::multi_array<double,2> generateIncrements(args, Random * random) {

        // random->normal();
        // random->uniform();
    }
}
*/
double covariance(long i, double alpha) {
    if (i == 0) return 1;
    else return (pow(i-1,alpha)-2*pow(i,alpha)+pow(i+1,alpha))/2;
}
void Increments1::generateIncrements1 (double D ,double tau, double alpha)
{
        int dimensions= 3;
        int Nextended = N*2;

        fftw_plan planforward, planbackward;

        //std::complex <double> *r_x = new std::complex<double>[Nextended];
        std::complex <double> *innew = new std::complex<double>[Nextended];
        std::complex <double> *in = new std::complex<double>[Nextended];
        std::complex <double> *out = new std::complex<double>[Nextended];
        std::complex <double> *outnew = new std::complex<double>[Nextended];

        planforward = fftw_plan_dft_1d(Nextended, reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out),FFTW_FORWARD, FFTW_ESTIMATE);
        planbackward = fftw_plan_dft_1d(Nextended, reinterpret_cast<fftw_complex*>(innew), reinterpret_cast<fftw_complex*>(outnew),FFTW_BACKWARD, FFTW_ESTIMATE );

        gsl_rng * r;
        const gsl_rng_type * T;
        gsl_rng_env_setup();
        T = gsl_rng_ranlxs0;
        r = gsl_rng_alloc (T);
        gsl_rng_set(r, 0);
        std::complex <double> icomplex (0,1);
        double factor= sqrt(D* std::pow(tau ,alpha))/N;










        for(int ii = 0; ii < Nextended; ++ii)

        {
                    if(ii<=N)
                    {
                     in[ii] =  covariance(ii,alpha);
                    }

                    if(ii>N)
                    {
                     in[ii] = in[2*N-(ii)];
                    }
        }
         in[N]=0;
        /*
        for (ipart=0; ipart < particles; ++ipart)
       {
            for (idim=0; idim < dimensions; ++idim)
                {
                //std::cout << idim;
                for(i = 0; i < Nextended; ++i)

                    {
                    increments[ipart][idim][i]=in[i].real();
                    }
                }
       }
       */
       fftw_execute(planforward);



       for (int ipart=0; ipart < particles; ++ipart)
       {
            for (int idim=0; idim < dimensions; ++idim)
                {
                //std::cout << idim;
                for(int i = 0; i < Nextended; ++i)
                {
                    if( i<N)
                    {
                        innew[i] = sqrt(out[i].real()*(double)N)*gsl_ran_gaussian(r,1)*std::exp(2*M_PI*gsl_rng_uniform(r)*icomplex);
                    }


                    if(i>N)
                    {
                        innew[i] = std::conj(innew[2*N-i]);
                    }

                }
                    innew[0]=sqrt(out[0].real()*(double)N)*gsl_ran_gaussian(r,1);
                    innew[N]=sqrt(out[N].real()*(double)N)*gsl_ran_gaussian(r,1);

                    fftw_execute(planbackward);


                 for(int i=0; i < N; i++)
                 {
                    increments[ipart][idim][i]=(outnew[i].real())*factor;
                     //increments[ipart][idim][i]=in[i].real();
                 }

            }
       }

       fftw_destroy_plan(planforward);
       fftw_destroy_plan(planbackward);
       gsl_rng_free(r);

}
