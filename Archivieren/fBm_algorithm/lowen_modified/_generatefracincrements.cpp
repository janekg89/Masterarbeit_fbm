#include "_generatefracincrements.h"

Increments1::Increments1(int N, int particles) {
    this->N = N;
    this->particles = particles;
    increments.resize(boost::extents[this->particles][3][this->N]);
}


void Increments1::generateIncrements1 (double D ,double tau, double alpha)
{
        int dimensions= 3;
        int Nextended = N*2;

        fftw_plan planforward, planbackward;


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
        double factor= sqrt(D* std::pow((double) N ,alpha)/(N));
        for(int ii = 0; ii < Nextended; ++ii)
        {
                    if(ii<=N)
                    {
                     in[ii] =  std::pow(tau,alpha)*( 1- std::pow((double) ii/N,alpha))/2;
                    }

                    if(ii>N)
                    {
                     in[ii] = in[2*N-(ii)];
                    }
        }

       fftw_execute(planforward);

       innew[0]= std::complex <double> (0,0);
       innew[N]=sqrt(2*D* std::pow(tau ,alpha)*out[N])*gsl_ran_gaussian(r,1);

       for (int ipart=0; ipart < particles; ++ipart)
       {
            for (int idim=0; idim < dimensions; ++idim)
                {
                //std::cout << idim;
                for(int i = 0; i < Nextended; ++i)
                {
                    if( i>0  and  i<N)
                    {
                        innew[i] = sqrt(out[i]/2.0)*(gsl_ran_gaussian(r,1)+gsl_ran_gaussian(r,1)*icomplex);
                    }


                    if(i>N)
                    {
                        innew[i] = std::conj(innew[2*N-i]);
                    }

                }
                fftw_execute(planbackward);

                 for(int i=0; i < N; i++)
                 {
                    increments[ipart][idim][i]=(-outnew[i].real()+outnew[i+1].real())*factor;
                 }

            }
       }

       fftw_destroy_plan(planforward);
       fftw_destroy_plan(planbackward);
       gsl_rng_free(r);

}
