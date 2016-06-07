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
        T = gsl_rng_default;
        r = gsl_rng_alloc (T);
        std::complex <double> icomplex (0,1);
        double factor= sqrt(N*4.*D* std::pow((double) N ,alpha));
        int i;
        int ipart;
        int idim;
        for(i = 0; i <= Nextended; ++i)
        {
                    if(i<=N)
                    {
                     in[i] =  std::pow(tau,alpha)*( 1- std::pow((double) i/N,alpha))/2;
                    }

                    if(i>N)
                    {
                     in[i] = in[2*N-i];
                    }
                    //in[i] = gsl_ran_gaussian(r, sigma);
        }
       fftw_execute(planforward);

       in[0]= std::complex <double> (0,0);

       for (ipart=0; ipart < particles; ++ipart)
       {
            for (idim=0; idim < dimensions; ++idim)
                {
                //std::cout << idim;
                for(i = 0; i < Nextended; ++i)

                {
                    if( i>0  and  i<N)
                    {
                        innew[i] = sqrt(out[i])*gsl_ran_gaussian(r,1)*std::exp(2*M_PI*gsl_rng_uniform(r)*icomplex);
                    }

                    if(i>=N)
                    {
                        innew[i] = std::conj(in[2*N-i]);
                    }

                }
                fftw_execute(planbackward);

                for(i=0; i < N; i++)
                {
                    increments[ipart][idim][i]=(outnew[i].real()-outnew[i+1].real())*factor;
                    //increments[ipart][idim][i]=i+ipart+idim;
                    //std::cout << ipart <<"," <<  idim <<","<< i << "__";
                }

            }
       }

       fftw_destroy_plan(planforward);
       fftw_destroy_plan(planbackward);
       gsl_rng_free(r);

}
