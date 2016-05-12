#include "_generatefracincrements.h"

Increments::Increments(int N, int particles) {
    this->N = N;
    this->particles = particles;
    increments.resize(boost::extents[this->particles][3][this->N]);
}

void Increments::generateIncrements (double D ,double tau, double alpha)
{
        int dimensions= 3;
        int x = 1;
        int Nextended =N*x;

        fftw_plan planforward, planbackward;

        std::complex <double> *zw = new std::complex<double>[Nextended];
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
        double sigma = sqrt(tau);
        std::complex <double> pot_alpha1 (0,1);
        int i;
        int ipart;
        int idim;

       for (ipart=0; ipart < particles; ++ipart)
       {
            for (idim=0; idim < dimensions; ++idim)
                {
                //std::cout << idim;
                for(i = 0; i < Nextended; ++i)
                {
                    if( i<Nextended/2)
                    {

                        zw[i] = std::pow(((2*i*M_PI)/(Nextended*tau))*pot_alpha1, 1. - alpha)* D * tgamma(1.+ alpha);
                    }

                    if(i>=Nextended/2)
                    {

                        zw[i] = std::pow((2*((Nextended-i)*M_PI)/ (Nextended*tau))*pot_alpha1, 1. - alpha)* D * tgamma(1.+ alpha);
                    }

                    in[i] = gsl_ran_gaussian(r, sigma);



                }
                //std::cout << in[1] <<"," ;
                fftw_execute(planforward);

                innew[0]=gsl_ran_gaussian(r, sqrt(2.*D*pow(Nextended*tau,alpha)));

                for(i=1; i < Nextended; i++)
                {
                    innew[i]=out[i]*sqrt(2.0*zw[i].real());
                }
                fftw_execute(planbackward);

                for(i=0; i < N; i++)
                {
                    increments[ipart][idim][i]=outnew[i].real()/Nextended;
                    //increments[ipart][idim][i]=i+ipart+idim;
                    //std::cout << ipart <<"," <<  idim <<","<< i << "__";

                }

            }
       }
       fftw_destroy_plan(planforward);
       fftw_destroy_plan(planbackward);
       gsl_rng_free(r);

}
