#include<Python.h>
#include<stdio.h>
#include<vector>
#include<complex.h>
#include</home/mi/janekg89/Dokumente/Masterarbeit/test_cython/include/fftw3.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<stdlib.h>




std:: vector <double> _generateIncrements (int N, double* D ,double* tau, double* alpha,int particles)
{

        int dimensions= 3;
        double d3increments [particles][N][dimensions];
        //vector<vector<vector<double>>> d3increments;
        //d3increments.resize(particles);
        //for (int a=0;a<particles;a++)
        //{
	    //    d3increments[a].resize(dimensions);
	    //    for (int b=0;b<dimensions;b++)
		//            d3increments[a][b].resize(N);
        //}
        fftw_complex  *out, *zw, *innew, *outnew, *in;
        int x = 3;
        int Nextended =N*x;
        fftw_plan planforward, planbackward;
        zw = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nextended);
        innew = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nextended);
        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nextended);
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nextended);
        outnew = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nextended);
        planforward = fftw_plan_dft_1d(Nextended, in, out,FFTW_FORWARD, FFTW_ESTIMATE);
        planbackward = fftw_plan_dft_1d(Nextended, innew, outnew,FFTW_BACKWARD, FFTW_ESTIMATE );
        gsl_rng * r;
        const gsl_rng_type * T;
        gsl_rng_env_setup();
        gsl_rng_default_seed = time(NULL);
        T = gsl_rng_default;
        r = gsl_rng_alloc (T);
        double K_alpha = *D;
        double dt= *tau;
        double sigma = sqrt(*tau);
        std::vector <double> result(N);
        double complex pot_alpha;
        int i;
        int ipart;
        int idim;

       //for (ipart=0;ipart < particles;i++)
       //{

            //for (idim=0; idim < dimensions;i++)
                //{

                for(i = 0; i < Nextended; i++)
                {


                    if( i<Nextended/2)
                    {
                        pot_alpha = cpow(((2*i*M_PI)/(Nextended*dt))*I, 1. - *alpha);
                        zw[i] = pot_alpha*K_alpha*tgamma(1.+ *alpha);
                    }

                    if(i>=Nextended/2)
                    {
                        pot_alpha = cpow((2*((Nextended-i)*M_PI)/ (Nextended*dt))*I, 1. - *alpha);
                        zw[i] = pot_alpha*K_alpha*tgamma(1.+ *alpha);
                    }

                    in[i] = gsl_ran_gaussian(r, sigma)+0*I;
                }
                fftw_execute(planforward);
                innew[0]=gsl_ran_gaussian(r, sqrt(2.*K_alpha*pow(Nextended*dt,*alpha)))+0*I;
                //innew[0]=np.random.normal(0,np.sqrt(2.*self.K_alpha*(self.ki*self.dt)**self.alpha))
                for(i=1; i < Nextended; i++)
                {
                        innew[i]=out[i]*sqrt(2.0*creal(zw[i]));

                        //result.at(i)=creal(innew[i]);
                        //result.at(i)=creal(zw[i]);
                }
                fftw_execute(planbackward);

                for(i=0; i < N; i++)
                {
                    //result.at(i)=creal(zw[i]);
                    //result.at(i)=creal(cpow((M_PI*i)/(double) N*I,1.-*alpha))*K_alpha;
                    //result.at(i)=tgamma(*alpha+1);
                    //result.at(i)=((M_PI* i)/(double) N);
                    //result.at(i)=*tau;
                    //result.at(i)=*D;
                    //result.at(i)=*alpha;
                    result.at(i)=creal(outnew[i]/Nextended);
                    //d3increments[ipart][i][idim]=creal(outnew[i]/Nextended);

                }

            //}
       //}
       fftw_destroy_plan(planforward);
       fftw_destroy_plan(planbackward);
       fftw_free(in);
       fftw_free(out);
       fftw_free(innew);
       fftw_free(outnew);
       gsl_rng_free(r);
       return result;
}
    