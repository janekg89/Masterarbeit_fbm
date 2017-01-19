#include<stdio.h>
#include<math.h>
#include<complex.h>									//This library is declared before fftw3.h
#include<fftw3.h>

int main(void)
{
	int i;
	int Npoints = 10;
	fftw_complex *in, *out;
	fftw_plan plan;

	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Npoints);			//allocating memory
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Npoints);		//allocating memory
	plan = fftw_plan_dft_1d(Npoints, in, out, FFTW_FORWARD, FFTW_ESTIMATE); 	//Here we set which kind of transformation we want to perform

	printf("\nCoefficcients of the expansion:\n\n");
	for(i = 0; i < Npoints; i++)
	{
		in[i] = (i+1.)+(3.*i-1.)*I;
		printf("%d %11.7f %11.7f\n", i, creal(in[i]), cimag(in[i]));		//creal and cimag are functions of complex.h 
	}
	printf("\n");

	fftw_execute(plan); 								//Execution of FFT

	printf("Output:\n\n");
	for(i = 0; i < Npoints; i++)
	{
		printf("%d %11.7f %11.7f\n", i, creal(out[i]), cimag(out[i]));
	}

	printf("\nCoefficcients of the expansion:\n\n");
	for(i = 0; i < Npoints; i++)
	{
		in[i] = exp(-i);
		printf("%d %11.7f %11.7f\n", i, creal(in[i]), cimag(in[i]));		//creal and cimag are functions of complex.h 
	}
	printf("\n");

	fftw_execute(plan); 								//Execution of FFT

	printf("Output:\n\n");
	for(i = 0; i < Npoints; i++)
	{
		printf("%d %11.7f %11.7f\n", i, creal(out[i]), cimag(out[i]));
	}

	
	fftw_destroy_plan(plan);							//Destroy plan
	fftw_free(in);			 						//Free memory
	fftw_free(out);			 						//Free memory
	return 0;
}

