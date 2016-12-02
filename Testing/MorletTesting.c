#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

#define DATA_SIZE 8000
#define QUAD_ROOT_PI 0.7511255444649425 //Precalculated to machine precision
#define start -4.0
#define FS 1000

double CompleteFourierMorlet(double w, double scale)
{
	double w_0 = 6.0;
	double c_sigma = pow( (1.0 + exp(-w_0*w_0) - 2.0 * exp(-0.75 * w_0 * w_0)), -0.5 );
	double k_sigma = exp( -0.5 * w_0 * w_0 );

	double out = exp( -0.5 * ( w_0 - scale * w ) * (w_0 - scale * w) ) 
					- k_sigma * (exp ( -0.5 * scale * w * w));
	out = c_sigma * QUAD_ROOT_PI * out;
	return(out);
}

void CompleteMorlet (fftw_complex* x, double scale)
{
	double w_0 = 6.0;
	double c_sigma = pow( (1.0 + exp(-w_0*w_0) - 2.0 * exp(-0.75 * w_0 * w_0)), -0.5 );
	double k_sigma = exp( -0.5 * w_0 * w_0 );

	// double realMorlet = c_sigma * QUAD_ROOT_PI * exp(-0.5*t*t) * ( cos( w_0 * t ) - k_sigma);
	// double complexMorlet = c_sigma * QUAD_ROOT_PI * exp(-0.5*t*t) * ( sin( w_0 * t ) - k_sigma);
	double dt = (2.0 * abs(start))/DATA_SIZE;
	double t = start;
	for (int i = 0; i < DATA_SIZE; ++i)
	{
		
		x[i][0] = c_sigma * QUAD_ROOT_PI * exp(-0.5*t*t) * ( cos( w_0 * t ) - k_sigma);
		x[i][1] = c_sigma * QUAD_ROOT_PI * exp(-0.5*t*t) * ( sin( w_0 * t ) - k_sigma);
		t += dt;
	}

}


void main(void)
{
	//Initialize Variables
	fftw_complex *morlet_time, *morlet_frequency;
	fftw_plan plan_forward;
	FILE * output_file = fopen("data.txt", "w");

	//Allocate the memory
	morlet_time = (fftw_complex *) fftw_malloc (DATA_SIZE * sizeof(fftw_complex));
	morlet_frequency = (fftw_complex *) fftw_malloc (DATA_SIZE * sizeof(fftw_complex));

	//Make the FFTW plan
	plan_forward  = fftw_plan_dft_1d( DATA_SIZE, morlet_time, morlet_frequency, FFTW_FORWARD,
						FFTW_ESTIMATE|FFTW_PRESERVE_INPUT );

	CompleteMorlet(morlet_time, 1);

	fftw_execute(plan_forward);

	double dt = (2.0 * abs(start))/DATA_SIZE;
	printf("dt = %f\n", dt);
	double t = start;
	for (int i = 0; i < DATA_SIZE; ++i)
	{
		fprintf(output_file, "%d\t%f\t%f\n", i, morlet_time[i][0], morlet_time[i][1]);
		t += dt;
	}

	//Clear up memory
	fclose(output_file);
	fftw_free(morlet_time); fftw_free(morlet_frequency);
	fftw_destroy_plan(plan_forward);

}