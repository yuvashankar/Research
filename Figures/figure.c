#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define QUAD_ROOT_PI 0.75112554446494248
#define C_SIGMA 1.0000000000018794
#define K_SIGMA 1.522997974471262843e-08

#define W_0 6.0

#define MAGNITUDE(x,y) fabs( sqrt( (x * x) + (y * y) ) )

double CompleteFourierMorlet(double w, const double scale);
double Morlet(double t, double scale);

void main(void)
{
	int fs = 2048;
	double dt = 1.0/fs;
	int n = 6 * fs;
	
	double scale = 0.125;

	double* morlet   = (double*) malloc(n * sizeof(double));
	double* f_morlet = (double*) malloc(n * sizeof(double));

	double* morlet1   = (double*) malloc(n * sizeof(double));
	double* f_morlet1 = (double*) malloc(n * sizeof(double));

	double* morlet2   = (double*) malloc(n * sizeof(double));
	double* f_morlet2 = (double*) malloc(n * sizeof(double));

	double* morlet3   = (double*) malloc(n * sizeof(double));
	double* f_morlet3 = (double*) malloc(n * sizeof(double));

	const double dw = (2 * M_PI * fs)/(n);

	FILE* out = fopen("output.log", "w");
	for (int i = 0; i < n; ++i)
	{
		morlet[i]   = Morlet(i * dt - 3.0, scale);
		f_morlet[i] = CompleteFourierMorlet(i * dw, scale);

		morlet1[i]   = Morlet(i * dt - 3.0, scale * 2.0) + 3.0;
		f_morlet1[i] = CompleteFourierMorlet(i * dw, scale/2.0);
		
		morlet2[i]   = Morlet(i * dt - 3.0, scale * 4.0) + 2 * 3.0;
		f_morlet2[i] = CompleteFourierMorlet(i * dw, scale/4.0);

		morlet3[i]   = Morlet(i * dt - 3.0, scale * 8.0) + 3 * 3.0;
		f_morlet3[i] = CompleteFourierMorlet(i * dw, scale/8.0);

		fprintf(out, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", (double) i * dt, (double) (i * dw)/ (2.0 * M_PI), morlet[i], morlet1[i], morlet2[i], morlet3[i], f_morlet[i], f_morlet1[i], f_morlet2[i], f_morlet3[i]);
		// fprintf(out, "%f\t%f\t%f\t%f\n", (double) i * dt, morlet[i], (double) i * dw, f_morlet[i]);
	}

	free(morlet);
	free(f_morlet);

	free(morlet1);
	free(f_morlet1);

	free(morlet2);
	free(f_morlet2);

	free(morlet3);
	free(f_morlet3);
}


double Morlet(double t, double scale)
{
	double norm = 1.0/sqrt(scale);
	t = t/scale;
	double real_wave = exp(-(t * t)/2)*( cos(W_0 * t)- K_SIGMA);
	// double complex_wave = exp(-(t * t)/2)*( sin(W_0 * t)- K_SIGMA);

	// double wave = MAGNITUDE(real_wave, complex_wave);

	double out = norm * C_SIGMA * QUAD_ROOT_PI * real_wave;
	return(out);
}

double CompleteFourierMorlet(double w, const double scale)
{
	// double norm = 1.0/sqrt(scale);
	double norm = sqrt(scale);
	w = w * scale; 
	double out = exp( -0.5 * ( W_0 - w ) * (W_0 - w) )
                                       - K_SIGMA * (exp ( -0.5 * w * w));

	out = norm * C_SIGMA * QUAD_ROOT_PI * out;
	return(out);
}