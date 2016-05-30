//Wavelet.c
//Starting over and trying to implement the Wavelet Function in this new file. 
#include "Morlet.h"

int Wavelet(double* raw_data, double dt, int n, double dj, double s0, int J, double* result)
{
	
	//Variable Declarations
	double value;

	double* filter;
	fftw_plan plan_forward, plan_backward;
	fftw_complex *data_in, *fft_data, *fftw_result;

	//An ouptut file for debugging. 
	FILE* debug_file=fopen("debug.log","w");
    assert(debug_file != NULL);

	//Calculate Padding Required
	//Where did that 0.4999 come from I don't know
	int pad = floor(log(DATA_SIZE)/log(2.0) + 0.499);
    double PADDED_SIZE = pow(2, pad + 1);

    int oldN = n;
    //Reallocate int to the padded size;
    n = (int) PADDED_SIZE;

    //Memory Allocations
    filter = malloc(n * J * sizeof(double));
    assert(filter != NULL);

    data_in = fftw_alloc_complex(PADDED_SIZE); 
    fft_data = fftw_alloc_complex(PADDED_SIZE);
    fftw_result = fftw_alloc_complex(PADDED_SIZE);

    //populate the data vector. 
    for (int i = 0; i < oldN; ++i)
    {
    	data_in[i][0] = raw_data[i];
    	data_in[i][1] = 0.0;
    }

	//Making omega steps 
	double k[oldN];
	double dw =  (2 * M_PI) / (oldN * 1);
	// value = -(oldN/2) * 2 * M_PI/ (oldN * 1);
	
	value = M_PI/2;
    for (int i = 0; i < oldN; ++i)
    {
    	k[i] = value;
    	value += dw;
    }

	//populate the scales array ... not used right now, i'm just doing one scale
	double scale[J];
	for (int i = 0; i < J; ++i)
	{
		
		scale[i] = s0 * pow(2, i*dj);
		// printf("scale = %f\n", scale[i]);
		
	}

	//Zeroing out everything just because of that fluctuation... 
	for (int i = 0; i < PADDED_SIZE; ++i)
	{
		fft_data[i][0] = 0.0;
		fft_data[i][1] = 0.0;
	}

	//Compute the fourier transform of the data signal
	plan_forward = fftw_plan_dft_1d(PADDED_SIZE, data_in, fft_data, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan_forward);

	double df = 1.0/oldN/dt;
	double sign = 1.0;
	for (int i = 0; i < oldN/2; ++i)
	{
		filter[i] = sign * NewFourierMorlet(i*df, 5.0, 22.0, n);
		fft_data[i][0] = fft_data[i][0] * filter[i];
		fft_data[i][1] = 0.0;

		fft_data[oldN - i - 1][0] = 0.0;
		fft_data[oldN - i - 1][1] = 0.0;

		sign *= -1.0;
	}


	//Take the inverse FFT. 
	plan_backward = fftw_plan_dft_1d(PADDED_SIZE, fft_data, fftw_result, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(plan_backward);

	//Write to debug file.
	for (int i = 0; i < oldN; ++i)
	{
		value = Magnitude(fftw_result[i][0], fftw_result[i][1]);

		fprintf(debug_file, "%d\t%f\t%f\t%f\t%f\n", i, fftw_result[i][0], fftw_result[i][1], value, filter[i]);
	}

	//Clean things up... this may not be needed because this is a function
	//FFTW sanitation. 
    fftw_destroy_plan(plan_forward); fftw_destroy_plan(plan_backward);
    fftw_free(data_in); fftw_free(fft_data); fftw_free(fftw_result);

    free(filter);
    fclose(debug_file);

    return(0);
}

double NewFourierMorlet(double w, double w0, double scale, int n)
{
	//Wikipedia's Definition
	const double w02 = w0 * w0;
	const double k = exp(-0.5 * w02);
	const double cSigma = sqrt((1. + exp(-w02) - 2*exp(-0.75*w02)));

	double out = exp( -0.5 * (w0 - w)*(w0 - w)) - k * exp(-0.5 * w*w);
	out = cSigma * out;
	return(out);
}

