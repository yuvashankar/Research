//Wavelet.c
//Starting over and trying to implement the Wavelet Function in this new file. 
#include "Morlet.h"

int Wavelet(double* raw_data, double dt, int n, double dj, double s0, int J, 
	double* result, double* frequency)
{
	
	//Variable Declarations
	double value;

	double *filter;
	fftw_plan plan_forward, plan_backward;
	fftw_complex *data_in, *fft_data, *filter_convolution, *fftw_result;

	//An ouptut file for debugging. 
	FILE *debug_file=fopen("debug.log","w");
    assert(debug_file != NULL);

	//Calculate Padding Required
	//Where did that 0.4999 come from I don't know
	int pad = floor(log(DATA_SIZE)/log(2.0) + 0.5);
    double PADDED_SIZE = pow(2, pad + 1);

    //Memory Allocations
    filter = malloc(PADDED_SIZE * J * sizeof(double));
    assert(filter != NULL);

    data_in = fftw_alloc_complex(PADDED_SIZE); 
    fft_data = fftw_alloc_complex(PADDED_SIZE);
    filter_convolution = fftw_alloc_complex(PADDED_SIZE);
    fftw_result = fftw_alloc_complex(PADDED_SIZE);

    //populate the data vector. 
    for (int i = 0; i < n; ++i)
    {
    	data_in[i][0] = raw_data[i];
    	data_in[i][1] = 0.0;
    }

	//Compute the fourier transform of the data signal
	plan_forward = fftw_plan_dft_1d(PADDED_SIZE, data_in, fft_data, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan_forward);

	plan_backward = fftw_plan_dft_1d(PADDED_SIZE, filter_convolution, fftw_result, FFTW_BACKWARD, FFTW_ESTIMATE);


	double df = 1.0/n/dt;
	double scale; 

	double fourier_wavelength_factor = (4 * M_PI)/(W_0 + sqrt(2 + W_0_2));
	// double fourier_wavelength_factor = (5.0/FS) * 2 * M_PI;
	for (int i = 0; i < J; ++i)
	{
		scale = s0 * pow(2, i*dj);
		frequency[i] = scale * fourier_wavelength_factor;
		printf("i is: %d, Scale is: %f, frequency is: %f\n", i, scale, frequency[i]);

		for (int j = 0; j < n; ++j)
		{
			filter[i*n + j] = NewFourierMorlet(j*df, W_0, scale, n);
			filter_convolution[j][0] = fft_data[j][0] * filter[i * n + j];
			filter_convolution[j][1] = 0.0;
		}

		//copy the rest of fft_data into filter_convolution
		for (int j = n/2; j < PADDED_SIZE; ++j)
		{
			filter_convolution[j][0] = fft_data[j][0];
			filter_convolution[j][1] = fft_data[j][1];
		}

		//Take the inverse FFT. 
		fftw_execute(plan_backward);

		//Copy to result array
		for (int j = 0; j < n; ++j)
		{
			value = Magnitude(fftw_result[j][0], fftw_result[j][1]);
			result[i * n + j] = value;
		}
	}

	//Clean things up... this may not be needed because this is a function
	//FFTW sanitation. 
    fftw_destroy_plan(plan_forward); fftw_destroy_plan(plan_backward);
    fftw_free(data_in); fftw_free(fft_data); fftw_free(fftw_result);
    fftw_free(filter_convolution);

    free(filter);
    fclose(debug_file);

    return(0);
}

double NewFourierMorlet(double w, double w0, double scale, int n)
{
	w = w/scale;

	//Wikipedia's Definition
	const double w02 = w0 * w0;
	const double k = exp(-0.5 * w02);
	double cSigma = pow(1. + exp(-w02) - 2*exp(-0.75*w02), -0.5);


	const double normal = sqrt((2 * M_PI * scale)/(1.0/FS));

	double out = exp( -0.5 * (w0 - w)*(w0 - w)) - k * exp(-0.5 * w*w);
	out = cSigma * quadRootPi * normal * out;
	return(out);
}

double Magnitude (double x, double y)
{
	double output = x * x + y * y;
	output = sqrt(output);
	return (output);
}