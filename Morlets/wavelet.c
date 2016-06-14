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
	//The 5.0 should be changed to w0 ASAP.
	double fourier_wavelength_factor = (4 * M_PI)/(5.0 + sqrt(2 + 5.0 * 5.0));
	for (int i = 0; i < J; ++i)
	{
		scale = s0 * pow(2, i);
		// frequency[i] = scale * fourier_wavelength_factor;
		// printf("Scale is: %f, frequency is: %f\n", scale, frequency[i]);
		for (int j = 0; j < n; ++j)
		{
			filter[i*n + j] = NewFourierMorlet(j*df, 5.0, scale, n);
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

	// //Write to debug file.
	// for (int i = 0; i < n; ++i)
	// {
		// fprintf(debug_file, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", i,
		// 	filter[oldN*0 + i] + 0., filter[oldN*1 + i] + 5., filter[oldN*2 + i] + 10., 
		// 	filter[oldN*3 + i] + 15, filter[oldN*4 + i] + 20, filter[oldN*5 + i] + 25, 
		// 	filter[oldN*6 + i] + 30, filter[oldN*7 + i] + 35, filter[oldN*8 + i] + 40, 
		// 	filter[oldN*9 + i] + 45);

		// fprintf(debug_file, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", i,
		// 	result[oldN*0 + i] + 0., result[oldN*1 + i] + 5., result[oldN*2 + i] + 10., 
		// 	result[oldN*3 + i] + 15, result[oldN*4 + i] + 20, result[oldN*5 + i] + 25, 
		// 	result[oldN*6 + i] + 30, result[oldN*7 + i] + 35, result[oldN*8 + i] + 40, 
		// 	result[oldN*9 + i] + 45);
	// }

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