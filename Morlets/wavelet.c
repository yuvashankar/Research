//Wavelet.c
#include "Morlet.h"

int Wavelet(double* raw_data, double dt, int n, double dj, double s0, int J, 
	double* result, double* frequency)
{
	
	//Variable Declarations
	double value;

	// double *filter; //Un-comment to look at each filter
	fftw_plan plan_forward, plan_backward;
	fftw_complex *data_in, *fft_data, *filter_convolution, *fftw_result;

	//An ouptut file for debugging. 
	// FILE *debug_file=fopen("debug.log","w");
 //    assert(debug_file != NULL);

	//Calculate Padding Required
	//Where did that 0.4999 come from I don't know
	const int pad = floor(log(DATA_SIZE)/log(2.0) + 0.499);
    const double PADDED_SIZE = pow(2, pad + 1);

    //Memory Allocations
    // filter = malloc(PADDED_SIZE * J * sizeof(double));
    // assert(filter != NULL);

    //FFTW allocations.
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
    int plan_forward_flag = fftw_import_wisdom_from_filename(PLAN_FORWARD);

    //If there is no wisdom, then run the fft and save the wisedom. 
    if (plan_forward_flag == 0)
    {
    	printf("FFTW Plan forward does not exist, running tests and determining best FFT method...\n");
		plan_forward = fftw_plan_dft_1d(PADDED_SIZE, data_in, fft_data, FFTW_FORWARD, FFTW_EXHAUSTIVE);
		fftw_execute(plan_forward);
		plan_forward_flag = fftw_export_wisdom_to_filename(PLAN_FORWARD);
		assert(plan_forward_flag != 0);
    }

    int plan_backward_flag = fftw_import_wisdom_from_filename(PLAN_BACKWARD);
    if (plan_backward_flag == 0)
    {
    	printf("FFTW Plan Backward does not exist, running tests and determining best FFT method...\n");
    	plan_backward = fftw_plan_dft_1d(PADDED_SIZE, filter_convolution, fftw_result, FFTW_BACKWARD, FFTW_EXHAUSTIVE);
    	plan_backward_flag = fftw_export_wisdom_to_filename(PLAN_BACKWARD);
    	assert(plan_backward_flag != 0);
    }
	// plan_forward = fftw_plan_dft_1d(PADDED_SIZE, data_in, fft_data, FFTW_FORWARD, FFTW_ESTIMATE);
	// fftw_execute(plan_forward);

	


	const double df = FS/n;
	double scale;

	//This works but I don't know why... Where'd this number come from?
	const static double fourier_wavelength_factor = 1.8827925275534296252520792527491;

	// printf("Fourier Wavelength Factor = %f\n", fourier_wavelength_factor);

	for (int i = 0; i < J; ++i)
	{
		scale = s0 * pow(2, i * dj);
		frequency[i] = scale * fourier_wavelength_factor;

		//Copy the fft_data into a seperate filter_convolution 
		memcpy(filter_convolution, fft_data, (PADDED_SIZE * sizeof(fftw_complex)));

		//Caluclate the Fourier Morlet at the specific scale. 
		for (int j = 0; j < n/2; ++j)
		{
			filter_convolution[j][0] *= FourierMorlet(j*df, scale, n);
			filter_convolution[j][1] = 0.0;
		}

		//Take the inverse FFT. 
		fftw_execute(plan_backward);

		//Calculate the power and store it in result
		for (int j = 0; j < n; ++j)
		{
			result[i * n + j] = Magnitude(fftw_result[j][0], fftw_result[j][1]);
		}
	}

	//FFTW sanitation engineering. 
    fftw_destroy_plan(plan_forward); fftw_destroy_plan(plan_backward);
    fftw_free(data_in); fftw_free(fft_data); fftw_free(fftw_result);
    fftw_free(filter_convolution);

    // free(filter);
    // fclose(debug_file);

    return(0);
}

double FourierMorlet(double w, double scale, int n)
{
	//Wikipedia's Definition
	w = w/scale;
	const double w2 = w * w;
	const double k = exp(-0.5 * W_0_2);
	const double normal = sqrt((2 * M_PI * scale)/(1.0/FS));

	const double cSigma = pow(1.0 + exp(-W_0_2) - 2*exp(-0.75*W_0_2), -0.5);
	// double out = exp( -0.5 * (W_0 - w)*(W_0 - w)) - k * exp(-0.5 * w*w);
	double out = exp( -0.5 * (W_0_2 - 2*W_0*w + w2)) - k * exp(-0.5 * w2);

	out = cSigma * QUAD_ROOT_PI * normal * out;
	return(out);
}

double Magnitude (double x, double y)
{
	double output = (x * x) + (y * y);
	output = sqrt(output);
	return (output);
}