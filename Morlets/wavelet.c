//Wavelet.c
#include "Morlet.h"
#include <omp.h>

int Wavelet(double* raw_data, double sampling_frequency, int n, double dj, double s0, int J, 
	double* result, double* frequency)
{
	
	//Variable Declarations
	double value;
	

	// double *filter; //Un-comment to look at each filter
	fftw_plan plan_forward, plan_backward;
	fftw_complex *data_in, *fft_data, *filter_convolution, *fftw_result;

	//An ouptut file for debugging. 
	// FILE *debug_file=fopen("debug.log","w");
	// assert(debug_file != NULL);
	//Memory Allocations
    // filter = malloc(PADDED_SIZE * J * sizeof(double));
    // assert(filter != NULL);

	//Calculate Padding Required
	const int pad = floor(log(n)/log(2.0) + 0.499);
    const double PADDED_SIZE = pow(2, pad + 1);

    //FFTW allocations.
    data_in = 			 fftw_alloc_complex(PADDED_SIZE); 
    fft_data = 			 fftw_alloc_complex(PADDED_SIZE);
    filter_convolution = fftw_alloc_complex(PADDED_SIZE);
    fftw_result = 		 fftw_alloc_complex(PADDED_SIZE);

    //populate the data vector. 
    for (int i = 0; i < n; ++i)
    {
    	data_in[i][0] = raw_data[i];
    	data_in[i][1] = 0.0;
    }

    const double df = sampling_frequency/PADDED_SIZE;
    const double k = exp(-0.5 * W_0_2);
    const double cSigma = pow(1.0 + exp(-W_0_2) - 2*exp(-0.75*W_0_2), -0.5);
    
    //Calculate the FFT for the Data
	plan_forward = fftw_plan_dft_1d(PADDED_SIZE, data_in, fft_data, 
		FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan_forward);

	//Copy the data into filter Convolution
	memcpy(filter_convolution, fft_data, (PADDED_SIZE * sizeof(fftw_complex)));

	//Preapre for the plan backwards
	plan_backward = fftw_plan_dft_1d(PADDED_SIZE, filter_convolution, fftw_result, 
		FFTW_BACKWARD, FFTW_ESTIMATE);

	double FOURIER_WAVELENGTH_FACTOR = (8 * M_PI)/(W_0);
	// FOURIER_WAVELENGTH_FACTOR *= (n*dt)/pow(2, (n/sampling_frequency)/2 - 1);
	printf("Fourier Wavelength Factor = %f\n", FOURIER_WAVELENGTH_FACTOR);
	for (int i = 0; i < J; ++i)
	{
		//Calculate the scale and frequency at the specific Scale
		double scale = s0 * pow(2, i * dj);

		// frequency[i] = scale * (DATA_SIZE * M_PI)/(sampling_frequency * W_0);
		frequency[i] = scale * FOURIER_WAVELENGTH_FACTOR;

		//Normalization Factor needes to be recomputed at every scale.
		// double normal = sqrt((2 * M_PI * scale)/(dt));
		double normal = sqrt(2 * M_PI * scale * sampling_frequency);

		//Caluclate the Fourier Morlet at the specific scale. 

		for (int j = 0; j < PADDED_SIZE; ++j)
		{
			value = FourierMorlet(j*df, scale, k, cSigma, normal);

			filter_convolution[j][0] *= value;
			filter_convolution[j][1] *= value;
		}

		//Take the inverse FFT. 
		fftw_execute(plan_backward);

		//Calculate the power and store it in result
		for (int j = 0; j < n; ++j)
		{
			result[i * n + j] = Magnitude(fftw_result[j][0], fftw_result[j][1]);
		}

		//Copy the fft_data into a seperate filter_convolution 
		memcpy(filter_convolution, fft_data, (PADDED_SIZE * sizeof(fftw_complex)));
	}

	//FFTW sanitation engineering. 
    fftw_destroy_plan(plan_forward); fftw_destroy_plan(plan_backward);
    fftw_free(data_in); fftw_free(fft_data); fftw_free(fftw_result);
    fftw_free(filter_convolution);

    // free(filter);
    // fclose(debug_file);

    return(0);
}

double FourierMorlet(double w, double scale, double k, double cSigma,
	double normal)
{
	w = w/scale;
	const double w2 = w * w;
	////These are all needed by Fourier Morlet i'm going to move 
	////them out to optimize the code
	// const double k = exp(-0.5 * W_0_2);
	// const double normal = sqrt((2 * M_PI * scale)/(1.0/sampling_frequency));
	// const double cSigma = pow(1.0 + exp(-W_0_2) - 2*exp(-0.75*W_0_2), -0.5);
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
