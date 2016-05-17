//Wavelet.c
//Starting over and trying to implement the Wavelet Function in this new file. 
#include "Morlet.h"

int Wavelet(double* raw_data, double dt, int n, double dj, double s0, int J, double* result)
{
	
	//Variable Declarations
	double* filter;
	fftw_plan plan_forward, plan_backward;
	fftw_complex *data_in, *fft_data, *fftw_result;

	//An ouptut file for debugging. 
	FILE* debug_file=fopen("debug.log","w");
    assert(debug_file != NULL);

	//Calculate Padding Required
	int pad = floor(log(DATA_SIZE)/log(2.0) + 0.4999);
    double PADDED_SIZE = pow(2, pad + 1);

    //Memory Allocations
    filter = malloc(n * J * sizeof(double));
    assert(filter != NULL);

    data_in = fftw_alloc_complex(PADDED_SIZE); 
    fft_data = fftw_alloc_complex(PADDED_SIZE);
    fftw_result = fftw_alloc_complex(PADDED_SIZE);

    //Prepare the plans. 
    plan_forward = fftw_plan_dft_1d(PADDED_SIZE, data_in, fft_data, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_backward = fftw_plan_dft_1d(PADDED_SIZE, fft_data, fftw_result, FFTW_BACKWARD, FFTW_ESTIMATE);

	//Making w step. 
	double k[n];
    for (int i = 0; i < n/2 + 1; ++i)
    {
        k[i] = (i * 2 * M_PI / (n * 1));
        // fprintf(debug_file, "%f\n", k[i]);
    }

    int counterVariable = n/2 - 1;
    for (int i = n/2 + 1; i < n; ++i)
    {
        k[i] = -k[counterVariable];
        counterVariable -- ;
        // fprintf(debug_file, "%f\n", k[i]);
    }

	//populate the scales array
	double scale[J];
	for (int i = 0; i < J; ++i)
	{
		scale[J] = s0 * pow(2, i * dj);
	}

	//Compute the fourier transform of the data signal
	fftw_execute(plan_forward);

	//Populate the filter array
	for (int i = 0; i < J; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			filter[i * n + j] = NewFourierMorlet(k[j], 5.0, scale[i], n);
			//Only do the real part, not the complex part. 
			fft_data[j][0] *= filter[i * n + j];
		}
		//Compute the Fourier transform of xhat and sigmaHat this should put results in 
		fftw_execute(plan_backward);

		for (int j = 0; j < n; ++j)
		{
			//Take the magnitude of the real and the complex part, this is a guess. 
			result[i * n + j] = Magnitude(fftw_result[j][0], fftw_result[j][1]);
		}
	}


	//Clean things up... this may not be needed because this is a function
	//FFTW sanitation. 
    fftw_destroy_plan(plan_forward); fftw_destroy_plan(plan_backward);
    fftw_free(data_in); fftw_free(fft_data); fftw_free(fftw_result);

    fclose(debug_file);

    return(0);
}

double NewFourierMorlet(double w, double w0, double scale, int n)
{
	float kplus = 0.;
	if (w > 0) kplus = 1;
	double dw = (1 * 2 * M_PI / (n * 1));

	double exponent = -0.5 * kplus * (scale * w - w0) * (scale * w - w0);
	double norm = sqrt(scale * dw) * quadRootPi * sqrt(n);

	double out = norm * exp(exponent);
	out = out * kplus;
	return(out);
}