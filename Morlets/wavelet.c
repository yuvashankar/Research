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

    int oldN = n;
    //Reallocate int to the padded size;
    n = (int) PADDED_SIZE;

    //Memory Allocations
    filter = malloc(n * J * sizeof(double));
    assert(filter != NULL);

    data_in = fftw_alloc_complex(PADDED_SIZE); 
    fft_data = fftw_alloc_complex(PADDED_SIZE);
    fftw_result = fftw_alloc_complex(PADDED_SIZE);

    //Prepare the plans. 
    plan_forward = fftw_plan_dft_1d(PADDED_SIZE, data_in, fft_data, FFTW_FORWARD, FFTW_ESTIMATE);
    

    //populate the data vector. 
    for (int i = 0; i < oldN; ++i)
    {
    	data_in[i][0] = raw_data[i];
    }

	//Making w step. 
	double k[oldN];
    for (int i = 0; i < oldN/2 + 1; ++i)
    {
        k[i] = (i * 2 * M_PI / (oldN * 1));
        // fprintf(debug_file, "%f\n", k[i]);
    }

    int counterVariable = oldN/2 - 1;
    for (int i = oldN/2 + 1; i < oldN; ++i)
    {
        k[i] = -k[counterVariable];
        counterVariable -- ;
        // fprintf(debug_file, "%f\n", k[i]);
    }

	//populate the scales array
	double scale[J];
	for (int i = 0; i < J; ++i)
	{
		// scale[i] = s0 * pow(2, i * dj);
		scale[i] = s0 * pow(2, i*dj);
		printf("scale = %f\n", scale[i]);
		// fprintf(debug_file, "%f\n", scale[i]);
	}

	//Compute the fourier transform of the data signal
	fftw_execute(plan_forward);

	for (int i = 0; i < oldN; ++i)
	{
		filter[i] = NewFourierMorlet(k[i], 5.0, 22.0, n);
		fft_data[i][0] = fft_data[i][0] * filter[i];
	}

	//Take theinverse FFT. 
	plan_backward = fftw_plan_dft_1d(PADDED_SIZE, fft_data, fftw_result, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(plan_backward);

	double value;
	for (int i = 0; i < oldN; ++i)
	{
		value = Magnitude(fftw_result[i][0], fftw_result[i][1]);
		fprintf(debug_file, "%d\t%f\t%f\n", i, raw_data[i], value);
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
	float kplus = 0.;
	if (w > 0) kplus = 1;
	double dw = (1 * 2 * M_PI / (n * 1));

	double exponent = -0.5 * kplus * (scale * w - w0) * (scale * w - w0);
	double norm = sqrt(scale * dw) * quadRootPi * sqrt(n);

	double out = norm * exp(exponent);
	out = out * kplus;
	return(out);
}