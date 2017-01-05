#include "wavelet.h"

// double CompleteRealMorlet (double x, double scale);
// double CompleteComplexMorlet(double x, double scale);

double CompleteRealMorlet (double x, double scale)
{
	// double w_0 = 6.0;
	// double c_sigma = pow( (1.0 + exp(-w_0*w_0) - 2.0 * exp(-0.75 * w_0 * w_0)), -0.5 );
	// double k_sigma = exp( -0.5 * w_0 * w_0 );

	x = x * scale;

	double out = exp( -0.5 * x * x) * ( sin(W_0 * x) - K_SIGMA );
	out = C_SIGMA * QUAD_ROOT_PI * out;
	return(out);
}

double CompleteComplexMorlet(double x, double scale)
{
	// double w_0 = 6.0;
	// double c_sigma = pow( (1.0 + exp(-w_0*w_0) - 2.0 * exp(-0.75 * w_0 * w_0)), -0.5 );
	// double k_sigma = exp( -0.5 * w_0 * w_0 );

	x = x * scale; 
	double out = exp( -0.5 * x * x) * ( cos(W_0 * x) - K_SIGMA );

	out = C_SIGMA * QUAD_ROOT_PI * out;
	return(out);
}


void Convolute(double *data, double *conWindow, double * complexWindow, int data_size, int conSize,
	double* realResult, double* complexResult)
{
	for (int i = 0; i < data_size; ++i) //For every element in the data file
	{
		realResult[i] = 0.0;
		for (int j = -conSize; j < conSize; ++j) 
		{
			if ( (i - j) > 0)
			{
				if (j >= 0)
				{
					realResult[i] += data[i - j] * conWindow[j];
					complexResult[i] += data[i - j] * complexWindow[j];
				}
					

				if (j < 0)
				{
					realResult[i] += data[i - j] * conWindow[-j];
					complexResult[i] -= data[i - j] * complexWindow[-j];
				}
					
			}
		}
	}
}