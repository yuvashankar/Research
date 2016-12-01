#include "wavelet.h"

double CompleteRealMorlet (double x, double scale)
{
	double w_0 = 6.0;
	double c_sigma = pow( (1.0 + exp(-w_0*w_0) - 2.0 * exp(-0.75 * w_0 * w_0)), -0.5 );
	double k_sigma = exp( -0.5 * w_0 * w_0 );

	x = x * scale; 

	double out = exp( -0.5 * x * x) * ( sin(w_0 * x) - k_sigma );
	out = c_sigma * QUAD_ROOT_PI * out;
	return(out);
}

double CompleteComplexMorlet(double x, double scale)
{
	double w_0 = 6.0;
	double c_sigma = pow( (1.0 + exp(-w_0*w_0) - 2.0 * exp(-0.75 * w_0 * w_0)), -0.5 );
	double k_sigma = exp( -0.5 * w_0 * w_0 );
	x = x * scale; 

	double out = exp( -0.5 * x * x) * ( cos(w_0 * x) - k_sigma );

	out = c_sigma * QUAD_ROOT_PI * out;
	return(out);
}

double CreateFilter(double * conWindow, double* complexWindow, double w0)
{
	double signalFrequency = w0/FS;
    double dw = 2 * M_PI * signalFrequency;

    int conSize = (int) 1./signalFrequency;
    conSize *=4;
    
    // printf("CompleteComplexMorlet = %f\n", CompleteComplexMorlet(0.0, 1.0));
    // FILE* debug_out = fopen("debug.log", "w");
    double mTime = 0.0;
    
    for (int i = 0; i < DATA_SIZE; ++i)
    {
        conWindow[i] = CompleteRealMorlet(mTime, 1.0);
        complexWindow[i] = CompleteComplexMorlet(mTime, 1.0);
        // fprintf(debug_out, "%f\t%f\t%f\t%f\n", mTime, val1, val2, data[i]);
        mTime += DT;
    }

    return(conSize);
}

void Convolute(double *data, double *conWindow, double * complexWindow, double conSize,
	double* result, double* complexResult)
{
	for (int i = 0; i < DATA_SIZE; ++i) //For every element in the data file
	{
		result[i] = 0.0;
		for (int j = -conSize; j < conSize; ++j) 
		{
			if ( (i - j) > 0)
			{
				if (j >= 0)
				{
					result[i] += data[i - j] * conWindow[j];
					complexResult[i] += data[i - j] * complexWindow[j];
				}
					

				if (j < 0)
				{
					result[i] += data[i - j] * conWindow[-j];
					complexResult[i] -= data[i - j] * complexWindow[-j];
				}
					
			}
		}
	}
}



