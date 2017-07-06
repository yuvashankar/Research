#include "wavelet.h"

// double CompleteRealMorlet (double x, double scale);
// double CompleteComplexMorlet(double x, double scale);

int CWT_Convolution(double *data, double * scales, int data_size, int num_of_scales, double d_t,
    double* result)
{
	double *conWindow, *complexWindow;
    double *realResult, *complexResult;

	conWindow     =  (double*) malloc(data_size * sizeof(double));
    complexWindow =  (double*) malloc(data_size * sizeof(double));
    realResult    =  (double*) malloc(data_size * sizeof(double));
    complexResult =  (double*) malloc(data_size * sizeof(double));

	int conSize = 0;
    for (int i = 0; i < num_of_scales; ++i)
    {
        memset(conWindow,     0.0, data_size * sizeof(double));
        memset(complexWindow, 0.0, data_size * sizeof(double));
        memset(realResult,    0.0, data_size * sizeof(double));
        memset(complexResult, 0.0, data_size * sizeof(double));
        
        conSize = (int) W_0/(2 * M_PI * scales[i]);
        // printf("Scale[%d] = %f, conSize = %d\n", i, scales[i], conSize);
        
        double temp = 0.0;
        //Populate Convolution windows.
        for (int j = 0; j < conSize; ++j)
        {
            conWindow[j] =       CompleteRealMorlet(temp, scales[i]);
            complexWindow[j] = - CompleteComplexMorlet(temp, scales[i]); //Complex Conjucate
            temp += d_t;
        }

        Convolute(data, conWindow, complexWindow, data_size, conSize,
            	  realResult, complexResult);

        for (int j = 0; j < data_size; ++j)
        {
            result[i * data_size + j] = MAGNITUDE(realResult[j], complexResult[j]);
        }
    }
    
    free(conWindow); free(complexWindow);
    free(realResult); free(complexResult);
    return 0;
}

void Convolute(double *data, double *conWindow, double * complexWindow, int data_size, int conSize,
	double* realResult, double* complexResult)
{
	for (int i = 0; i < data_size; ++i) //For every element in the data file
    {
        realResult[i]    = 0.0;
        complexResult[i] = 0.0;
        for (int j = -conSize + 1; j < conSize; ++j) 
        {
            if ( (i - j) >= 0 )
            {
                // printf("data[%d - %d] = %d, conWindow[%d] = %d\n", i,j,  data[i - j], j, conWindow[-j]);   
                if (j >= 0)
                {
                    realResult[i]    += data[i - j] * conWindow[j];
                    complexResult[i] += data[i - j] * complexWindow[j];
                }
                    
                if (j < 0)
                {
                    
                    realResult[i]     += data[i - j] * conWindow[-j];
                    complexResult[i] -= data[i - j] * complexWindow[-j];
                }
            }
            // else
            // {
            //     // int count = (i - j);
            //     int count = ( data_size + (i - j) )%data_size;
            //     if (j >= 0)
            //     {
            //         realResult[i]    += data[count] * conWindow[j];
            //         complexResult[i] += data[count] * complexWindow[j];
            //     }

            //     if (j < 0)
            //     {
            //         realResult[i]  += data[count] * conWindow[-j];
            //         complexResult[i] -= data[count] * conWindow[-j];
            //     }
            // }
        }
    }

}

double CompleteRealMorlet (double x, double scale)
{
	// double w_0 = 6.0;
	// double c_sigma = pow( (1.0 + exp(-w_0*w_0) - 2.0 * exp(-0.75 * w_0 * w_0)), -0.5 );
	// double k_sigma = exp( -0.5 * w_0 * w_0 );

	x = x / scale;

	// double normal = 1.0/sqrt(scale);

	double out = exp( -0.5 * x * x) * ( cos(W_0 * x) - K_SIGMA );
	out = C_SIGMA * QUAD_ROOT_PI * out;
	return(out);
}

double CompleteComplexMorlet(double x, double scale)
{
	// double w_0 = 6.0;
	// double c_sigma = pow( (1.0 + exp(-w_0*w_0) - 2.0 * exp(-0.75 * w_0 * w_0)), -0.5 );
	// double k_sigma = exp( -0.5 * w_0 * w_0 );

	x = x /scale; 
	// double normal = 1.0/sqrt(scale);
	double out = exp( -0.5 * x * x) * ( sin(W_0 * x) - K_SIGMA );

	out = C_SIGMA * QUAD_ROOT_PI * out;
	return(out);
}