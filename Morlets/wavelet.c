//Wavelet.c
#include "wavelet.h"
#include <omp.h>

int Wavelet(double* raw_data,  double* frequency, 
	double sampling_frequency, int n, double dj, double s0, int J, double minimum_frequency,
	double* result)
{
	
	//Variable Declarations
	int i;
	int start = (int)floor( log2( (W_0 * minimum_frequency)
									 /(8 * M_PI * s0) )
									 /dj);
	// printf("outside J = %d\n", J);
	#pragma omp parallel num_threads(2) private(i) shared (result, frequency, raw_data, dj, s0, sampling_frequency, minimum_frequency, J, n, start) default(none)
	{
		int j;
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

	    //Things needed for FourierMorlet Calculated only once.
	    const double df = sampling_frequency/PADDED_SIZE;
	    const double k = exp(-0.5 * W_0_2);
	    const double cSigma = pow(1.0 + exp(-W_0_2) - 2*exp(-0.75*W_0_2), -0.5);
	    const double FOURIER_WAVELENGTH_FACTOR = (8 * M_PI)/(W_0);


	    //FFTW allocations.
	    data_in = 			 fftw_alloc_complex(PADDED_SIZE); 
	    fft_data = 			 fftw_alloc_complex(PADDED_SIZE);
	    filter_convolution = fftw_alloc_complex(PADDED_SIZE);
	    fftw_result = 		 fftw_alloc_complex(PADDED_SIZE);


		
		//populate the data vector. 
	    #pragma omp for
	    for (i = 0; i < n; ++i)
	    {
	    	data_in[i][0] = raw_data[i];
	    	data_in[i][1] = 0.0;
	    }

	    //FFTW allocations only one thread can do this. 
	    #pragma omp critical (make_plan)
		{
		//Calculate the FFT for the Data
		plan_forward = fftw_plan_dft_1d(PADDED_SIZE, data_in, fft_data, 
			FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(plan_forward);

		//Copy the data into filter Convolution
		memcpy(filter_convolution, fft_data, (PADDED_SIZE * sizeof(fftw_complex)));

		//Preapre for the plan backwards
		plan_backward = fftw_plan_dft_1d(PADDED_SIZE, filter_convolution, fftw_result, 
			FFTW_BACKWARD, FFTW_ESTIMATE);
	
		}
	    // printf("J = %d\n", J);
		#pragma omp for
		for (i = 0; i < J; ++i)
		{
			//Calculate the scale and frequency at the specific Scale
			double scale = s0 * pow(2, i * dj);
			frequency[i] = scale * FOURIER_WAVELENGTH_FACTOR;

			//Normalization Factor needes to be recomputed at every scale.
			double normal = sqrt(2 * M_PI * scale * sampling_frequency);

			//Caluclate the Fourier Morlet at the specific scale. 
			for (j = 0; j < PADDED_SIZE; ++j)
			{
				value = FourierMorlet(j*df, scale, k, cSigma, normal);

				filter_convolution[j][0] *= value;
				filter_convolution[j][1] *= value;
			}

			//Take the inverse FFT. 
			fftw_execute(plan_backward);

			//Calculate the power and store it in result
			for (j = 0; j < n; ++j)
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
	}
	

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

void FillData(double * data)
{
	// Fit a FREQ signal at two points
	// double dt = 1./FS;
	double fsig = FREQ/FS;
	double dw = 2*M_PI*fsig;
	double w0 =  0.01; // A SMALL PHASE SHIFT SO ITS NOT ALL INTERGER ALIGNED
	int one_peri = (int)1./fsig;
	printf("FS  %.2f   Pitch %.f   Discrete Period = %d \n",FS,FREQ,one_peri);

	for (int i = 0; i < DATA_SIZE; ++i)
	{
		data[i] = 0.0;
	}

	// //Impulse Sample
	// data[2000] = 1.0;

	///Sine Wave Sample
	int i;
	// double t=0;
	for(i=0;i<DATA_SIZE;i++){
		// data[i] = sin(i*dw) + sin(i*dw*4);
		data[i]=0.;
		// if((i>200)&(i<400))data[i]=sin( (i-200)*dw+w0);
		if((i>0.25*DATA_SIZE)&(i<0.25*DATA_SIZE+one_peri)) data[i]=sin( (i-200)*dw+w0);
		if((i>0.5*DATA_SIZE)&(i<0.5*DATA_SIZE+2*one_peri))data[i]=sin( (i-1000)*dw+w0);
		if((i>0.75*DATA_SIZE)&(i<0.75*DATA_SIZE+3*one_peri))data[i]=sin( (i-2000)*dw+w0);
	}
}

void TestCases(double *data, int flag)
{

	// Fit a FREQ signal at two points
	// double dt = 1./FS;
	double fsig = FREQ/FS;
	double dw = 2*M_PI*fsig;
	double w0 =  0.01; // A SMALL PHASE SHIFT SO ITS NOT ALL INTERGER ALIGNED
	int one_peri = (int)1./fsig;
	printf("FS  %.2f   Pitch %.f   Discrete Priode = %d \n",FS,FREQ,one_peri);
	
	switch(flag)
	{
		//Impulse
		case 1:
			for (int i = 0; i < DATA_SIZE; ++i)
			{
				data[i] = 0.0;
			}

			data[1500] = 1.0;
			break;
		
		//Multiple Sines at t = 1500
		case 2:
			for (int i = 1500; i < 1500 + 2*one_peri; ++i)
			{
				data[i] = sin((i - 1500)* dw + w0) + sin((i - 1500)* 2* dw + w0);
			}
			break;
		//Multiple Sines at all times
		case 3:
			for (int i = 0; i < DATA_SIZE; ++i)
			{
				data[i] = sin(i*dw + w0) + sin(i*2*dw + w0);
			}
	}
}

int ReadFile(double data[], char filename[])
{
	FILE* signalFile = fopen(filename, "r");
	assert(signalFile != NULL);
	// obtain file size:
	fseek (signalFile , 0 , SEEK_END);
	long lSize = ftell (signalFile);
	rewind (signalFile);

	char * buffer = (char*) malloc(sizeof(char)*lSize);
	assert(buffer != NULL);

	size_t result = fread (buffer, 1, lSize, signalFile);
	assert(result == lSize);
	// puts(buffer);


	char * token = strtok(buffer, "\n");
	
    //Get input from text.
	int counterVariable = 0;
	while (token !=NULL)
    {
    	data[counterVariable] = atof(token);
    	counterVariable++;
        token = strtok (NULL, "\n");

    }
    fclose(signalFile);

    return (counterVariable);
}

int WriteFile(double *data, double *frequency, int x, int y, char filename[])
{

    FILE* out_file=fopen(filename,"w");
    if (out_file == NULL) return -1;

    //Xticks
    fprintf(out_file, "%d\t", x);
    for (int i = 0; i < y; ++i)
    {
    	fprintf(out_file, "%f\t", i/FS);
    }
    fprintf(out_file, "\n");


	for (int i = 0; i < x; ++i)
    {
    	//Feed Frequency
    	fprintf(out_file, "%f\t", frequency[i]);

    	//Feed Data
        for (int j = 0; j < y; ++j)
        {
            // value = Magnitude(result[i*n + j], result[i*n + j]);
            fprintf(out_file, "%f\t", data[i*y + j]);
        }
        //Ready for the next line.
        fprintf(out_file, "\n");

    }

    fclose(out_file);
    return(0);
}

int WriteTestCases(double *data, int length, char filename[])
{
	FILE* out_file=fopen(filename,"w");
    if (out_file == NULL) return -1;

	for (int i = 0; i < length; ++i)
    {
    	fprintf(out_file, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", i,
    		data[length*0 + i] + 0., data[length*1 + i] + 5., data[length*2 + i] + 10, 
    		data[length*3 + i] + 15, data[length*4 + i] + 20, data[length*5 + i] + 25, 
    		data[length*6 + i] + 30, data[length*7 + i] + 35, data[length*8 + i] + 40, 
    		data[length*9 + i] + 45, data[length*10+ i] + 50);
    }
    
    fclose(out_file);
    return 0;
}