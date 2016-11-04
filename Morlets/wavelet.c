//Wavelet.c
#include "wavelet.h"
#include <omp.h>

int Wavelet(double* raw_data,  double* period, 
	double sampling_frequency, int n, double dj, double s0, int J, double maximum_frequency,
	double* result)
{
	
	//Variable Declarations
	int i, j;
	fftw_complex *data_in, *fft_data;
	fftw_plan plan_forward;

	int start = (int) floor( log2( 1.0/(s0 * maximum_frequency * FOURIER_WAVELENGTH_FACTOR) ) /dj);

	//Calculate Padding Required
	const int pad = floor(log(n)/log(2.0) + 0.499);
    const double PADDED_SIZE = pow(2, pad + 1);

    const double dw = (2 * M_PI * sampling_frequency)/(PADDED_SIZE); //NOT IN RAD/SEC in Hz

    data_in  = (fftw_complex *) fftw_malloc( sizeof( fftw_complex )*PADDED_SIZE );
	fft_data = (fftw_complex *) fftw_malloc( sizeof( fftw_complex )*PADDED_SIZE );

	//populate the data vector. 
	for (i = 0; i < n; ++i)
    {
    	data_in[i][0] = raw_data[i];
    	data_in[i][1] = 0.0;
    }

	//Calculate the FFT for the Data
	plan_forward = fftw_plan_dft_1d(PADDED_SIZE, data_in, fft_data, 
		FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan_forward);

	// // Heaviside Step Function whenever w0 < 0 the function is zero
	// for (int j = PADDED_SIZE/2; j < PADDED_SIZE; ++j)
	// {
	// 	fft_data[j][0] = 0.0;
	// 	fft_data[j][1] = 0.0;
	// }
		
	double value;
	
	fftw_plan plan_backward;
	fftw_complex *filter_convolution, *fftw_result;

	filter_convolution = (fftw_complex *) fftw_malloc( sizeof( fftw_complex )*PADDED_SIZE );
	fftw_result  = (fftw_complex *) fftw_malloc( sizeof( fftw_complex )*PADDED_SIZE );
	
		//Preapre for the plan backwards
		plan_backward = fftw_plan_dft_1d(PADDED_SIZE, filter_convolution, fftw_result, 
			FFTW_BACKWARD, FFTW_ESTIMATE);	
    
	for (i = start; i < J; ++i)
	{
		//Calculate the scale and corrosponding frequency at the specific Scale
		double scale = s0 * pow(2, i * dj);
		period[i] = (W_0)/(scale * 2 * M_PI);

		//Normalization Factor needes to be recomputed at every scale.
		double normal = sqrt(2 * M_PI * scale * sampling_frequency);

		//Caluclate the Fourier Morlet at the specific scale. 
		for (j = 0; j < PADDED_SIZE; ++j)
		{
			value = FourierMorlet(j*dw, scale, normal);
			filter_convolution[j][0] = fft_data[j][0] * value;
			filter_convolution[j][1] = fft_data[j][1] * value;
		}

		//Take the inverse FFT. 
		fftw_execute(plan_backward);

		//Calculate the power and store it in result
		for (j = 0; j < n; ++j)
		{
			result[i * n + j] = Magnitude(fftw_result[j][0], fftw_result[j][1]);
		}
	}

	//FFTW sanitation engineering. 
	fftw_destroy_plan(plan_backward);
    fftw_free(fftw_result);
    fftw_free(filter_convolution);

	fftw_destroy_plan(plan_forward); 
	fftw_free(fft_data); fftw_free(data_in);  
    return(0);
} /*Wavelet */

double FourierMorlet(double w, double scale, double normal)
{
	// w = w/scale;
	// const double w2 = w * w;
	////These are all needed by Fourier Morlet i'm going to move 
	////them out to optimize the code
	// const double k = exp(-0.5 * W_0_2);
	// const double normal = sqrt((2 * M_PI * scale)/(1.0/sampling_frequency));
	// const double cSigma = pow(1.0 + exp(-W_0_2) - 2*exp(-0.75*W_0_2), -0.5);
	// double out = exp( -0.5 * (W_0 - w) * (W_0 - w)) - k * exp(-0.5 * w * w);
	// double out = exp( -0.5 * (W_0_2 - 2*W_0*w + w2)) - k * exp(-0.5 * w2);
	// out = cSigma * QUAD_ROOT_PI * out;

	double exponent = -0.5 * (scale * w - W_0) * (scale * w - W_0);
	double out = QUAD_ROOT_PI* normal * exp(exponent);
	return(out);
}

double CompleteFourierMorlet(double w, double scale)
{
	double w_0 = 6.0;
	double c_sigma = pow( (1.0 + exp(-w_0*w_0) - 2.0 * exp(-0.75 * w_0 * w_0)), -0.5 );
	double k_sigma = exp( -0.5 * w_0 * w_0 );

	double out = exp( -0.5 * ( w_0 - scale * w ) * (w_0 - scale * w) ) 
					- k_sigma * (exp ( -0.5 * scale * w * w));
	out = c_sigma * QUAD_ROOT_PI * out;
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
	int t = 2 * FS;
	switch(flag)
	{
		//Impulse at T = 2 seconds
		case 1:
			
			for (int i = 0; i < DATA_SIZE; ++i)
			{
				data[i] = 0.0;
			}

			data[t] = 1.0;
			break;
		
		//Multiple Sines at t = 1500
		case 2:
			for (int i = DATA_SIZE/2; i < DATA_SIZE/2 + 2*one_peri; ++i)
			{
				data[i] = sin((i - DATA_SIZE/2)* dw + w0) + sin((i - DATA_SIZE/2)* 2* dw + w0);
			}
			break;
		//Multiple Sines at all times
		case 3:
			for (int i = 0; i < DATA_SIZE; ++i)
			{
				data[i] = sin(i*dw + w0) + sin(i*2*dw + w0);
			}
			break;
		//Single sine at t = 1.0s;
		case 4: 
			for (int i = DATA_SIZE/3; i < DATA_SIZE/3 + 2 * one_peri; ++i)
			{
				data[i] = sin( (i - DATA_SIZE/2) * dw + w0 );
			}
			break;
		//Single sine all the way through. 
		case 5:
			for (int i = 0; i < DATA_SIZE; ++i)
			{
				data[i] = cos( i * dw + w0 );
				if (i >= DATA_SIZE/3 && i <= DATA_SIZE/2)
				{
					data[i] = 2 * cos(i * dw + w0);
				}
			}
			break;
		case 6:
			for (int i = 0; i < DATA_SIZE; ++i)
			{
				data[i] = cos(i * dw + w0 );
				if (i >= DATA_SIZE/3 && i <= DATA_SIZE/2)
				{
					data[i] = cos(i * (dw - 0.005) + w0);
				}
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

int WriteFile(double *data, double *period, int x, int y, char filename[])
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
    	fprintf(out_file, "%f\t", period[i]);

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
    	// fprintf(out_file, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", i,
    	// 	data[length*0 + i] + 0., data[length*1 + i] + 5., data[length*2 + i] + 10, 
    	// 	data[length*3 + i] + 15, data[length*4 + i] + 20, data[length*5 + i] + 25, 
    	// 	data[length*6 + i] + 30, data[length*7 + i] + 35, data[length*8 + i] + 40, 
    	// 	data[length*9 + i] + 45, data[length*10+ i] + 50);
    	fprintf(out_file, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", i,
    		data[length*40 + i] + 0., data[length*41 + i] + 5., data[length*42 + i] + 10, 
    		data[length*43 + i] + 15, data[length*44 + i] + 20, data[length*45 + i] + 25, 
    		data[length*46 + i] + 30, data[length*47 + i] + 35, data[length*48 + i] + 40, 
    		data[length*49 + i] + 45, data[length*50+ i] + 50);
    }
    
    fclose(out_file);
    return 0;
}
