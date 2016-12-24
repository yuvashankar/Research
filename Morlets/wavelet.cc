//Wavelet.c
#include "wavelet.h"
#include <omp.h>


#define TEST 0.00001

int Wavelet(double* raw_data,  double* period, double* scales, 
	double sampling_frequency, int n, int J,
	double* result)
{
	//Variable Declarations
	int i, j;
	fftw_complex *data_in, *fft_data;
	fftw_plan plan_forward;

	//Calculate Padding Required
    const int PADDED_SIZE = CalculatePaddingSize(n, 1);

    const double dw = (2 * M_PI * sampling_frequency)/(PADDED_SIZE); //NOT IN RAD/SEC in Hz

    data_in  = (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * PADDED_SIZE );
	fft_data = (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * PADDED_SIZE );

	//populate the FFTW data vector. 
	for (i = 0; i < n; ++i)
    {
    	data_in[i][0] = raw_data[i];
    	data_in[i][1] = 0.0;
    }

    //Force the rest of the data vector to zero just in case
    for (int i = n; i < PADDED_SIZE; ++i)
    {
    	data_in[i][0] = 0.0;
    	data_in[i][1] = 0.0;
    }

	//Calculate the FFT of the data and store it in fft_data
	plan_forward = fftw_plan_dft_1d(PADDED_SIZE, data_in, fft_data, 
									FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan_forward);

	// FILE* out_file = fopen("debug.log", "w");

	#pragma omp parallel num_threads(1) private(i, j) shared (result, period, sampling_frequency, J, n, scales,  fft_data) default(none)
	{
		double value;

		fftw_plan plan_backward;
		fftw_complex *filter_convolution, *fftw_result;

		filter_convolution = (fftw_complex *) fftw_malloc( sizeof( fftw_complex )* PADDED_SIZE );
		fftw_result  = 		 (fftw_complex *) fftw_malloc( sizeof( fftw_complex )* PADDED_SIZE );

		#pragma omp critical (make_plan)
		{
			//Preapre for the plan backwards
			plan_backward = fftw_plan_dft_1d(PADDED_SIZE, filter_convolution, fftw_result, 
				FFTW_BACKWARD, FFTW_ESTIMATE);
		}
		
	    #pragma omp for
		for (i = 0; i < J; ++i)
		{
			//Calculate the corrosponding frequency to the scale
			period[i] = (W_0)/(scales[i] * 2 * M_PI);

			//Compute the Fourier Morlet at 0 and N/2
			value = CompleteFourierMorlet(0.0, scales[i]);
		    

			filter_convolution[0][0] = fft_data[0][0] * value;
			filter_convolution[0][1] = fft_data[0][1] * value;
			
			filter_convolution[PADDED_SIZE/2][0] = 0.0;
			filter_convolution[PADDED_SIZE/2][1] = 0.0;
			//Compute the Fourier Morlet Convolution in between
			for (j = 1; j < PADDED_SIZE/2 - 1; ++j)
			{
				value = CompleteFourierMorlet( j * dw , scales[i]);
				filter_convolution[j][0] = fft_data[j][0] * value;
				filter_convolution[j][1] = fft_data[j][1] * value;

				filter_convolution[PADDED_SIZE- j][0] = 0.0;
				filter_convolution[PADDED_SIZE- j][1] = 0.0;
			}

			//Take the inverse FFT. 
			fftw_execute(plan_backward);
		    
			//Calculate the power and store it in result
			for (j = 0; j < n; ++j)
			{
				result[i * n + j] = MAGNITUDE(fftw_result[j][0], fftw_result[j][1]);
			}
		}

		//FFTW sanitation engineering. 
		fftw_destroy_plan(plan_backward);
	    fftw_free(fftw_result);
	    fftw_free(filter_convolution);
	}
	//Sanitation Engineering
	fftw_destroy_plan(plan_forward); 
	fftw_free(fft_data); fftw_free(data_in);
    return(0);
} /*Wavelet */

int CalculatePaddingSize(int array_size, int FLAG)
{
	const int pad = ceil(log2(array_size));
	int out = array_size;
	switch(FLAG)
	{
		case 0: //No Padding what so ever. 
			out = array_size;
			break;
		case 1: //Zero - Padding and preforming a Radix-2 Operation
			out = (int) pow(2, pad + 1);
    		break;
    	case 2: //Duplicate array and ramp up and ramp down output
    		out = 2 * array_size;
    		break;
    	default: //Else return the array size
    		out = array_size;
    		break;
	}
	return(out);
}

double* GenerateScales(double minimum_frequency, double maximum_frequency)
{
	int min_i = FREQ_TO_SCALE(maximum_frequency);
	int max_i = FREQ_TO_SCALE(minimum_frequency) + 1;

	double act_max_freq = SCALE_TO_FREQ(min_i);
	double act_min_freq = SCALE_TO_FREQ(max_i);
	printf("Max Frequency = %f, Min Frequency = %f\n", act_max_freq, act_min_freq);

	double * scales = (double*) malloc ( (max_i - min_i) * sizeof(double) );
	int count = ( max_i - min_i );

	//Populate the scales array
	for (int i = 0; i < count; ++i)
	{
		int counterVariable = min_i + i;
		scales[i] = S0 * pow(2, counterVariable * D_J);
	}

	return(scales);
}

double FourierMorlet(double w, double scale, double normal)
{
	double exponent = -0.5 * (scale * w - W_0) * (scale * w - W_0);
	double out = QUAD_ROOT_PI* normal * exp(exponent);
	return(out);
}

double CompleteFourierMorlet(double w, double scale)
{
	// double W_0 = 6.0;
	// double C_SIGMA = pow( (1.0 + exp(-W_0_2) - 2.0 * exp(-0.75 * W_0_2)), -0.5 );
	// double K_SIGMA = exp( -0.5 * W_0_2 );

	double norm = pow(scale, -0.5);
	double out = exp( -0.5 * ( W_0 - scale * w ) * (W_0 - scale * w) ) 
					- K_SIGMA * (exp ( -0.5 * scale * w * w));
	out = norm * C_SIGMA * QUAD_ROOT_PI * out;
	return(out);
}

void FillData(double * data)
{
	// Fit a FREQ signal at two points
	// double DT = 1./FS;
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
	// double DT = 1./FS;
	double fsig = FREQ/FS;
	double dw = 2*M_PI*fsig;
	double w0 =  0.01; // A SMALL PHASE SHIFT SO ITS NOT ALL INTERGER ALIGNED
	int one_peri = (int)1./fsig;
	printf("FS  %.2f   Pitch %.f   Discrete Priode = %d \n",FS,FREQ,one_peri);

	// double frequency = MIN_FREQUENCY;
	double frequency_increment = (MAX_FREQUENCY - MIN_FREQUENCY)/ 3.0; //3.0 seconds. 
	
	int t = 2 * FS; //At 2 seconds. 

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
				// if (i >= DATA_SIZE/2 && i <= 2 * (DATA_SIZE)/3)
				// {
				// 	data[i] = 2 * cos(i * dw + w0);
				// }
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
			break;
		//Frequency Sweep
		case 7:
			for (int i = 0; i < DATA_SIZE; ++i)
			{
				// fsig = frequency/FS;
				// dw = 2*M_PI*fsig;

				data[i] = sin(w0 + 2 * M_PI * (MIN_FREQUENCY + (frequency_increment/2) * pow(i/FS, 2)) );

				// frequency += frequency_increment; 
			}
			break;
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

	int result = fread (buffer, 1, lSize, signalFile);
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

int WriteFile(double *data, double *period, int x, int y, const char* filename)
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
