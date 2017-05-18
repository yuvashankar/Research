/**
	\file "Wavelet.cc"
	\brief This file contains all of the functions that support the ERSP and CWT functions
*/
#include "wavelet.h"
#include <omp.h>
#include <assert.h>
#include <gsl/gsl_statistics.h>

#define TEST 0.00001
#define SETTLING_PERCENTAGE 0.02
#define NORMALIZATION_FACTOR 0.375402913609157562

int Wavelet(double* raw_data, double* scales, 
	double sampling_frequency, int n, int J,
	double* result)
{
	//Variable Declarations
	int i, j;
	fftw_complex *data_in, *fft_data;
	fftw_plan plan_forward;


	//Calculate Padding Required
    // const int PADDED_SIZE = CalculatePaddingSize(n, 1);
    const int PADDED_SIZE = n;

    const double dw = (2 * M_PI * sampling_frequency)/(PADDED_SIZE); //NOT IN RAD/SAMPLE in RAD/SEC

    data_in  = (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * PADDED_SIZE );
	fft_data = (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * PADDED_SIZE );


	//populate the FFTW data vector
	for (i = 0; i < n; ++i)
    {
    	data_in[i][0] = raw_data[i];
    	// data_in[i + n][0] = raw_data[i];

    	// data_in[i + n][1] = 0.0;
    	data_in[i    ][1] = 0.0;
    }

    // //Force the rest of the data vector to zero just in case
    // for (int i = n; i < PADDED_SIZE; ++i)
    // {
    // 	data_in[i][0] = 0.0;
    // 	data_in[i][1] = 0.0;
    // }

    double *temp = (double*) malloc(n * sizeof(double));
    FILE* debug_file = fopen("debug.log", "w");

	//Calculate the FFT of the data and store it in fft_data
	plan_forward = fftw_plan_dft_1d(PADDED_SIZE, data_in, fft_data, 
									FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan_forward);

	// #pragma omp parallel num_threads(1) private(i, j) shared (result, sampling_frequency, J, n, scales,  fft_data) default(none)
	// {
		double value;

		fftw_plan plan_backward;
		fftw_complex *filter_convolution, *fftw_result;

		filter_convolution = (fftw_complex *) fftw_malloc( sizeof( fftw_complex )* PADDED_SIZE );
		fftw_result  = 		 (fftw_complex *) fftw_malloc( sizeof( fftw_complex )* PADDED_SIZE );

		// #pragma omp critical (make_plan)
		// {
			//Preapre for the plan backwards
			plan_backward = fftw_plan_dft_1d(PADDED_SIZE, filter_convolution, fftw_result, 
				FFTW_BACKWARD, FFTW_ESTIMATE);
		// }
		
	    // #pragma omp for
		for (i = 0; i < J; ++i)
		{
			//Force the arrays to zero
			memset(filter_convolution, 0.0, sizeof( fftw_complex ) * PADDED_SIZE);
			memset(fftw_result,        0.0, sizeof( fftw_complex ) * PADDED_SIZE);

			//Compute the Fourier Morlet at 0 and N/2
			double norm = sqrt(scales[i]);

			value = CompleteFourierMorlet(0.0, scales[i], norm);

			filter_convolution[0][0] = ( fft_data[0][0] / PADDED_SIZE ) * value;
			filter_convolution[0][1] = ( fft_data[0][1] / PADDED_SIZE ) * value;
			
			filter_convolution[PADDED_SIZE/2][0] = 0.0;
			filter_convolution[PADDED_SIZE/2][1] = 0.0;

			//Compute the Fourier Morlet Convolution in between
			for (j = 1; j < PADDED_SIZE/2 - 1; ++j)
			{
				value = CompleteFourierMorlet( j * dw , scales[i], norm);

				filter_convolution[j][0] = ( fft_data[j][0] / PADDED_SIZE ) * value;
				filter_convolution[j][1] = ( fft_data[j][1] / PADDED_SIZE ) * value;

				filter_convolution[PADDED_SIZE- j][0] = 0.0;
				filter_convolution[PADDED_SIZE- j][1] = 0.0;
			}

			//Take the inverse FFT. 
			fftw_execute(plan_backward);
		    
			//Calculate the power and store it in result
			for (j = 0; j < n; ++j)
			{
				result[i * n + j] = MAGNITUDE(fftw_result[j][0], fftw_result[j][1]) / ( NORMALIZATION_FACTOR * sqrt(scales[i]) );
				temp[j] = result[i * n + j];
			}
			
			double total = 0;
			for (int j = 0; j < n; ++j)
			{
				total += temp[j];
			}

			double variance = gsl_stats_variance(temp, 1, n);
			fprintf(debug_file, "%.16f\t%.16f\t%.16f\n", SCALE_TO_FREQ(scales[i]), variance, total);

		}

		free(temp);
		fclose(debug_file);

		//FFTW sanitation engineering. 
		fftw_destroy_plan(plan_backward);
	    fftw_free(fftw_result);
	    fftw_free(filter_convolution);
	// }
	//Sanitation Engineering
	fftw_destroy_plan(plan_forward); 
	fftw_free(fft_data); fftw_free(data_in);
    return(0);
} /*Wavelet */


double* ShortTimeFourierTransform(double * raw_data, double sampling_frequency, int n, int window_size)
{
	//Variable Declarations
	int i, j;
	fftw_complex *data_in, *fft_data;
	fftw_plan plan_forward;
	double* result;

	// FILE * STFT_FILE = fopen("STFT_Result.log", "w");

	int num_windows = ceil((double) n / window_size); 
	assert(num_windows > 0);

	// const int PADDED_SIZE = CalculatePaddingSize(window_size, 1);
	const int PADDED_SIZE = window_size;


	result =   (double*)        malloc(num_windows * PADDED_SIZE * sizeof(double));
	data_in  = (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * PADDED_SIZE );
	fft_data = (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * PADDED_SIZE );

	//Calculate the FFT of the data and store it in fft_data
	plan_forward = fftw_plan_dft_1d(PADDED_SIZE, data_in, fft_data, 
									FFTW_FORWARD, FFTW_ESTIMATE);

	for (i = 0; i < num_windows; ++i)
	{
		memset(data_in,  0.0, sizeof( fftw_complex ) * PADDED_SIZE);
		memset(fft_data, 0.0, sizeof( fftw_complex ) * PADDED_SIZE);

		//Fill Data into FFT Array
		for (j = 0; j < window_size; ++j)
		{
			if (i * window_size + j < n)
			{
				data_in[j][0] = raw_data[i * window_size + j];
				data_in[j][1] = 0.0;
			}
			else
			{
				data_in[j][0] = 0.0;
				data_in[j][1] = 0.0;
			}
		}
		
		fftw_execute(plan_forward);

		for (j = 0; j < window_size / 2; ++j)
		{
			result[j * num_windows + i] = MAGNITUDE(fft_data[j][0], fft_data[j][1]);
			// fprintf(STFT_FILE, "%d\t%f\t%d\n", j, result[j * num_windows + i], j * num_windows + i);
		}

	}

	// fclose(STFT_FILE);
	fftw_destroy_plan(plan_forward); 
	fftw_free(data_in); fftw_free(fft_data);

	return(result);
}


int Find_Peaks(double* array, double* frequency, int n, int J)
{
	FILE*       maximum_file  = fopen("maximum.log", "w");
	ARRAY_DATA *maximum_array = (ARRAY_DATA*) malloc (J * sizeof(ARRAY_DATA));
	int local_maximum_location[J];
	double* temp = (double*) malloc(n * sizeof(double));

	//Find the local maximum at every frequency
	for (int i = 0; i < J; ++i)
	{
		ARRAY_DATA max;
		max.value = array[i * n];
		max.index = i * n;

		for (int j = 0; j < n; ++j)
		{
			if (array[i * n + j] > max.value)
			{
				max.value = array[i * n + j];
				max.index = i * n + j;
			}
		}

		maximum_array[i] = max;
		fprintf(maximum_file, "%f\t%f\n", frequency[i], maximum_array[i].value);
	}

	//Calculate the deravitive of the signal and isolate the peaks
	double sign = (maximum_array[1].value - maximum_array[0].value) / (frequency[1] - frequency[0]);
	int    max_count = 0;
	for (int i = 0; i < J - 1; ++i)
	{
		double slope = (maximum_array[i + 1].value - maximum_array[i].value) / (frequency[i + 1] - frequency[i]);
		if (signbit(slope) != signbit(sign) && sign < 0)
		{
			local_maximum_location[max_count] = i;
			max_count++;
		}

		sign = slope;
	}

	
	for (int i = 0; i < max_count; ++i)
	{
		int arr_index = local_maximum_location[i];

		//Copy data into memory block
		for (int j = 0; j < n; ++j)
		{
			temp[j] = array[arr_index * n + j];
			// if (i == 1)
			// 	fprintf(maximum_file, "%f\t%.16f\n", (double) j/FS, array[arr_index * n + j]);
			
		}

		ARRAY_DATA impact_site = Max(temp, n);

		int system_setteled = 0;
		int setteled_index = 0;
		for (int j = impact_site.index; j < n; ++j)
		{
			if (temp[j] < SETTLING_PERCENTAGE * impact_site.value && system_setteled == 0)
			{
				setteled_index = j;
				double setteled_time = (double) (setteled_index - impact_site.index)/FS;
				printf("Frequency[%d]: %f, Settled Time = %f\n", i, frequency[arr_index], setteled_time);
				system_setteled = 1;
			}
		}
	}


	fclose(maximum_file);
	free(maximum_array);
	free(temp);
	return(0);
}



int CalculatePaddingSize(const int array_size, const int pad_flag)
{
	const int pad = ceil(log2(array_size));
	int out = array_size;
	switch(pad_flag)
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

    	case 3: //Repeat the array once. 
    		out = 2 * array_size;
    		break;

    	default: //Else return the array size
    		out = array_size;
    		break;
	}
	return(out);
}

void GenerateScalesAndFrequency(const int min_i, const int max_i, const double s_0, 
	double* scales, double* frequency)
{
	int count = ( max_i - min_i ) + 1;

	//Populate the scales array
	for (int i = 0; i < count; ++i)
	{
		int counterVariable = min_i + i;
		scales[i] = s_0 * pow(2, counterVariable * D_J);
	}

	IdentifyFrequencies(scales, count, frequency);

}

void IdentifyFrequencies(double* scales, int count, double* frequency)
{
	// double * frequency = (double*) malloc( count * sizeof(double) );

	for (int i = 0; i < count; ++i)
	{
		frequency[i] = (W_0)/(scales[i] * 2 * M_PI);
	}

}

double CompleteFourierMorlet(double w, const double scale, double norm)
{
	// double norm = 1.0/sqrt(scale);
	// double norm = sqrt(scale);
	
	w = w * scale; 
	double out = exp( -0.5 * ( W_0 - w ) * (W_0 - w) )
                                       - K_SIGMA * (exp ( -0.5 * w * w));

	out = norm * C_SIGMA * QUAD_ROOT_PI * out;
	return(out);
}

void TestCases(double *data, const int flag, double freq, double sampling_frequency, int data_size)
{
	// Fit a freq signal at two points
	// double DT = 1./sampling_frequency;
	double fsig = freq/sampling_frequency;
	double dw = 2 * M_PI * fsig;
	double w0 =  0.001; // A SMALL PHASE SHIFT SO ITS NOT ALL INTERGER ALIGNED
	int one_peri = (int)1./fsig;

	double frequency_increment = (MAX_FREQUENCY - MIN_FREQUENCY)/ 3.0; //3.0 seconds. 
	
	int t = 2 * sampling_frequency; //At 2 seconds. 

	switch(flag)
	{
		//Impulse at T = 2 seconds
		case 1:		
			for (int i = 0; i < data_size; ++i)
			{
				data[i] = 0.0;
			}

			data[t] = 1.0;
			break;
		
		//Multiple Sines at t = 1500
		case 2:
			for (int i = data_size/2; i < data_size/2 + 2*one_peri; ++i)
			{
				data[i] = sin((i - data_size/2)* dw + w0) + sin((i - data_size/2)* 2* dw + w0);
			}
			break;

		//Multiple Sines at all times
		case 3:
			for (int i = 0; i < data_size; ++i)
			{
				data[i] = sin(i*dw + w0) + sin(i*2*dw + w0);
			}
			break;

		//Single sine at t = 1.0s;
		case 4: 
			for (int i = data_size/3; i < data_size/3 + 2 * one_peri; ++i)
			{
				data[i] = sin( (i - data_size/2) * dw + w0 );
			}
			break;

		
		case 5:
			for (int i = 0; i < data_size; ++i)
			{
				data[i] = cos( i * dw + w0 );
				
				if (i >= data_size/2 && i <= 2 * (data_size)/3)
				{
					data[i] = 2 * cos(i * dw + w0);
				}
			}
			break;

		case 6:
			for (int i = 0; i < data_size; ++i)
			{
				data[i] = cos(i * dw + w0 );
				if (i >= data_size/3 && i <= data_size/2)
				{
					data[i] = cos(i * (dw - 0.005) + w0);
				}
			}
			break;

		//Frequency Sweep
		case 7:
			for (int i = 0; i < data_size; ++i)
			{
				// fsig = frequency/sampling_frequency;
				// dw = 2*M_PI*fsig;

				data[i] = sin(w0 + 2 * M_PI * (MIN_FREQUENCY + (frequency_increment/2) * pow(i/sampling_frequency, 2)) );

				// frequency += frequency_increment; 
			}
			break;
			
		//Single sine all the way through. 
		case 8:
			for (int i = 0; i < data_size; ++i)
			{
				data[i] = cos(i*dw + w0);
			}
			break;
	}
}

int WriteDebug(const double *data, const int length, const int sampling_frequency,
	const char* filename)
{
	FILE* out_file=fopen(filename,"w");
    if (out_file == NULL) return -1;

    double t = 0.0;
    double dt = 1.0/sampling_frequency;

	for (int i = 0; i < length; ++i)
    {
    	// double value = (double) i/length;
    	fprintf(out_file, "%f\t%.16e\n", t, data[i]);
    	t += dt;
    }
    
    fclose(out_file);
    return 0;
}

void FillData(double * data)
{
	// Fit a FREQ signal at two points
	// double DT = 1./FS;
	double fsig = FREQ/FS;
	double dw = 2*M_PI*fsig;
	double w0 =  0.01; // A SMALL PHASE SHIFT SO ITS NOT ALL INTERGER ALIGNED
	int one_peri = (int)1./fsig;
	// printf("FS  %.2d   Pitch %.f   Discrete Period = %d \n",FS,FREQ,one_peri);

	for (int i = 0; i < DATA_SIZE; ++i)
	{
		data[i] = 0.0;
	}

	// //Impulse Sample
	// data[2000] = 1.0;
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

int GetFileSize(char filename[])
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

	char * token = strtok(buffer, "\n");
	
    //Get input from text.
	int counterVariable = 0;
	while (token !=NULL)
    {
    	// data[counterVariable] = atof(token);
    	counterVariable++;
        token = strtok (NULL, "\n");

    }
    fclose(signalFile);

    return (counterVariable);
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


double CWT_Cosine_Real(double time, double scale)
{
    double norm = sqrt(scale);
    double w_o = FREQ * 2 * M_PI;
    w_o *= scale;
    double mor = exp( - 0.5 * (W_0 - w_o) * (W_0 - w_o) ) + K_SIGMA * exp( -0.5 * w_o * w_o );
    double cosine = cos(w_o * time);

    double out = 0.5 * C_SIGMA * QUAD_ROOT_PI * norm * mor * cosine;
    return(out);
}

double CWT_Cosine_Complex(double time, double scale)
{
    double norm = sqrt(scale);
    double w_o = FREQ * 2 * M_PI;
    w_o *= scale;
    double mor = exp( - 0.5 * (W_0 - w_o) * (W_0 - w_o) ) + K_SIGMA * exp( -0.5 * w_o * w_o );
    double sine = sin(w_o * time);

    double out = 0.5 * C_SIGMA * QUAD_ROOT_PI * norm * mor * sine;
    return(out);
}

double CWT_Dirac_Real(double time, double scale)
{
    double impulse = 2.0;
    time = (impulse - time )/scale;
    double norm = 1.0/sqrt(scale);
    double out = exp( - 0.5 * time * time ) * ( cos( W_0 * time ) - K_SIGMA );

    out = norm * C_SIGMA * QUAD_ROOT_PI * out;
    return(out);
}

double CWT_Dirac_Complex(double time, double scale)
{
    double impulse = 2.0;
    time = (impulse - time )/scale;
    double norm = 1.0/sqrt(scale);
    double out = exp( - 0.5 * time * time ) * ( sin( W_0 * time ) - K_SIGMA );
    out = norm * C_SIGMA * QUAD_ROOT_PI * out;
    return(out);
}

ARRAY_DATA Max(double * array, int size)
{
	double max = array[0];
	int array_index = 0;

	ARRAY_DATA out;

	for (int i = 0; i < size; ++i)
	{
		if (array[i] > max)
		{
			max = array[i];
			// if (max != max)
			// {
			// 	printf("Naan Alert!\n");
			// }
			array_index = i;
		}
	}

	out.value = max;
	out.index = array_index;

	// printf("Max: Array[%d] = %.17f\n", array_index, array[array_index]);
	return(out);
}

/**
    \fn int OpenFile(const char* fileName, struct edf_hdr_struct *header)
    
    \brief Openes a .BDF file and allocates it to an edf_hdr_struct.

    \param fileName The name and location of the file to be opened
    \param header The pointer to the edf header structure

    \return 0 if file is opened successfully
    \return -1 if there is an error

*/
int OpenFile(const char* fileName, struct edf_hdr_struct *header)
{
    if(edfopen_file_readonly(fileName, header, EDFLIB_READ_ALL_ANNOTATIONS))
    {
        switch(header->filetype)
        {
            case EDFLIB_MALLOC_ERROR                : printf("\nmalloc error\n\n");
                                                    break;
            case EDFLIB_NO_SUCH_FILE_OR_DIRECTORY   : printf("\ncan not open file, no such file or directory\n\n");
                                                    break;
            case EDFLIB_FILE_CONTAINS_FORMAT_ERRORS : printf("\nthe file is not EDF(+) or BDF(+) compliant\n"
                                                           "(it contains format errors)\n\n");
                                                    break;
            case EDFLIB_MAXFILES_REACHED            : printf("\nto many files opened\n\n");
                                                    break;
            case EDFLIB_FILE_READ_ERROR             : printf("\na read error occurred\n\n");
                                                    break;
            case EDFLIB_FILE_ALREADY_OPENED         : printf("\nfile has already been opened\n\n");
                                                    break;
            default                                 : printf("\nunknown error\n\n");
                                                    break;
        }

    return(-1);
    }
    return(0);
}

/**
    \fn int64_t FindTriggers(const int * statusInput, const int64_t numberOfRecords,
                        int64_t * outputBuffer)

    \brief This function should take an array input and return the rising and falling edges of the triggers. 

    \param statusInput: The Status Channel Input from the BDF or EDF flie. use the edfread_digital_samples
    \param numberOfRecords: The size of statusInput
    \param outputBuffer: a 1 x 2 * MAXIMUM_TRIGGERS int64_t array with the odd entries being the
                rising edge and the even entries being the falling edges. 

    \return counterVariable The number of triggers that were found.                
*/


int FindTriggers(const int * statusInput, const int64_t numberOfRecords,
                        int64_t * outputBuffer)
{
    int64_t counterVariable = 0;
    int counterVariable_int = 0;
    //int i needs to be int64_t because we are recording it into a int64_t arrray. 
    int edge = 0;
    for (int64_t i = 0; i < numberOfRecords; ++i)
    {
        //Bit and the lower 16 bits to see if any of the triggers have been triggered to on. 
        if ( ((statusInput[i] & 0x0000FFFF) > 0) && (edge == 0) && (statusInput[i-1] != statusInput[i]) ) //Rising Edge Detected.
        {
            outputBuffer[counterVariable] = i;
            // printf("outputBuffer[%" PRId64 "] = %" PRId64"\n", counterVariable, outputBuffer[counterVariable]);
            
            counterVariable++;
            counterVariable_int++;
            edge = 1;
        }

        if (statusInput[i-1] != statusInput[i] && edge == 1) //Falling Edge Detected.
        {
            edge = 0;
        }

        if (counterVariable > MAXIMUM_TRIGGERS)
            return -1;
    }

    return(counterVariable_int);
}

/**
    \fn int FilterTriggers(const int code, const int button, const int numberOfRecords, 
            const int64_t * triggerList,
            const int * readBuffer, 
            int * outputBuffer)

    \brief Filteres the triggers coming in, and finds the specified events

    \param code The code of the trigger list that is needed 
                Possible Inputs: 1, or 2
    \param button The code of the button that is needed
                Possible Inputs 1, or 2
    \param triggerList The list of all of the possible triggers
    \param readBuffer The Status Channel input from the file
    
    \param outputBuffer The buffer that FilterTriggers will populate with the location of the location of the triggers that we're looking for

    \return counterVariable The number of triggers found.
*/

int FilterTriggers(const int code, 
    const int button, 
    const int numberOfRecords, 
    const int64_t * triggerList,
    const int * readBuffer, 
    int * outputBuffer)
{
    int readCode;
    int buttonCode;
    int counterVariable = 0;
    for (int i = 0; i < numberOfRecords; ++i)
    {
        readCode = readBuffer[i] & 255;
        buttonCode = (readBuffer[i] >> 8) & 3;

        if ( (readCode == code) && (buttonCode == button) )
        {
            outputBuffer[counterVariable] = triggerList[i];
            counterVariable++;
        }
    }
    return counterVariable;
}

void CleanData(double * data, double n)
{
    double mean = gsl_stats_mean(data, 1, n);
    double sDeviation = gsl_stats_sd_m(data, 1, n, mean);
    // printf("Mean: %f, SD: %f\n", mean, sDeviation);

    //Compute the Z-Score or Standard Score
    for (int i = 0; i < n; ++i)
    {
        data[i] = (data[i] - mean)/sDeviation;
    }
}

int  WriteFile(const double *data, const double *frequency, const int x, const int y, int sampling_frequency,
    const char filename[])
{

    FILE* out_file=fopen(filename,"w");
    if (out_file == NULL) return -1;

    //Xticks
    fprintf(out_file, "%d\t", x);
    for (int i = 0; i < y; ++i)
    {
        fprintf(out_file, "%f\t", (double) i/sampling_frequency);
    }
    fprintf(out_file, "\n");

    // double small_eps = 0.00001; //Add a small eps so that logs of zero don't happen. 
    for (int i = 0; i < x; ++i)
    {
        //Feed Frequency
        fprintf(out_file, "%f\t", frequency[i]);

        //Feed Data
        for (int j = 0; j < y; ++j)
        {
            // value = Magnitude(result[i*n + j], result[i*n + j]);
            fprintf(out_file, "%.16e\t", data[i*y + j]);
        }
        //Ready for the next line.
        fprintf(out_file, "\n");

    }

    fclose(out_file);
    return(0);
}

int WriteGnuplotScript(const char graph_title[], const char filename[])
{
    FILE* gnuplot_file = fopen("script.gplot", "w");
    if (gnuplot_file == NULL) return -1;
    
    // fprintf(gnuplot_file, "set term x11\n");
    
    // fprintf(gnuplot_file, "set logscale z 10\n");

    fprintf(gnuplot_file, "set term pngcairo enhanced font 'arial,12'\n");
    fprintf(gnuplot_file, "%s%s%s","set output '", graph_title,  ".png' \n");
    fprintf(gnuplot_file, "set pm3d map\n");
    fprintf(gnuplot_file, "set logscale y 2\n");
    fprintf(gnuplot_file, "set ticslevel 0\n");
    fprintf(gnuplot_file, "set xlabel \"time (s)\"\n");
    fprintf(gnuplot_file, "set ylabel \"Frequency (Hz)\"\n");

    //Input Graph Title
    fprintf(gnuplot_file, "%s%s%s\n", "set title \"", graph_title, "\"");
    //Plot filename
    // plot "DATA.log" matrix nonuniform with pm3d t ''
    fprintf(gnuplot_file, "%s%s%s\n", "splot \"", filename, "\" matrix nonuniform with pm3d t ''");
    // fprintf(gnuplot_file, "pause -1 \"Hit Return to continue\"");

    return(0);
}