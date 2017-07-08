#include "wavelet.h"
#include <assert.h>
#include <gsl/gsl_statistics.h>

#define SETTLING_PERCENTAGE 1.02

double Determine_Dampening_Ratio(ARRAY_DATA impact_location, ARRAY_DATA settled_location, 
								int sampling_frequency);

int Find_Peaks(double* array, double* frequency, int sampling_frequency, int n, int J)
{
	FILE* debug_file = fopen("debug.log", "w");
	ARRAY_DATA *maximum_array = (ARRAY_DATA*) malloc(J * sizeof(ARRAY_DATA));
	double *temp              = (double*)     malloc(n * sizeof(double)    );
	int local_maximum_location[J];
	

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
		fprintf(debug_file, "%f\t%f\n", frequency[i], maximum_array[i].value);
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
			if (i == 1)
			{
				// fprintf(debug_file, "%f\t%.16f\n", (double) j/sampling_frequency, array[arr_index * n + j] );
			}
		}

		ARRAY_DATA impact_site = Max(temp, n);

		double local_mean = gsl_stats_mean(temp, 1, impact_site.index);
		
		int system_setteled = 0;
		ARRAY_DATA setteled_site;
		for (int j = impact_site.index; j < n; ++j)
		{
			if (temp[j] < SETTLING_PERCENTAGE * local_mean && system_setteled == 0)
			{
				setteled_site.index = j;
				setteled_site.value = temp[j];

				double dampening_ratio = Determine_Dampening_Ratio(impact_site, setteled_site, sampling_frequency);
				double setteled_time = (double) (setteled_site.index - impact_site.index)/sampling_frequency;
				printf("Frequency[%d]: %f, Settled Time = %f, Dampening Ratio = %f\n", i, frequency[arr_index], setteled_time, dampening_ratio);
				system_setteled = 1;
			}
		}
	}
	
	fclose(debug_file);
	// fclose(maximum_file);
	free(maximum_array);
	free(temp);
	return(0);
}

double Determine_Dampening_Ratio(ARRAY_DATA impact_location, ARRAY_DATA settled_location, int sampling_frequency)
{
	double dampening_ratio = 0.0;

	double rise = (double) log2( settled_location.value ) - log2( impact_location.value );
	double run  = (double)     ( settled_location.index   -       impact_location.index ) / sampling_frequency;

	// printf("Rise = %f, Run = %f\n", rise, run);
	if (run > 0.0)
	{
		dampening_ratio = rise / run;
	}
	else
	{
		return (-1);
	}

	return(dampening_ratio);
}