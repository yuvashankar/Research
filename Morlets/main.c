#include "Morlet.h"

int main(void)
{
    //Memory Allocation
    double* filter = malloc(DATA_SIZE * sizeof(double));

    // double raw_data[MAX_DATA_SIZE];
    double * raw_data = malloc(MAX_DATA_SIZE * sizeof(double));
    char fileName[] = "sst_nino3.dat";
    FILE* out_file = fopen("DATA.log", "w");
    

    // int n = ReadFile(raw_data, fileName);
    fillData(raw_data);
    int n = DATA_SIZE;
    
    double mean = gsl_stats_mean(raw_data, 1, n);
    double deviation = gsl_stats_variance(raw_data, 1, n);

    //Normalizing with the standard deviation because we're getting some
    //crazy numbers on the other end. 
    for (int i = 0; i < n; ++i)
    {
        raw_data[i] = (raw_data[i] - mean)/deviation;
    }

    //Now we're done with GSL, move on to FFTW
    //FFTW allocations
    //Find the closest 2 power of DATA_SIZE, and pad the arrays to the left with that. 
    int pad = floor(log(n)/log(2.0) + 0.4999);
    double PADDED_SIZE = pow(2, pad + 1);
    fftw_complex *data_in, *fft_data, *result;
    fftw_plan plan_forward, plan_backwards;

    data_in = fftw_alloc_complex(PADDED_SIZE); 
    fft_data = fftw_alloc_complex(PADDED_SIZE);
    result = fftw_alloc_complex(PADDED_SIZE);
    
    //Copy memory
    for (int i = 0; i < n; ++i)
    {
        data_in[i][0] = raw_data[i];
    }

    //Plan and calculate the Fourier Transform of the data, the fft is stored in fft_data
    plan_forward = fftw_plan_dft_1d(PADDED_SIZE, data_in, fft_data, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan_forward);

    CreateComplexFilter(filter);

    for (int i = 0; i < DATA_SIZE; ++i)
    {
        fft_data[i][0] *= filter[i];
    }

    plan_backwards = fftw_plan_dft_1d(PADDED_SIZE, fft_data, result, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan_backwards);

    for (int i = 0; i < DATA_SIZE; ++i)
    {
        // fprintf(out_file, "%i\t%f\t%f\n", i, data_in[i][0], result[i][0]);
        // fprintf(out_file, "%i\t%f\n", i, data_in[i][0]);
        fprintf(out_file, "%i\t%f\n", i, raw_data[i]);
    }

    fclose(out_file);

    free(raw_data); //Freeing the data because it's HUGE!
    free(filter);
    //FFTW sanitation. 
    fftw_destroy_plan(plan_forward); fftw_destroy_plan(plan_backwards);
    fftw_free(data_in); fftw_free(fft_data); fftw_free(result);

    return 0;
}