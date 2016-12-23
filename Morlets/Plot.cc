//Plot.c
#include "wavelet.h"
#include <pngwriter.h>


void Plot(double * data,int num_x,int num_y)
{
	int i, j, k;
	// int y;

	int lines_size =2;

	// printf("X= %d num Lines= %d \n",num_x,num_y);
	CalculateLog( data, num_x * num_y );
	minimum = Min(data, num_x*num_y);
	maximum = Max(data, num_x*num_y);
	printf("Maximum = %f, Minimum = %f after Log\n", maximum, minimum);


	pngwriter png(num_x+ 2*PLOT_OX,lines_size*num_y+2*PLOT_OY,0,"test.png");

	for(i = 1; i <= num_x; i++) 
	{
		for(j = 1; j<=num_y; j++)
		{
			for(k = 0; k<lines_size; k++)
			{
				// double value= scale_log( (double) (data[i+j*num_x]));
				
				// png.plot(PLOT_OX + i, PLOT_OY + j*lines_size + k, 
					// getR(value), getG(value), getB(value));
			}
		}
		
		// ADD A TIME MARKER
		if(i%50==0)
		{
	       png.plot(PLOT_OX+i,PLOT_OY,
	       	1.,1.,1.);
		}
	}
	png.close();
}

void CalculateLog(double * array, int size)
{
	double val;
	for (int i = 0; i < size; ++i)
	{
		//If the number is greater than machine precision of zero
		if (array[i] > 2.220446049250313e-16)
		{
			val = log10(array[i]);
			array[i] = val;
		}
		//Te number is the logarithm of the lowest possible number.
		else
		{
			array[i] = log10(2.220446049250313e-16);
		}
			

	}
}


double Max(double * array, int size)
{
	double max = array[0];
	int array_index = 0;

	for (int i = 0; i < size; ++i)
	{
		if (array[i] > max)
		{
			max = array[i];
			array_index = i;
		}
	}

	printf("Max: Array[%d] = %.17f\n", array_index, array[array_index]);
	return(max);
}

double Min(double* array, int size)
{
	int array_index = 0;
	double min = array[0];
	for (int i = 0; i < size; ++i)
	{
		if (array[i] < min)
		{
			min = array[i];
			array_index = i;
		}
	}
	printf("Min: Array[%d] = %.17f\n", array_index, array[array_index]);
	return(min);
}