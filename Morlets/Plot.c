//Plot.c
#include <pngwriter.h>
#include "wavelet.h"

void plot(double * data,int num_x,int num_y)
{
	int i, j, k;
	int y;

	int lines_size =2;

	printf("X= %d num Lines= %d \n",num_x,num_y);
	// get_normalizer(data,num_x,num_y);
	double minimum = min(data, num_x*num_y);
	double maximum = max(data, num_x*num_y);

	pngwriter png(num_x+ 2*PLOT_OX,lines_size*num_y+2*PLOT_OY,0,"test.png");

	for(i = 1; i <= num_x; i++) 
	{
		for(j=1; j<=num_y; j++)
		{
			for(k=0; k<lines_size; k++)
			{
				double value= scale_log( (double) (data[i+j*num_x]));
				
				png.plot(PLOT_OX + i, PLOT_OY + j*lines_size + k, 
					getR(value), getG(value), getB(value));
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


double max(double * array, int size)
{
	double max = 0.0;
	for (int i = 0; i < size; ++i)
	{
		if (array[i] > max)
		{
			max = array[i];
		}
	}
	return(max);
}

double min(double* array, int size)
{
	double min = array[0];
	for (int i = 0; i < size; ++i)
	{
		if (array[i] < min)
		{
			min = array[i];
		}
	}
	return(min);
}