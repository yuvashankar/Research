//Plot.c
#include "wavelet.h"
#include <pngwriter.h>


void Plot(double * data,int num_x,int num_y)
{
	int i, j, k;
	double maximum, minimum, range;
	// int y;

	int lines_size =2;

	// printf("X= %d num Lines= %d \n",num_x,num_y);
	CalculateLog( data, num_x * num_y );
	minimum = Min(data, num_x*num_y);
	maximum = Max(data, num_x*num_y);

	range = maximum - minimum;

	// printf("Maximum = %f, Minimum = %f after Log\n", maximum, minimum);


	pngwriter png(num_x+ 2*PLOT_OX,lines_size*num_y+2*PLOT_OY,0,"test.png");

	for ( i = 1; i <= num_x; ++i)
	{
		for ( j = 1; j <= num_y; ++j)
		{
			for ( k = 0; k < lines_size; ++k)
			{
				double value = (double) (data[i + j*num_x]);
				COLOUR c = GetColour(value, minimum, maximum);
				
				png.plot(PLOT_OX + i , PLOT_OY + j*lines_size + k, c.r, c.g, c.b);
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

COLOUR GetColour(double v,double vmin,double vmax)
{
   COLOUR c = {1.0,1.0,1.0}; // white
   double dv;

   if (v < vmin)
      v = vmin;
   if (v > vmax)
      v = vmax;
   dv = vmax - vmin;

   if (v < (vmin + 0.25 * dv)) {
      c.r = 0;
      c.g = 4 * (v - vmin) / dv;
   } else if (v < (vmin + 0.5 * dv)) {
      c.r = 0;
      c.b = 1 + 4 * (vmin + 0.25 * dv - v) / dv;
   } else if (v < (vmin + 0.75 * dv)) {
      c.r = 4 * (v - vmin - 0.5 * dv) / dv;
      c.b = 0;
   } else {
      c.g = 1 + 4 * (vmin + 0.75 * dv - v) / dv;
      c.b = 0;
   }

   return(c);
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
		//The number is the logarithm of the lowest possible number.
		else
		{
			array[i] = log10(2.220446049250313e-16);
		}
			

	}
}