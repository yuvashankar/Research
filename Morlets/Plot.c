//Plot.c
#include <pngwriter.h>
#include "wavelet.h"

double base( double val );
double interpolate( double val, double y0, double x0, double y1, double x1 );

void plot(double * data,int num_x,int num_y)
{
	int i, j, k;
	int y;

	int lines_size =2;

	printf("X= %d num Lines= %d \n",num_x,num_y);
	// get_normalizer(data,num_x,num_y);

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

double getR(double gray) 
{
    return base( gray - 0.5 );
}
double getG(double gray) 
{
    return base( gray );
}
double getB(double gray) 
{
	return base( gray + 0.5 );
}

double base(double val) 
{
	if ( val <= -0.75 ) return 0;
	else if ( val <= -0.25 ) return interpolate( val, 0.0, -0.75, 1.0, -0.25 );
	else if ( val <= 0.25 ) return 1.0;
	else if ( val <= 0.75 ) return interpolate( val, 1.0, 0.25, 0.0, 0.75 );
	else return 0.0;
}

double interpolate( double val, double y0, double x0, double y1, double x1 ) 
{
	return (val-x0)*(y1-y0)/(x1-x0) + y0;
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