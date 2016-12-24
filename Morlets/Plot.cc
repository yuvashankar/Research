//Plot.c
//All the functions needed to plot the png. 
#include "wavelet.h"
#include <pngwriter.h>

void Plot(double * data,int num_x,int num_y)
{
	int i, j, k;
	
	int lines_size =10;
	int image_width  = num_x + 2*PLOT_OX;
	int image_height = lines_size*num_y+2*PLOT_OY;
	int height = lines_size*num_y;

	char* font_location = "../lib/VeraMono.ttf";
	int label_font_size = 20.0;
	int title_font_size = 25.0;

	char* x_label = "Time (s)";
	char* y_label = "Frequency (Hz)";
	char* graph_title = "Time Frequency Graph of an Impulse";
	char t_in_string[4];

	CalculateLog( data, num_x * num_y );

	RANGE r = GetRange(data, num_x*num_y);
	printf("New Max = %f, New Min = %f\n", r.maximum, r.minimum);

	pngwriter png( image_width, image_height , 0 , "test.png");

	//X_label
	int x_text_width = png.get_text_width(font_location, label_font_size, x_label);
	png.plot_text( font_location, label_font_size,
				(0.5*image_width - 0.5*x_text_width), 0.5*PLOT_OY, 0.0,
				x_label,
				1.0, 1.0, 1.0);
	
	//Y_Label
	int y_text_width = png.get_text_width(font_location, label_font_size, x_label);
	png.plot_text( font_location, label_font_size,
				0.5*PLOT_OY, (0.5*image_height - 0.5*y_text_width), M_PI*0.5,
				y_label,
				1.0, 1.0, 1.0);

	//Title
	int title_text_width = png.get_text_width(font_location, label_font_size, x_label);
	png.plot_text( font_location, title_font_size,
				(0.5*image_width - title_text_width), (image_height - 0.5*PLOT_OY), 0.0,
				graph_title,
				1.0, 1.0, 1.0);

	// int counterVariable = 0;
	for ( i = 1; i <= num_x; ++i)
	{
		for ( j = 1; j <= num_y; ++j)
		{
			for ( k = 0; k < lines_size; ++k)
			{
				double value = data[i + j*num_x];
				COLOUR c = GetColour(value, r.minimum, r.maximum);
				
				png.plot(PLOT_OX + i , PLOT_OY + (height - j*lines_size + k), 
					c.r, c.g, c.b);
			}
		}
		// ADD A TIME MARKER
		if( i % int(FS) == 0 )
		{
			png.filledsquare( PLOT_OX + i - 4, PLOT_OY,
				PLOT_OY + i + 4, PLOT_OY - 20,
				1.0, 1.0, 1.0);

			sprintf(t_in_string, "%.1f", i/FS);
			int time_text_width = png.get_text_width(font_location, label_font_size, t_in_string);
			png.plot_text( font_location, label_font_size,
				PLOT_OY + i - 4 - 0.5* time_text_width, PLOT_OY - 50, 0.0,
				t_in_string,
				1.0, 1.0, 1.0);

		}

	}
	png.close();
}


RANGE GetRange(double* array, int size)
{
	RANGE r = {array[0], array[0]};


	for (int i = 0; i < size; ++i)
	{
		if (array[i] > r.maximum)
		{
			r.maximum = array[i];
		}

		if (array[i] < r.minimum)
		{
			r.minimum = array[i];
		}
	}
	return(r);
}



void CalculateLog(double * array, int size)
{
	// double val;
	for (int i = 0; i < size; ++i)
	{
		array[i] = log10(array[i]);
		// array[i] = val;
	}
}


COLOUR GetColour(double v,double vmin,double vmax)
{
   COLOUR c = {1.0, 1.0, 1.0}; // black
   double dv = vmax - vmin;

   if (v < (vmin + 0.25 * dv)) 
   {
      c.r = 0.0;
      c.g = 4 * (v - vmin) / dv;
   } 
   else if (v < (vmin + 0.5 * dv)) 
   {
      c.r = 0.0;
      c.b = 1 + 4 * (vmin + 0.25 * dv - v) / dv;
   } 
   else if (v < (vmin + 0.75 * dv)) 
   {
      c.r = 4 * (v - vmin - 0.5 * dv) / dv;
      c.b = 0.0;
   } 
   else 
   {
      c.g = 1 + 4 * (vmin + 0.75 * dv - v) / dv;
      c.b = 0.0;
   }

   return(c);
}

double Max(double * array, int size)
{
	double max = array[0];
	// int array_index = 0;

	for (int i = 0; i < size; ++i)
	{
		if (array[i] > max)
		{
			max = array[i];
			// array_index = i;
		}
	}

	// printf("Max: Array[%d] = %.17f\n", array_index, array[array_index]);
	return(max);
}

double Min(double* array, int size)
{
	// int array_index = 0;
	double min = array[0];
	for (int i = 0; i < size; ++i)
	{
		if (array[i] < min)
		{
			min = array[i];
			// array_index = i;
		}
	}
	// printf("Min: Array[%d] = %.17f\n", array_index, array[array_index]);
	return(min);
}