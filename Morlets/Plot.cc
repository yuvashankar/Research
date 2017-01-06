/**
	\file "Plot.cc"
	\brief All the functions needed to plot the png. 
*/
#include "wavelet.h"
#include <pngwriter.h>
#include <float.h>
#include <cmath>

//Plotting Constants
/**
	\var PLOT_OY
	\brief The amount of vertical black space in the plot
*/
#define PLOT_OY 200

/**
	\var PLOT_OX
	\brief The amount of horizontal black space in the plot
*/
#define PLOT_OX 200


//Data Structures
/**
	\struct COLOUR
	\brief the RGB representation of every pixel
*/
typedef struct 
{
	double r,g,b;
} COLOUR;

/**
	\struct RANGE
	\brief The lower and upper bound of the data matrix
*/
typedef struct
{
	double minimum, maximum;
} RANGE;


RANGE GetRange(double* array, int size);
void CalculateLog(double * array, int size);
COLOUR GetColour(double v,RANGE data_range);
double Max(double * array, int size);
double Min(double* array, int size);

int Plot(double * data, double * frequency, int num_x, int num_y, int plot_type,
	char graph_title[],
	char filename[])
{
	switch(plot_type)
	{
		case 0: //Plot_PNG
			Plot_PNG(data, frequency, num_x, num_y, graph_title, filename);
			break;
		case 1:
			WriteFile(data, frequency, num_x, num_y, filename);
	}
	return(0);
}

void Plot_PNG(double * data, double * periods, int num_x, int num_y, char graph_title[], 
	char filename[])
{
	int i, j, k;
	
	const int lines_size = 10;
	const int stride = 4; //Stride needs to be even. 
	const int image_width  = (num_x/stride) + 2 * PLOT_OX;
	const int image_height = lines_size * num_y+ 2 * PLOT_OY;
	const int height = lines_size * num_y;

	const int label_font_size = 30;
	const int tic_font_size   = 15;
	const int title_font_size = 40;
	const int freq_text_width = 15;
	const int time_text_width = 15;


	char  font_location[] = "../lib/VeraMono.ttf";
	char x_label[] = "Time (s)";
	char y_label[] = "Frequency (Hz)";
	// char graph_title[] = "Time Frequency Graph of an Impulse";
	char temp_string[4];

	CalculateLog( data, num_x * num_y );

	RANGE r = GetRange(data, num_x*num_y);
	// printf("New Max = %f, New Min = %f\n", r.maximum, r.minimum);

	pngwriter png( image_width, image_height , 0 , filename);

	// X_label
	int x_text_width = png.get_text_width(font_location, label_font_size, x_label);
	png.plot_text( font_location, label_font_size,
				(0.5*image_width - 0.5*x_text_width), 0.5*PLOT_OY, 0.0,
				x_label,
				1.0, 1.0, 1.0);
	
	//Y_Label
	int y_text_width = png.get_text_width(font_location, label_font_size, y_label);
	png.plot_text( font_location, label_font_size,
				0.25*PLOT_OY, (0.5*image_height - 0.5*y_text_width), M_PI*0.5,
				y_label,
				1.0, 1.0, 1.0);

	//Title
	int title_text_width = png.get_text_width(font_location, label_font_size, graph_title);
	png.plot_text( font_location, title_font_size,
				(0.5*image_width - title_text_width), (image_height - 0.5*PLOT_OY), 0.0,
				graph_title,
				1.0, 1.0, 1.0);

	//Plot the Graph itself.
	for ( i = 1; i <= num_x/stride; ++i)
	{
		for ( j = 1; j <= num_y; ++j)
		{
			for ( k = 0; k < lines_size; ++k)
			{
				int counterVar = stride * (i-1) + (j-1) *num_x;
				double value = data[counterVar];


				
				//Ensure that the indexes do not exceed the data limit.
				assert( counterVar <= num_x*num_y );

				COLOUR c = GetColour(value, r);
				
				png.plot(PLOT_OX + i , PLOT_OY + (height - j*lines_size + k), 
					c.r, c.g, c.b);
					// 0.0, 0.0, 0.0);
			}

			//Add Frequency Markers
			if (j % 15 == 0)
			{
				png.filledsquare( PLOT_OX - 20, PLOT_OY + (height - j*lines_size) - 4,
						PLOT_OX, PLOT_OY + (height - j*lines_size) + 4,
						1.0, 1.0, 1.0);

				sprintf(temp_string, "%.1f", periods[j]);

				// int freq_text_width = png.get_text_width(font_location, tic_font_size, temp_string);
				// int freq_text_width = 15;
				png.plot_text(font_location, tic_font_size, 
					PLOT_OX - 30 - freq_text_width, PLOT_OY + (height - j*lines_size) - 0.5*tic_font_size, 0.0, 
					temp_string, 
					1.0, 1.0, 1.0);
			}
		}
		// ADD A TIME MARKER
		if( (stride * i) % FS == 0 )
		{
			png.filledsquare( PLOT_OX + i - 4, PLOT_OY,
							  PLOT_OX + i + 4, PLOT_OY - 20,
							1.0, 1.0, 1.0);

			sprintf(temp_string, "%.1d", (stride * i)/FS);
			
			// int time_text_width = png.get_text_width(font_location, tic_font_size, temp_string);
			// int time_text_width = 15;
			png.plot_text( font_location, tic_font_size,
				PLOT_OX + i - 4 - 0.5* time_text_width, PLOT_OY - 50, 0.0,
				temp_string,
				1.0, 1.0, 1.0);
		}

	}
	png.close();
}

/**
	\fn int WriteFile(const double *data, const double *period, const int x, const int y, const char* filename)

	\brief A function that writes the Wavelet Results to the disk. 

	\param data A x x y array with the  data that is going to be written 
	\param period A 1 x y array with the frequencies that were analyzed
	\param x The number of samples in the signal
	\param y The number of frequencies analyzed
	\param filename The name of the file that will be written

	\return 0 if successful
	\return -1 if unsuccessful

	This function will write the resultant data computed by Wavelet() and ERSP() into the disk so that it can be graphed by Gnuplot. 
	One can plot the output of this function using the matrix.gplot file. 
*/
int WriteFile(const double *data, const double *frequency, const int x, const int y, 
	char filename[])
{

    FILE* out_file=fopen(filename,"w");
    if (out_file == NULL) return -1;

    //Xticks
    fprintf(out_file, "%d\t", x);
    for (int i = 0; i < y; ++i)
    {
    	fprintf(out_file, "%f\t", (double) i/FS);
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
	
	fprintf(gnuplot_file, "set term x11\n");
	fprintf(gnuplot_file, "set pm3d map\n");
	fprintf(gnuplot_file, "set logscale z 10\n");
	fprintf(gnuplot_file, "set logscale y 2\n");
	fprintf(gnuplot_file, "set ticslevel 0\n");
	fprintf(gnuplot_file, "set xlabel \"time (s)\"\n");
	fprintf(gnuplot_file, "set ylabel \"Frequency (Hz)\"\n");

	//Input Graph Title
	fprintf(gnuplot_file, "%s%s%s\n", "set title \"", graph_title, "\"");
	//Plot filename
	// plot "DATA.log" matrix nonuniform with pm3d t ''
	fprintf(gnuplot_file, "%s%s%s\n", "splot \"", filename, "\" matrix nonuniform with pm3d t ''");

	return(0);
}

/**
	\fn RANGE GetRange(double* array, int size)

	\param array The data array that will be plotted
	\param size The total size of the contiguous memory

	\return RANGE The maximum and minimum of the data array.

	This function will iterate through the entire function and returns the maximum and minimum of the entire array.
*/
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


/**
	\fn void CalculateLog(double * array, int size)

	\param array The array that needs to be computed
	\param size The size of the contiguous block of memory

	This function will iterate through every element in the array and compute the logarithm. 
	This will override the array.
*/
void CalculateLog(double * array, int size)
{
	for (int i = 0; i < size; ++i)
	{
		array[i] = log10(array[i]);

	}
}

/**
	\fn COLOUR GetColour(double v,RANGE data_range)

	\param v the value of the pixel to be plotted
	\param data_range The range of the given data

	\return pixel_colour The colour of the pixel that will be plotted

	This function takes a double and maps to a colour map. High values are closer to the red colour spectrum, and low values are mapped to the blue colour spectrum. 
*/

COLOUR GetColour(double v, RANGE data_range)
{
   COLOUR c = {1.0, 1.0, 1.0}; // white

   if (v > data_range.maximum)
   	v = data_range.maximum;
   
   if (v < data_range.minimum)
   	v = data_range.minimum;

   double dv = data_range.maximum - data_range.minimum;

   if (v < (data_range.minimum + 0.25 * dv)) 
   {
      c.r = 0.0;
      c.g = 4 * (v - data_range.minimum) / dv;
   } 
   else if (v < (data_range.minimum + 0.5 * dv)) 
   {
      c.r = 0.0;
      c.b = 1 + 4 * (data_range.minimum + 0.25 * dv - v) / dv;
   } 
   else if (v < (data_range.minimum + 0.75 * dv)) 
   {
      c.r = 4 * (v - data_range.minimum - 0.5 * dv) / dv;
      c.b = 0.0;
   } 
   else 
   {
      c.g = 1 + 4 * (data_range.minimum + 0.75 * dv - v) / dv;
      c.b = 0.0;
   }

   return(c);
}


/**
	\fn double Max(double* array, int size)

	\brief A function that finds the Maximum of a given array
	\param array The array to be analyzed
	\param size The size of the array
*/
double Max(double * array, int size)
{
	double max = array[0];
	int array_index = 0;

	for (int i = 0; i < size; ++i)
	{
		if (array[i] > max)
		{
			max = array[i];
			if (max != max)
			{
				printf("Naan Alert!\n");
			}
			array_index = i;
		}
	}

	printf("Max: Array[%d] = %.17f\n", array_index, array[array_index]);
	return(max);
}

/**
	\fn double Min(double* array, int size)

	\brief A function that finds the minimum of a given array
	\param array The array to be analyzed
	\param size The size of the array
*/
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