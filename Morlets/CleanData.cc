 // MIT License

 // Copyright (c) [2017] [Vinay Yuvashankar]

 // Permission is hereby granted, free of charge, to any person obtaining a copy
 // of this software and associated documentation files (the "Software"), to deal
 // in the Software without restriction, including without limitation the rights
 // to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 // copies of the Software, and to permit persons to whom the Software is
 // furnished to do so, subject to the following conditions:

 // The above copyright notice and this permission notice shall be included in all
 // copies or substantial portions of the Software.

 // THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 // IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 // FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 // AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 // LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 // OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 // SOFTWARE.

/** 
    \file CleanData.cc
    \brief This file contains the code to remove any noise in the input data. 
*/


#include "processEEG.h"
#include <gsl/gsl_statistics.h>

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
