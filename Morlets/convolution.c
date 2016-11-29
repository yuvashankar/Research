#include "wavelet.h"

double CompleteRealMorlet (double x, double scale)
{
	double w_0 = 6.0;
	double c_sigma = pow( (1.0 + exp(-w_0*w_0) - 2.0 * exp(-0.75 * w_0 * w_0)), -0.5 );
	double k_sigma = exp( -0.5 * w_0 * w_0 );

	x = x * scale; 

	double out = exp( -0.5 * x * x) * ( sin(w_0 * x) - k_sigma );
	out = c_sigma * QUAD_ROOT_PI * out;
	return(out);
}

double CompleteComplexMorlet(double x, double scale)
{
	double w_0 = 6.0;
	double c_sigma = pow( (1.0 + exp(-w_0*w_0) - 2.0 * exp(-0.75 * w_0 * w_0)), -0.5 );
	double k_sigma = exp( -0.5 * w_0 * w_0 );
	x = x * scale; 

	double out = exp( -0.5 * x * x) * ( cos(w_0 * x) - k_sigma );

	out = c_sigma * QUAD_ROOT_PI * out;
	return(out);
}