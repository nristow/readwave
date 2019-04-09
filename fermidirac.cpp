#include "fermidirac.h"

double fermi_dirac(double *x, double *par)
{
	double t = x[0];
	double alpha = par[0];
	double t0 = par[1];
	double A = par[2];

	double z = 1+exp((t0-t)/alpha);
	z = A/z;
	return z;
}
