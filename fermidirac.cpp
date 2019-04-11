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

double fermi_dirac_parabola(double *x, double *par)
{
	double t = x[0];
	double alpha = par[0];
	double t0 = par[1];
	double A = par[2];
        double c1 = par[3];
        double c2 = par[4];
        //double c3 = par[5];
        
        double z = 1+exp((t0-t)/alpha);
        return (A+ c2*t + c1*t*t)/z;
}
