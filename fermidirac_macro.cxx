
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
        double B = par[3];
        double C = par[4];
        double D = par[5];

	double z = 1+exp((t0-t)/alpha);
	z = A/z;
	return z;
}

TF1 *f_fit = new TF1("fermi_dirac", fermi_dirac, 0, 9, 3);

