
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
        double c3 = par[5];
        
        double z = 1+exp((t0-t)/alpha);
        return (A+c3 + c2*t + c1*t*t)/z;
}

TF1 *ttttt = new TF1("fermi_dirac", fermi_dirac, 0, 9, 3);
TF1 *f1 = new TF1("fermi_dirac_parabola", fermi_dirac_parabola, 0, 9.5, 6);
//f1->SetParNames("Alpha","t0","A", "c1", "c2", "c3");
//f1->SetParameters(1.0,8,108, 0, 0, 0);

