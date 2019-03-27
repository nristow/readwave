#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <deque>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <root/TGraph.h>
#include <root/TApplication.h>
#include <root/TH2D.h>
#include <root/TH1D.h>
#include <root/TProfile.h>
#include <root/TProfile2D.h>
#include <root/TCanvas.h>
#include <root/THStack.h>
#include <root/TF1.h>
#include <root/TFile.h>

struct package
{
	std::vector<double> time;
	std::vector<double> ch1;
	std::vector<double> ch2;
	std::vector<double> ch3;
};

/* Function from
 * http://www.nbi.dk/~petersen/Teaching/Stat2016/PythonRootIntro/ROOT_TipsAndTricks.pdf
 * */


double crystalball(double* x, double *par)
{
	//http://en.wikipedia.org/wiki/Crystal_Ball_function   
	double xcur = x[0];   
	double alpha = par[0];   
	double n = par[1];   
	double mu = par[2];   
	double sigma = par[3];   
	double N = par[4];   
	
	TF1* exp = new TF1("exp","exp(x)",1e-20,1e20);   
	
	double A; 
	double B;   
	if (alpha < 0)
	{     
		A = pow((n/(-1*alpha)),n)*exp->Eval((-1)*alpha*alpha/2);     
		B = n/(-1*alpha) + alpha;
	}   
	else 
	{     
		A = pow((n/alpha),n)*exp->Eval((-1)*alpha*alpha/2);     
		B = n/alpha - alpha;
	}

	double f;   
	if ((xcur-mu)/sigma > (-1)*alpha)     
		f = N*exp->Eval((-1)*(xcur-mu)*(xcur-mu)/(2*sigma*sigma));   
	else     
		f = N*A*pow((B- (xcur-mu)/sigma),(-1*n));   
	
	delete exp;   
	return f; 
}


double Gaus_2sigma( double *x, double*par )
{
	double xx = x[0];
	double height = par[0];
	double mean = par[1];
	double sigma_1 = par[2];
	double sigma_2 = par[3];
	double center = par[4]; 
	if( xx < center ) 
	{
		return height*exp( -(xx - mean)*(xx - mean)/(2*sigma_1*sigma_1) );
	}
	else
	{
		return height*exp( -(xx - mean)*(xx - mean)/(2*sigma_2*sigma_2) );
	}
}

Double_t gaus_lau (Double_t *x, Double_t *par) 
{ 
  static Double_t pi = 3.1415926535; 

  const Double_t xx =x[0]; 

  const Double_t width = (xx > par[0]) ? par[1] : par[2];
  const Double_t arg    = pow(((xx-par[0])/width),2);
  const Double_t ampl  = par[3];

  return ampl*exp(-arg/2); 
}

int findUnique(std::vector<double> v)
{
	std::sort(v.begin(), v.end());
	return std::unique(v.begin(), v.end()) - v.begin();
}

void scaleTime(std::vector<double> &v)
{
	for (int i = 0; i < v.size(); i++)
	{
		v[i] *= 1000000000;
	}
}

double findMean(std::vector<double> v)
{
	return std::accumulate(v.begin(), v.end(), 0.0)/ v.size();
}

Double_t myfunc(Double_t *x, Double_t *par)
{
	//double f=par[0] + (par[1]*x[0]);
	double xx = x[0];
	double x2 = xx*xx;
	double x3 = x2*xx;
	double x4 = x3*xx;
	double x5 = x4*xx;
	double x6 = x5*xx;
	return par[0] + par[1]*xx + par[2]*x2 + par[3]*x3 + par[4]*x4 + par[5]*x5 + par[6]*x6;
}

double linear(Double_t *x, Double_t *par)
{
	double slope = par[0];
	double b = par[1];
	double xx = x[0];

	return slope*xx+b;
}

//double newtonRaphson(double x, tf1* function, double EPSILON = 0.001) 
//{ 
//    double h = tf1->Eval(x) / tf1->Derivative(x); 
//    while (abs(h) >= EPSILON) 
//    { 
//        h = func(x)/derivFunc(x); 
//   
//        // x(i+1) = x(i) - f(x) / f'(x)   
//        x = x - h; 
//    } 
//  
//    cout << "The value of the root is : " << x;
//   return x; 
//} 

int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		std::cout << "Usage: readwave <path/to/files>" << std::endl;
		return 0;
	}

	const int num = 1000;
	std::string filebase=argv[1];
        filebase += "/waveform"; 
	std::vector<package> data; 

	double mina(0),minb(0),minc(0), maxa(0), maxb(0), maxc(0);
	double pmaxa(0),pmaxb(0),pmaxc(0);
	int posa(0),posb(0),posc(0);
	
	for(int i=0; i < num ; i++)
	{
		std::ifstream file;	
		std::stringstream ss;
		ss << filebase << i << ".txt";
		file.open(ss.str());
		
		std::vector<double> time;
		std::vector<double> ch1;
		std::vector<double> ch2;
		std::vector<double> ch3;


		for(int i = 0; i < 10; i++)
		{	
			std::string line;
			std::getline(file, line);
		}

		package temp;
		const double scale = 1000; // scale to mV level
		for (std::string line; std::getline(file, line); )
		{
			std::stringstream ss(line);
			double d, a, b, c;
			ss >> d >> std::ws >> a >> std::ws >> b >> std::ws >> c;
			/* all 3 channels are inverted */
			a = -a * scale;
			b = -b * scale;
			c = -c * scale;
			
			temp.time.push_back(d);
			temp.ch1.push_back(a);
			temp.ch2.push_back(b);
			temp.ch3.push_back(c);

			if(a < mina)
			{
				mina = a;
				posa = temp.ch1.size();
			}	
			if(a > maxa)
			{
				maxa = a;
				pmaxa = temp.ch1.size();
			}
			if(b < minb)
			{
				minb = b;
				posb = temp.ch2.size();
			}	
			if(b > maxb)
			{
				maxb = b;
				pmaxb = temp.ch2.size();
			}
			if(c < minc)
			{
				minc = c;
				posc = temp.ch3.size();
			}	
			if(c > maxc)
			{
				maxc = c;
				pmaxc = temp.ch3.size();
			}


		}
		data.push_back(temp);
		file.close();
	}


	int uniquea(0), uniqueb(0), uniquec(0);
	for(int i = 0; i < num; i++)
	{
			int temp = findUnique(data[i].ch1);
			if (temp > uniquea)
				uniquea = temp;
			
			temp = findUnique(data[i].ch2);
			if (temp > uniqueb)
				uniqueb = temp;
			
			temp = findUnique(data[i].ch3);
			if (temp > uniquec)
				uniquec = temp;
	}		


	std::cout << "Unique counts: " << uniquea << "," << uniqueb << "," << uniquec << std::endl;		
	std::cout << "mina: " << mina << " minb: " << minb << " minc: " << minc << std::endl;
	std::cout << "maxa: " << maxa << " maxb: " << maxb << " maxc: " << maxc << std::endl;
	std::cout << "minimum position a: " << posa << " minimum position b: " << posb << " minimum position c: " << posc << std::endl;
	std::cout << "maximum position a: " << pmaxa << " maximum position b: " << pmaxb << " maximum position c: " << pmaxc << std::endl;

	int maxunique = 0;
	if (uniquea > uniqueb)
		maxunique = uniquea;
	else maxunique = uniqueb;
	if (maxunique < uniquec)
		maxunique = uniquec;
	
	for(int i =0; i < num; i++)
		scaleTime(data[i].time);

	double timestart = data[0].time[0];
	double timeend = data[0].time.back();
	std::cout << "Time start: " << data[0].time[0] << std::endl;
	std::cout << "Time end: " << data[0].time.back() << std::endl;

	std::vector<float> timech1;
	std::vector<float> ch1scaled;

	TApplication* rootapp = new TApplication("PMT timing", &argc, argv);
	TFile *file  = new TFile("PMT_fit.root", "RECREATE");
	
	TCanvas *c2 = new TCanvas("c2","test");	
	TH2D *h2 = new TH2D("h2", "PMT-1A response from LED pulse", data[0].time.size(), data[0].time.front(), data[0].time.back(), maxunique , mina, maxa);

	for (int j = 0; j < num; j++)
	{
		for (int i = 0; i < data[j].time.size(); i++)
		{
			h2->Fill(data[j].time[i], data[j].ch1[i]);
		}
	}	

	h2->SetMarkerColor(kGray);
	h2->SetOption("scat");
	h2->GetYaxis()->SetTitle("Voltage [mV]");
	h2->GetXaxis()->SetTitle("Time [nS]");
	h2->Draw();	
	TCanvas *c1 = new TCanvas("c1","test");	
	TProfile *p1 = new TProfile("p1", "Profile of CH1", data[0].time.size(), data[0].time.front(), data[0].time.back());

	for (int j = 0; j < num; j++)
	{
		for (int i = 0; i < data[j].time.size(); i++)
		{
			p1->Fill(data[j].time[i], data[j].ch1[i]);
		}
	}	


      	
	p1->SetOption("hist");
	p1->SetMarkerColor(kGreen);
	p1->Draw("P");
	
	double mintimepos = ((TAxis*)p1->GetXaxis())->GetBinCenter(posa);//-1.076E-007+((posa+1)*(4E-10));
	double maxtimepos = ((TAxis*)p1->GetXaxis())->GetBinCenter(pmaxa);//-1.076E-007+((posa+1)*(4E-10));
	double maximum(0);
	for(int i = 1; i < num; i++)
	{
		double temp = p1->GetBinContent(i);
		if (temp > maximum)
			maximum = temp;
	}
	
	
	std::cout << "Time at which the function is a maximum: " << maxtimepos << std::endl << "Value of the maximum: " << maximum << std::endl;

	TF1 *f_Gaus_2sigma = new TF1("f_Gaus_2sigma", Gaus_2sigma, 80, 130, 5 );
	f_Gaus_2sigma->SetParameters(60, 110, 10, 1, 100.3);
	//f_Gaus_2sigma->FixParameter(5,maxtimepos);
	//f_Gaus_2sigma->FixParameter(0,maximum);
	//p1->Fit("f_Gaus_2sigma");
	//f_Gaus_2sigma->Draw("SAME");

	TF1 *f = new TF1("sixthorderpoly", myfunc, 80, 130, 7 );
	f->SetParameters(-285,3.3,1,1,1,1,1);
	//p1->Fit("sixthorderpoly");
	//f->Draw("SAME");

	TF1 *f2 = new TF1("linearfit", linear, 0, 200, 2 );
	f->SetParameters(-285,3.3);
	//p1->Fit("linearfit");
	//f2->Draw("SAME");
	
	TF1 *f3 = new TF1("CrystalBall", crystalball, 80,130,5);
	rootapp->Run();	
	file->Write();	
	file->Close();
	delete file, h2, c1, p1,f ,f2;
	delete f_Gaus_2sigma;


	return 0;
}


