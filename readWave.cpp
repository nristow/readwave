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
#include <root/TSpline.h>

struct package
{
	std::vector<double> time;
	std::vector<double> ch1;
	std::vector<double> ch2;
	std::vector<double> ch3;
	std::vector<double> max;
	std::vector<double> min;
	std::vector<double> numunique;
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

double Gaus_2mean( double *x, double*par )
{
	double xx = x[0];
	double height = par[0];
	double mean1 = par[1];
	double mean2 = par[2];
	double sigma_1 = par[3];
	double sigma_2 = par[4];
	double center = par[5]; 
	if( xx < center ) 
	{
		return height*exp( -(xx - mean1)*(xx - mean1)/(2*sigma_1*sigma_1) );
	}
	else
	{
		return height*exp( -(xx - mean2)*(xx - mean2)/(2*sigma_2*sigma_2) );
	}
}

double cubicSpline(double *x, double *p)
{
	double xx = x[0];
	double x2 = xx*xx;
	double x3 = x2*xx;

	double split = p[0];

	double f0 = p[1]+p[2]*xx + p[3]*x2 + p[4]*x3;
	double f1 = p[5]+p[6]*xx + p[8]*x2 + p[8]*x3;
//	double f2 = p[9]+p[10]*xx + p[11]*x2 + p[12]*x3;
	
	if (xx < split)
		return f0; 
//	else if (xx == split)
//		return f0-f1;
	else
		return f1;

}

Double_t spline_4nodes(Double_t *x, Double_t *par)
{
   /*Fit parameters:
   par[0-3]=X of nodes (to be fixed in the fit!)
   par[4-7]=Y of nodes
   par[8-9]=first derivative at begin and end (to be fixed in the fit!)
   */
   Double_t xx = x[0];

   Double_t xn[6] = { par[0], par[1], par[2], par[3], par[4], par[5] };
   Double_t yn[6] = { par[6], par[7], par[8], par[9], par[10], par[11] };

   Double_t b1 = par[12];
   Double_t e1 = par[13];

   TSpline3 sp3("sp3", xn, yn, 6, "b1e1", b1, e1);

   return sp3.Eval(xx);
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

	const int num = 1;
	std::string filebase=argv[1];
        filebase += "/waveform"; 
	std::vector<package> data; 

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
		double mina(100),minb(100),minc(100), maxa(-100), maxb(-100), maxc(-100);
		double pmaxa(0),pmaxb(0),pmaxc(0);
		double uniquea(0), uniqueb(0), uniquec(0);
		for (std::string line; std::getline(file, line); )
		{
			std::stringstream ss(line);
			double d(0), a(0), b(0), c(0);
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

		int n=findUnique(temp.ch1);
		if (n > uniquea)
			uniquea = n;
	
		n=findUnique(temp.ch2);	
		if (n > uniqueb)
			uniqueb = n;
		n=findUnique(temp.ch3);
		if (n > uniquec)
			uniquec = n;

		temp.max.push_back(maxa);
		temp.max.push_back(maxb);
		temp.max.push_back(maxc);
		temp.min.push_back(mina);
		temp.min.push_back(minb);
		temp.min.push_back(minc);
		temp.numunique.push_back(uniquea);
		temp.numunique.push_back(uniqueb);
		temp.numunique.push_back(uniquec);
		data.push_back(temp);
		file.close();
	}

//	TApplication* rootapp = new TApplication("PMT timing", &argc, argv);
	TFile *file  = new TFile("PMT_muon_PEN.root", "RECREATE");
	TDirectory *ch1 = file->mkdir("ch1");
	TDirectory *ch2 = file->mkdir("ch2");
	TDirectory *ch3 = file->mkdir("ch3");
	TDirectory *directories[] = {ch1, ch2, ch3};
	TDirectory *stats = file->mkdir("stats");

	double uniquech1(0), uniquech2(0), maxch1(0), maxch2(0), smallestmaxch1(180), smallestmaxch2(180), minch1(180), minch2(180);
	for(int i = 0; i < num; i++)
	{
		if (data[i].numunique[0] > uniquech1)
			uniquech1 = data[i].numunique[0];
		if (data[i].numunique[1] > uniquech2)
			uniquech2 = data[i].numunique[1];
		if(data[i].max[0] > maxch1)
			maxch1 = data[i].max[0];
		if(data[i].max[1] > maxch2)
			maxch2 = data[i].max[1];
		if(data[i].max[0] < smallestmaxch1)
			smallestmaxch1 = data[i].max[0];
		if(data[i].max[1] <  smallestmaxch2)
			smallestmaxch2 = data[i].max[1];
		if(data[i].min[0] < minch1)
			minch1 = data[i].min[0];
		if(data[i].min[1] <  minch2)
			minch2 = data[i].min[1];
	}

	std::cout << "Number of unique heights in ch1: " << uniquech1 << " ch2: " << uniquech2 << std::endl;
	std::cout << "Max height in ch1: " << maxch1 << " ch2: " << maxch2 << std::endl;
	std::cout << "smallestmax height in ch1: " << smallestmaxch1 << " ch2: " << smallestmaxch2 << std::endl;
	
	/* construct histogram of pulse heights */
	
	stats->cd();
	std::string histname = "PH2";
	std::string info = "Pulse Height for PMT-1B";
	TH1D *h1 = new TH1D(histname.c_str(), info.c_str(), 254, smallestmaxch2, maxch2);
	for(int i = 0; i < num ; i++)
	{
		h1->Fill(data[i].max[1]);
	}
	h1->GetXaxis()->SetTitle("Maximum Waveform Height (mV)");
	h1->GetYaxis()->SetTitle("Number of waveforms");
	h1->Write();
	delete h1;
	
	histname = "PH1";
	info = "Pulse Height for PMT-1A";
	TH1D *h2 = new TH1D(histname.c_str(), info.c_str(), 254, smallestmaxch1, maxch1);
	for(int i = 0; i < num ; i++)
	{
		h2->Fill(data[i].max[0]);
	}
	h2->GetXaxis()->SetTitle("Maximum Waveform Height (mV)");
	h2->GetYaxis()->SetTitle("Number of waveforms");
	h2->Write();
	delete h2;

	/* Construct histogram of all pulses */	
	for (int i = 0; i < num; i++)
	{
		for( int j = 0; j < 2; j++) //per channel
		{
			directories[j]->cd();
			std::string histname = "h" + std::to_string(i+1);
			//std::cout << histname << std::endl;
			std::string info = "a Pulse";
			TH1D *h1 = new TH1D(histname.c_str(), info.c_str(), /*data[i].numunique[j]*/255, data[i].min[j],data[i].max[j]); 
			if (j == 0)
			{
				for(int k = 0; k < data[i].time.size(); k++)
					h1->Fill(data[i].ch1[k]);
			}
			else	
			{
				for(int k = 0; k < data[i].time.size(); k++)
					h1->Fill(data[i].ch2[k]);
			}

			h1->Write();
			delete h1;
		}
	}
	
	for (int i = 0; i < 1; i++)
	{
		
		TDirectory *bins = file->mkdir("bins");
		bins->cd();
		//std::cout << "minimum: " << data[1].min[0];
		for(int j = 1; j < 1024; j++)
		{
			std::string histname = "h" + std::to_string(j);
			//std::cout << histname << std::endl;
			std::string info = "a Pulse";
			TH1D *h1 = new TH1D(histname.c_str(), info.c_str(), j, data[0].min[0], data[0].max[0]); 
			for(int k = 0; k < data[0].time.size(); k++)
				h1->Fill(data[0].ch1[k]);
			h1->Write();
			delete h1;
		}
	}
//		
//		
//		
//		
//		
//		TCanvas *c2 = new TCanvas("c2","test");	
//		TH2D *h2 = new TH2D("h2", "A Pulse", data[0].time.size(), data[0].time.front(), data[0].time.back(), maxunique , mina, maxa);
//
//	for (int j = 0; j < num; j++)
//	{
//		for (int i = 0; i < data[j].time.size(); i++)
//		{
//			h2->Fill(data[j].time[i], data[j].ch1[i]);
//		}
//	}	
//
//	h2->SetMarkerColor(kGray);
//	h2->SetOption("scat");
//	h2->GetYaxis()->SetTitle("Voltage [mV]");
//	h2->GetXaxis()->SetTitle("Time [nS]");
//	h2->Draw();	
//	
//	//rootapp->Run();	
	file->Write();	
	file->Close();
	delete file;


	return 0;
}


