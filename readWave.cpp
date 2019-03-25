#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
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

struct package
{
	std::vector<double> time;
	std::vector<double> ch1;
	std::vector<double> ch2;
	std::vector<double> ch3;
};

double Gaus_2sigma( double *x, double*par )
{
	double xx = x[0];
	double height = par[0];
	double mean = par[1];
	double sigma_1 = par[2];
	double sigma_2 = par[3];
	 
	double Pol = par[4] + xx*(par[5] + xx*par[6] );
	  
	if( xx < mean ) 
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

  return ampl*exp(-arg/2) + par[4] + xx*(par[5] + xx * par[6]); 
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
	double f = par[0];
	return f;
}

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
	//std::cout << filebase << std::endl;
	std::vector<package> data; 

	double mina(0),minb(0),minc(0), maxa(0), maxb(0), maxc(0);
	int posa(0),posb(0),posc(0);
	
	for(int i=0; i < num ; i++)
	{
		std::ifstream file;	
		file.open(filebase + i + ".txt");
		
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
			a = a * scale;
			b = b * scale;
			c = c * scale;
			
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
				maxa = a;
			if(b < minb)
			{
				minb = b;
				posb = temp.ch2.size();
			}	
			if(b > maxb)
				maxb = b;
			if(c < minc)
			{
				minc = c;
				posc = temp.ch3.size();
			}	
			if(c > maxc)
				maxc = c;


		}
		data.push_back(temp);
		file.close();
	}


	/* move minimum to center of histogram */



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

	int maxunique = 0;
	if (uniquea > uniqueb)
		maxunique = uniquea;
	else maxunique = uniqueb;
	if (maxunique < uniquec)
		maxunique = uniquec;
	
	for(int i =0; i < num; i++)
		scaleTime(data[i].time);

	std::cout << "Time start: " << data[0].time[0] << std::endl;
	std::cout << "Time end: " << data[0].time.back() << std::endl;

	TApplication* rootapp = new TApplication("PMT timing", &argc, argv);
	TCanvas *c1 = new TCanvas("c1", "t");
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
	p1->Draw("SAME");

	TProfile *p2 = (TProfile*)p1->Clone();
	TCanvas* c2 = new TCanvas("c2", "Fit");
	p2->SetMarkerColor(kBlack);
	p2->Draw();


	//TF1 *f1 = new TF1("fit",myfunc, 88, 92, 2);
	//f1->SetLineColor(3);
	//p2->GetXaxis()->SetRange(p2->GetBin(88),p2->GetBin(92));
	//std::cout << p2->GetMean() << std::endl;
	//f1->SetParameters(p2->GetMean(),p2->GetRMS());
	//p2->GetXaxis()->SetRange();	
	//p2->Fit(f1);
	//f1->Draw("SAME");
		
	TF1 *f_Gaus_2sigma = new TF1("f_Gaus_2sigma", Gaus_2sigma, 0, 250, 7 );
	f_Gaus_2sigma->SetParameters(0, 100, 0.1, 0.5, 0, 0, 0);
	p2->Fit(f_Gaus_2sigma);
	f_Gaus_2sigma->Draw("SAME");	

	TF1 *f_Gaus_lau = new TF1("f_Gaus_lau", gaus_lau, 0, 250, 7 );
	f_Gaus_lau->SetParameters(10, 0.5, 0.5, 50, 0, 0, 0);
	f_Gaus_lau->SetLineColor(3);
	p2->Fit(f_Gaus_lau);
	f_Gaus_lau->Draw("SAME");	
	rootapp->Run();	

	return 0;
}


