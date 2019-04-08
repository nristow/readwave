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
#include "countfiles.h"

struct package
{
	std::vector<double> time;
	std::vector<double> ch1;
	std::vector<double> ch2;
	std::vector<double> ch3;
	std::vector<double> max; // 0=ch1 1=ch2 ...
	std::vector<double> min;
	std::vector<double> numunique;
	std::vector<double> ymulti;
	std::vector<double> filteredmax = {0,0};
	double maxdifference;
};


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

int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		std::cout << "Usage: readwave <path/to/files>" << std::endl;
		return 0;
	}

	std::string filebase=argv[1];
	const int num = numFiles(filebase);
	if (num == -1)
	{
		std::cout << "Directory doesn't exist" << std::endl;
		return -1;
	}

	std::cout << "Parsing " << num << " files" << std::endl;	
	filebase += "/waveform"; 
	std::vector<package> data; 

	int posa(0),posb(0),posc(0);
	bool oldheader = false;
	
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

		package temp;
		if(oldheader == false)
		{
			for(int i = 0; i < 12; i++)
			{	
				std::string line;
				std::getline(file, line);
				if(i == 3)
				{
					std::stringstream ss(line);
					std::string trash;
					std::string ym;
					ss >> trash >> trash >> trash >> ym; 
					temp.ymulti.push_back(std::stod(ym));
				}
				else if (i == 4)
				{
					std::stringstream ss(line);
					std::string trash;
					std::string ym;
					ss >> trash >> trash >> trash >> ym; 
					temp.ymulti.push_back(std::stod(ym));
				}

			}
		}
		else
		{
			for(int i = 0; i < 10; i++)
			{	
				std::string line;
				std::getline(file,line);
			}
	}
		const double scale = 1000; // scale to mV level
		double mina(100),minb(100),minc(100), maxa(-100), maxb(-100), maxc(-100);
		double pmaxa(0),pmaxb(0),pmaxc(0);
		double uniquea(0), uniqueb(0), uniquec(0);
		for (std::string line; std::getline(file, line); )
		{
			std::stringstream ss(line);
			double d(0), a(0), b(0), c(0);
			ss >> d >> std::ws >> a >> std::ws >> b;
			/* all 3 channels are inverted */
			a = -a * scale;
			b = -b * scale;
			//c = -c * scale;
			
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
			//if(c < minc)
			//{
			//	minc = c;
			//	posc = temp.ch3.size();
			//}
			//if(c > maxc)
			//{
			//	maxc = c;
			//	pmaxc = temp.ch3.size();
			//}

		}

		int n=findUnique(temp.ch1);
		if (n > uniquea)
			uniquea = n;
	
		n=findUnique(temp.ch2);	
		if (n > uniqueb)
			uniqueb = n;
	//	n=findUnique(temp.ch3);
	//	if (n > uniquec)
	//		uniquec = n;

		temp.max.push_back(maxa);
		temp.max.push_back(maxb);
		//temp.max.push_back(maxc);
		temp.min.push_back(mina);
		temp.min.push_back(minb);
		//temp.min.push_back(minc);
		temp.numunique.push_back(uniquea);
		temp.numunique.push_back(uniqueb);
		//temp.numunique.push_back(uniquec);
		data.push_back(temp);
		file.close();
	}

//	TApplication* rootapp = new TApplication("PMT timing", &argc, argv);
	
	std::string filename = argv[1];
	filename.erase(0,3);
	filename = "data/" + filename + ".root";

	TFile *file  = new TFile(filename.c_str(), "RECREATE");
	TDirectory *ch1 = file->mkdir("ch1");
	TDirectory *ch2 = file->mkdir("ch2");
	TDirectory *ch3 = file->mkdir("ch3");
	TDirectory *directories[] = {ch1, ch2, ch3};
	TDirectory *stats = file->mkdir("stats");
	TDirectory *ch1graphs = file->mkdir("ch1 graphs");
	TDirectory *timing = file->mkdir("Timing");

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
		data[i].maxdifference = abs(data[i].max[0] - data[i].max[1]);

		/* sliding mean */
		std::deque<double> window(11,0.0);
		for(int k = 0; k < 2; k++)
		{
			double tempmax{};
			for(int j = 0; j < data[i].time.size(); j++)
			{
				if(k == 0)
					window.push_back(data[i].ch1[j]);
				else if(k == 1)
					window.push_back(data[i].ch2[j]);
				window.pop_front();
				double temp = std::accumulate(window.begin(), window.end(), 0.0)/window.size();
				if(temp > tempmax)
				{
					tempmax = temp;
				}
			}

			data[i].filteredmax[k] = tempmax; 
			if(tempmax < 160)
				std::cout << "Maximum of " << tempmax << " at " << i << std::endl;
		}

	}

	std::cout << "Number of unique heights in ch1: " << uniquech1 << " ch2: " << uniquech2 << std::endl;
	std::cout << "Max height in ch1: " << maxch1 << " ch2: " << maxch2 << std::endl;
	std::cout << "smallestmax height in ch1: " << smallestmaxch1 << " ch2: " << smallestmaxch2 << std::endl;
	
	/* compute size of bins */
	
	if(oldheader == false)
	{
	double numbinsch1 = 1000*data[0].ymulti[0];
	double numbinsch2 = 1000*data[0].ymulti[1];
	std::cout << "Size of bins ch1: " << numbinsch1 << std::endl;
	std::cout << "Size of bins ch2: " << numbinsch2 << std::endl;
	}

	/* construct histogram of pulse heights */
	
	stats->cd();
	std::string histname = "PH2";
	std::string info = "Pulse Height for PMT-1B";
	TH1D *h1 = new TH1D(histname.c_str(), info.c_str(), 1000, -1, 400);
	for(int i = 0; i < num ; i++)
	{
		h1->Fill(data[i].filteredmax[1]);
	}
	h1->GetXaxis()->SetTitle("Maximum Waveform Height (mV)");
	h1->GetYaxis()->SetTitle("Number of waveforms");
	h1->Write();
	delete h1;
	
	histname = "PH1";
	info = "Pulse Height for PMT-1A";
	TH1D *h2 = new TH1D(histname.c_str(), info.c_str(), 1000, -1, 400);
	for(int i = 0; i < num ; i++)
	{
		h2->Fill(data[i].filteredmax[0]);
	}
	h2->GetXaxis()->SetTitle("Maximum Waveform Height (mV)");
	h2->GetYaxis()->SetTitle("Number of waveforms");
	h2->Write();
	delete h2;
	
	/* Plot difference in height between all pulses */

	std::string histname2 = "Difference";
	std::string info2 = "Absolute difference in pulse height";
	TH1D *h3 = new TH1D(histname2.c_str(), info2.c_str(), 254, 0, 204.4);
	for(int i = 0; i < num ; i++)
	{
		h3->Fill(data[i].maxdifference);
	}
	h3->GetXaxis()->SetTitle("Absolute difference (mV)");
	h3->GetYaxis()->SetTitle("Number of waveforms");
	h3->Write();
	delete h3;

	/* Graph 2d amplitudes of PMT's */
	std::vector<double> max1;
	std::vector<double> max2;

	for(int i = 0; i < data.size(); i++)
	{
		max1.push_back(data[i].filteredmax[0]);
		max2.push_back(data[i].filteredmax[1]);
	}
	
	TCanvas *c1 = new TCanvas("XY Plot", "PMT maximums XY plot");
	TGraph * g1 = new TGraph(data.size(),&max1[0],&max2[0]);
	g1->GetXaxis()->SetTitle("PMT-1A maximum (mV)");
	g1->GetYaxis()->SetTitle("PMT-1B maximum (mV)");
	g1->SetDrawOption("AP");
	g1->Draw("AP");
	c1->Write();
	
	ch1graphs->cd();

	for(int i = 0; i < num; i++)
	{
		std::vector<double> x;
		std::vector<double> y;
		for(int j = 0; j < data[i].time.size(); j++)
		{	
			y.push_back(data[i].ch1[j]);
			x.push_back(data[i].time[j]);
		}
		//TCanvas *c2 = new TCanvas("XY waveform", "PMT-1A waveform XY plot");
		TGraph * g2 = new TGraph(data[0].time.size(),&x[0],&y[0]);
		g2->GetXaxis()->SetTitle("Time (ns)");
		g2->GetYaxis()->SetTitle("Amplitude (mV)");
		g2->SetName(Form("g%d",i) );
		g2->SetDrawOption("AP");
		g2->Draw("AP");
		g2->Write();
		delete g2;
	}
	
	//g1->Write();

	/* Construct histogram of all pulses */	
//	for (int i = 0; i < num; i++)
//	{
//		for( int j = 0; j < 2; j++) //per channel
//		{
//			directories[j]->cd();
//			std::string histname = "h" + std::to_string(i+1);
//			//std::cout << histname << std::endl;
//			std::string info = "a Pulse";
//			TH2D *h1 = new TH2D(histname.c_str(), info.c_str(), data[i].time.size(), data[i].time.front(), data[i].time.back(), 255, data[i].min[j],data[i].max[j]); 
//			if (j == 0)
//			{
//				for(int k = 0; k < data[i].time.size(); k++)
//					h1->Fill(data[i].time[k], data[i].ch1[k]);
//			}
//			else	
//			{
//				for(int k = 0; k < data[i].time.size(); k++)
//					h1->Fill(data[i].time[k], data[i].ch2[k]);
//			}
//
//			h1->Write();
//			delete h1;
//		}
//	}


	file->Write();	
	file->Close();
	delete file;


	return 0;
}


