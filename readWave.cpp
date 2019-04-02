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
	std::vector<double> ymulti;
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

	const int num = 11000;
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

		package temp;
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
	
	/* compute number of bins */
	
	double numbinsch1 = 1000*data[0].ymulti[0];
	double numbinsch2 = 1000*data[0].ymulti[1];
	std::cout << "Number of bins ch1: " << numbinsch1 << std::endl;
	std::cout << "Number of bins ch2: " << numbinsch2 << std::endl;

	/* construct histogram of pulse heights */
	
	stats->cd();
	std::string histname = "PH2";
	std::string info = "Pulse Height for PMT-1B";
	TH1D *h1 = new TH1D(histname.c_str(), info.c_str(), 253, 0, 102);
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
	TH1D *h2 = new TH1D(histname.c_str(), info.c_str(), 253, 0, 204);
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
			TH2D *h1 = new TH2D(histname.c_str(), info.c_str(), data[i].time.size(), data[i].time.front(), data[i].time.back(), 255, data[i].min[j],data[i].max[j]); 
			if (j == 0)
			{
				for(int k = 0; k < data[i].time.size(); k++)
					h1->Fill(data[i].time[k], data[i].ch1[k]);
			}
			else	
			{
				for(int k = 0; k < data[i].time.size(); k++)
					h1->Fill(data[i].time[k], data[i].ch2[k]);
			}

			h1->Write();
			delete h1;
		}
	}
	file->Write();	
	file->Close();
	delete file;


	return 0;
}


