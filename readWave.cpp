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
#include "fermidirac.h"

struct package
{
	std::vector<double> time;
	std::vector<double> ch1;
	std::vector<double> ch2;
	std::vector<double> max; // 0=ch1 1=ch2 ...
        std::vector<double> maxtime;
	std::vector<double> min;
	std::vector<double> ymulti;
	std::vector<double> filteredmax = {0,0};
	std::vector<double> meanthreshold = {0,0};
	std::vector<double> meanthresholdtime = {0,0};
        std::vector<double> fitmaximum = {0,0};
        std::vector<double> timefitmaximum = {0,0};
        std::vector<double> fitthreshold = {0,0};
        std::vector<double> timefitthreshold = {0,0};
	std::vector<double> thresholds = {13,20,40,60,79};
	std::vector<double> differences={0,0,0,0,0};
        std::vector<double> thresholdtimes= {0,0,0,0,0};
};


int findUnique(std::vector<double> v)
{
	std::sort(v.begin(), v.end());
	return std::unique(v.begin(), v.end()) - v.begin();
}

void scaleTime(std::vector<double> &v)
{
	for (unsigned int i = 0; i < v.size(); i++)
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
	if(filebase.back() == '/')
		filebase.pop_back();
	
	int num = numFiles(filebase);
        int numcutoff(0);		//used to set num to how many files have passed threshold
	//const int numtocheck = 10;	//debugging small batch
	if (num == -1)
	{
		std::cout << "Directory doesn't exist" << std::endl;
		return -1;
	}

        const double cutoff = 80; //maximum height of waveforms considered [mV]

	std::cout << "Parsing " << num << " files" << std::endl;	
	filebase += "/waveform"; 
	std::vector<package> data; 

	int posa(0),posb(0);
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
                const double timescale = 1E9;
		double mina(100),minb(100), maxa(-100), maxb(-100), tmaxa(0),tmaxb(0);
		for (std::string line; std::getline(file, line); )
		{
			std::stringstream ss(line);
			double d(0), a(0), b(0);

			ss >> d >> std::ws >> a >> std::ws >> b;
			/* all 3 channels are inverted */
			a = -a * scale;
			b = -b * scale;
			d *= timescale;	
			temp.time.push_back(d);
			temp.ch1.push_back(a);
			temp.ch2.push_back(b);

			if(a < mina)
			{
				mina = a;
				posa = temp.ch1.size();
			}
			if(a > maxa)
			{
				maxa = a;
                                tmaxa = d;
			}
			if(b < minb)
			{
				minb = b;
				posb = temp.ch2.size();
			}
			if(b > maxb)
			{
				maxb = b;
                                tmaxb = d;
			}

		}

                if((maxa <= cutoff) && (maxb <= cutoff))
                {
                        temp.max.push_back(maxa);
                        temp.max.push_back(maxb);
                        temp.maxtime.push_back(tmaxa);
                        temp.maxtime.push_back(tmaxb);
                        temp.min.push_back(mina);
                        temp.min.push_back(minb);
                        data.push_back(temp);
                        numcutoff++;
                }

                file.close();
                
	}
        std::cout << numcutoff << " files under cutoff threshold\n";
        
        
        num = numcutoff; // total number of waveforms under threshold, set above
        
	std::string filename = argv[1];
	filename.erase(0,3);
	filename = "data/" + filename + ".root";

	TFile *file  = new TFile(filename.c_str(), "RECREATE");
	//TDirectory *ch1 = file->mkdir("ch1");
	//TDirectory *ch2 = file->mkdir("ch2");
	TDirectory *stats = file->mkdir("stats");
	TDirectory *ch1graphs = file->mkdir("Ch1 pulse graphs");
	TDirectory *ch2graphs = file->mkdir("Ch2 pulse graphs");
	TDirectory *timing = file->mkdir("Timing");


	/* Find maximum and minimums among all of the waveforms */

	double  maxch1(0), maxch2(0), smallestmaxch1(180), smallestmaxch2(180), minch1(180), minch2(180); 
	for(int i = 0; i < num; i++)
	{
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

		/* sliding mean */
		for(int k = 0; k < 2; k++)
		{
			double tempmax{};
			std::deque<double> window(3,0.0);
			for(unsigned int j = 0; j < data[i].time.size(); j++)
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
                        
			for(unsigned int n = 0; n < data[i].thresholds.size(); n++)
			{	
				/* determine amplitude of threshold as a scale of maximum
				 ** amplitude */
				
				//data[i].meanthreshold[k] = data[i].threshold[k] * tempmax;
				
				/* find first location of threshold in time */

				for(unsigned int j = 0; j < data[i].time.size()-1; j++)
				{
					if(data[i].ch1[j] < data[i].thresholds[n] && data[i].thresholds[n] < data[i].ch1[j+1])
					{
						data[i].thresholdtimes[n] = data[i].time[j];
					}
				}
			}


}
}


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

        /* Plot graphs of waveforms below cutoff and fit them */
        /* Find the fit maximum and threshold positions in time */

        const double ttimeoffset = 0.5; //[ns]  This is how far past the maximum in time the fit will continue
        
	for(int p = 0; p < 2; p++)
        {        

                ch1graphs->cd();
                if(p ==1 )
                        ch2graphs->cd();

                for(int i = 0; i < num; i++)
                {
                        std::vector<double> x;
                        std::vector<double> y;
                        for(unsigned int j = 0; j < data[i].time.size(); j++)
                        {	
                                if(p == 0)
                                {
                                        y.push_back(data[i].ch1[j]);
                                        x.push_back(data[i].time[j]);
                                }
                                else
                                {
                                        y.push_back(data[i].ch2[j]);
                                        x.push_back(data[i].time[j]);
                                }
                        }
			// TODO Create subdirectories in the graphs directories
			// one for each threshold. Name the subdirectory for
			// the threshold

                        std::string name = "Fit Waveform" + std::to_string(i) + " PMT " + std::to_string(p+1);
			TF1 *f1 = new TF1("fermi_dirac_poly", fermi_dirac_parabola, -100, data[i].maxtime[p]+ttimeoffset, 5);
			f1->SetParNames("Alpha","t0","A", "c1", "c2");
                        
			TCanvas *c2 = new TCanvas(name.c_str(), "Waveform fitting");
                        TGraph * g2 = new TGraph(data[i].time.size(),&x[0],&y[0]);
                        f1->SetParameters(1, data[i].meanthresholdtime[p], data[i].filteredmax[p],0,0,0);
                        g2->Fit(f1, "Q", "",-100,data[i].maxtime[p]+ttimeoffset);	
                        g2->GetXaxis()->SetTitle("Time (ns)");
                        std::string title = "Waveform " + std::to_string(i) + " fitting";
                        g2->SetTitle(title.c_str());
                        g2->GetYaxis()->SetTitle("Amplitude (mV)");
                        g2->SetName(Form("g%d",i) );
                        g2->Draw("AP");
                        f1->Draw("SAME");
                        c2->Write();

                        /* compute maximum from fit and threshold crossing time */

                        data[i].fitmaximum[p] = f1->GetMaximum(-100,data[i].maxtime[p]+ttimeoffset);
                        data[i].timefitmaximum[p] = f1->GetMaximumX(-100, data[i].maxtime[p]+ttimeoffset);
                        data[i].fitthresholds[p] = data[i].threshold;
                        data[i].timefitthresholds[p] = f1->GetX(data[i].fitthresholds[p], -100, data[i].maxtime[p]);

                        data[i].difference = data[i].timefitmaximum[1]-data[i].timefitmaximum[0];
                        delete g2;
			delete c2;
                }
        }
        /* construct histogram of pulse heights */
	//TODO compute number of bins and extremes for these histograms	
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

	/* Scatterplot of amplitudes of PMT's */
	std::vector<double> max1;
	std::vector<double> max2;

	for(unsigned int i = 0; i < data.size(); i++)
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
	

	/* Histogram of threshold timing difference */	
	timing->cd();
	for(unsigned int j =0; j < data[0].thresholds.size(); j++)
	{
		std::string histname3 = "Timing Difference " + std::to_string(data[0].thresholds[j]);
		std::string info3 = "Timing difference of PMT's at" + std::to_string(data[0].thresholds[j]) + "mV";
		TH1D *h4 = new TH1D(histname3.c_str(), info3.c_str(), 100000, -100, 100);
		for(int i = 0; i < num; i++)
		{
			h4->Fill(data[i].differences[j]);
		}
		h4->GetXaxis()->SetTitle("Timing difference PMT1B-PMT1A (ns)");
		h4->GetYaxis()->SetTitle("Number of waveforms");
		h4->Write();
		delete h4;
	}
	

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


