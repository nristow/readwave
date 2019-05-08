#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <deque>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <map>
#include <root/TGraph.h>
#include <root/TGraphErrors.h>
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
#include <root/TTree.h>
#include "countfiles.h"
#include "fermidirac.h"

// TODO read in file to set these parameters
class datarun
{
	public:
	bool invertwaveforms = true;
	std::vector<std::string> pmtnumber = {"1","2"};
	std::vector<std::string> pmtmanufacturer = {"OnSemi", "OnSemi"};
	std::vector<std::string> pmtpartnumber = {"C-Series","C-Series"};
	std::vector<int> voltages = {30,30};
	std::map<int,std::string> labels;
	std::string experimenttype = "Cosmics";
	std::string trigger = "Discriminator";
	double uppercutoffthreshold = 90;
	//std::vector<double> thresholds = {20,30,40,50,60,79};
	std::vector<double> thresholds = {5,6,8,10,12,9};
	int runnumber = 17;
	std::vector<std::string> sthresholds;
	int slidingwindowwidth = 7;
	double timeoffset = 3; // [ns] Time after maximum for fit to continue
	std::vector<double> fitmean = {0,0,0,0,0,0};
	std::vector<double> fitsigma = {0,0,0,0,0,0};
	std::vector<double> fitmaxX = {0,0,0,0,0,0};
	datarun();
	friend std::ostream& operator<<(std::ostream& os, const datarun& dr);
};

std::ostream& operator<<(std::ostream&os, const datarun& dr)
{
	os << "Run Number: " << dr.runnumber << std::endl;
	os << "Experiment type: " << dr.experimenttype << std::endl;
	os << "Trigger: " << dr.trigger << std::endl;
	os << "PMT " << dr.pmtnumber[0] << "A " << ": " << dr.pmtmanufacturer[0] << " " << dr.pmtpartnumber[0] << " at " << dr.voltages[0] << "V" << std::endl;
	os << "PMT " << dr.pmtnumber[1] << "A " << ": " << dr.pmtmanufacturer[1] << " " << dr.pmtpartnumber[1] << " at " << dr.voltages[1] << "V" << std::endl;
	os << "Time offset for waveform fit: " << dr.timeoffset << " ns" << std::endl;
	os << "Sliding window width: " << dr.slidingwindowwidth << std::endl;
	os << "Upper cutoff threshold: " << dr.uppercutoffthreshold << " mv" << std::endl;
	os << "Waveform inverted?: " << dr.invertwaveforms << std::endl;
       return os;	
}

datarun::datarun()
{
	labels.insert(std::pair<int,std::string>(0,"A"));
	labels.insert(std::pair<int,std::string>(1,"B"));
	
	for(auto&i : thresholds)
	{
		sthresholds.push_back(std::to_string(thresholds[i]));
	}

	fitmean.resize(thresholds.size(), 0);
	fitsigma.resize(thresholds.size(), 0);
	fitmaxX.resize(thresholds.size(), 0);
}
		
class package
{
	public:
		std::vector<double> time;
		std::vector<double> ch1;
		std::vector<double> ch2;
		std::vector<double> max; // 0=ch1 1=ch2 ... maximum of all datapoints
		std::vector<double> maxtime;
		std::vector<double> min;
		std::vector<double> ymulti;
		std::vector<double> filteredmax = {0,0}; //maximum of the sliding window function. 
		std::vector<double> halffilteredmax = {0,0};
		std::vector<double> hfiltmaxtime = {0,0};
		std::vector<double> differences = {0,0,0,0,0,0};
		std::vector<double> thresholdtimesch1 = {0,0,0,0,0,0};
		std::vector<double> thresholdtimesch2 = {0,0,0,0,0,0};
		std::vector<double> fitmax = {0,0};

		void setfilteredmax(double value, int channel);
		void setnumberofthresholds(unsigned int num);
};

void package::setnumberofthresholds(unsigned int num)
{
	differences.resize(num,0);
	thresholdtimesch1.resize(num,0);
	thresholdtimesch2.resize(num,0);
}

void package::setfilteredmax(double value, int channel)
{
	filteredmax[channel] = value;
	halffilteredmax[channel] = value/2;
}


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
	datarun thisrun;
	if (argc != 2)
	{
		std::cout << "Usage: readwave <path/to/files>" << std::endl;
		return 0;
	}

	std::string filebase=argv[1];
	if(filebase.back() == '/')
		filebase.pop_back();
	
	int num = numFiles(filebase);
	if (num == -1)
	{
		std::cout << "Directory doesn't exist" << std::endl;
		return -1;
	}

        //const double cutoff = 80; //maximum height of waveforms considered [mV]

	std::cout << "Parsing " << num << " files" << std::endl;	
	filebase += "/waveform"; 
	std::vector<package> data; 

	bool oldheader = false;
	int numcutoff = 0;	
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
		temp.setnumberofthresholds(thisrun.thresholds.size());

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
			if(thisrun.invertwaveforms == true)
			{
				a = -a * scale;
				b = -b * scale;
				d *= timescale;	
			}
			else
			{
		
				a = a * scale;
				b = b * scale;
				d *= timescale;	
			}
			temp.time.push_back(d);
			temp.ch1.push_back(a);
			temp.ch2.push_back(b);

			if(a < mina)
			{
				mina = a;
			}
			if(a > maxa)
			{
				maxa = a;
                                tmaxa = d;
			}
			if(b < minb)
			{
				minb = b;
			}
			if(b > maxb)
			{
				maxb = b;
                                tmaxb = d;
			}

		}

                if((maxa <= thisrun.uppercutoffthreshold) && (maxb <= thisrun.uppercutoffthreshold))
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
        std::cout << numcutoff << " files under cutoff threshold of " << thisrun.uppercutoffthreshold << " mv" << std::endl;
         
        num = numcutoff; // total number of waveforms under threshold, set above
        
	std::string filename = argv[1];
	filename.erase(0,3);
	std::string runfile = filename;
	filename = "data/" + filename + ".root";
	runfile = "data/" + runfile + ".txt";
	std::fstream of;
	of.open(runfile, std::fstream::trunc | std::fstream::out);
	of << thisrun;
	of.close();

	TFile *file  = new TFile(filename.c_str(), "RECREATE");
	//TDirectory *ch1 = file->mkdir("ch1");
	//TDirectory *ch2 = file->mkdir("ch2");
	TDirectory *stats = file->mkdir("stats");
	TDirectory *ch1graphs = file->mkdir("Ch1_pulse_graphs");
	TDirectory *ch2graphs = file->mkdir("Ch2_pulse_graphs");
	TDirectory *timing = file->mkdir("Timing");
	
	/* Save run parameters to root file */
	TTree tree("setup", "Run information");
	tree.Branch("PMT Manufacturers",&thisrun.pmtmanufacturer);
	tree.Branch("PMT Part Numbers",&thisrun.pmtpartnumber);
	tree.Branch("PMT numbers",&thisrun.pmtnumber);
	tree.Branch("Voltages",&thisrun.voltages);
	tree.Branch("Experiment Type",&thisrun.experimenttype);
	tree.Branch("Experiment Trigger",&thisrun.trigger);
	tree.Branch("Window width",&thisrun.slidingwindowwidth);
	tree.Branch("Fit time offset",&thisrun.timeoffset);




	std::cout << "reading waveforms" << std::endl;

	for(int i = 0; i < num; i++)
	{

		/* sliding mean */
		for(int k = 0; k < 2; k++)
		{
			double tempmax{};
			std::deque<double> window(thisrun.slidingwindowwidth,0.0);
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
                        data[i].setfilteredmax(tempmax, k); 
			

			//this half time of maximum not needed 
			for(unsigned int j = 0; j < data[i].time.size()-1; j++)
			{
				if(data[i].ch1[j] < data[i].halffilteredmax[k]	&& data[i].halffilteredmax[k] < data[i].ch1[j+1])
				{
					data[i].hfiltmaxtime[k] = data[i].time[j];
				}
			}
		}                       
	}

        /* Plot graphs of waveforms below cutoff and fit them */
        /* Find the fit maximum and threshold positions in time */

	std::cout << "Creating waveform graphs" << std::endl;
        
	for(int channel = 0; channel < 2; channel++)
        {        

                ch1graphs->cd();
                if(channel ==1 )
                        ch2graphs->cd();

                for(int numwaveform = 0; numwaveform < num; numwaveform++)
                {
                        std::vector<double> x;
                        std::vector<double> y;
                        for(unsigned int j = 0; j < data[numwaveform].time.size(); j++)
                        {	
                                if(channel == 0)
                                {
                                        y.push_back(data[numwaveform].ch1[j]);
                                        x.push_back(data[numwaveform].time[j]);
                                }
                                else
                                {
                                        y.push_back(data[numwaveform].ch2[j]);
                                        x.push_back(data[numwaveform].time[j]);
                                }
                        }

                        std::string name = "Fit_Waveform" + std::to_string(numwaveform) + "_PMT_" + thisrun.pmtnumber[channel] + thisrun.labels[channel];
			TF1 *f1 = new TF1("fermi_dirac_poly", fermi_dirac_parabola, -100, data[numwaveform].maxtime[channel]+thisrun.timeoffset, 5);
			f1->SetParNames("Alpha","t0","A", "c1", "c2");
                        
			TCanvas *c2 = new TCanvas(name.c_str(), "Waveform_fitting");
			TGraph * g2 = new TGraph(data[numwaveform].time.size(),&x[0],&y[0]);
                        f1->SetParameters(1, data[numwaveform].maxtime[channel], data[numwaveform].filteredmax[channel],0,0);
                        int status = g2->Fit(f1, "Q", "",-100,data[numwaveform].maxtime[channel]+thisrun.timeoffset);	
                        g2->GetXaxis()->SetTitle("Time (ns)");
                        std::string title = "Waveform " + std::to_string(numwaveform) + " fitting";
                        g2->SetTitle(title.c_str());
			g2->SetMarkerStyle(7);
                        g2->GetYaxis()->SetTitle("Amplitude (mV)");
                        g2->SetName(Form("g%d",numwaveform) );
                        g2->Draw("AP");
                        f1->Draw("SAME");
                        c2->Write();

			if(status < 0 || (status > 0 && status != 4))
			{
				std::cout << "Fit of channel " << channel << " waveform " << numwaveform << " returned " << status << std::endl;
			}
                        /* compute threshold crossing time for all thresholds*/
			/* make differences */

			data[numwaveform].fitmax[channel] = f1->GetMaximum(-100,data[numwaveform].maxtime[channel]+thisrun.timeoffset);

			for(unsigned int thresholdnum = 0; thresholdnum < thisrun.thresholds.size(); thresholdnum++)
			{
				if(channel == 0)
				{
					data[numwaveform].thresholdtimesch1[thresholdnum] = f1->GetX(thisrun.thresholds[thresholdnum], -100, data[numwaveform].maxtime[channel]);
				}
				else
				{
					data[numwaveform].thresholdtimesch2[thresholdnum] = f1->GetX(thisrun.thresholds[thresholdnum], -100, data[numwaveform].maxtime[channel]);

				}

			}

			for(unsigned int thresholdnum = 0; thresholdnum < thisrun.thresholds.size(); thresholdnum++)
			{	
				data[numwaveform].differences[thresholdnum] = data[numwaveform].thresholdtimesch2[thresholdnum] - data[numwaveform].thresholdtimesch1[thresholdnum];
			}
                        delete g2;
			delete c2;
			delete f1;
                }
        }


	std::cout << "Pulse height histograms" << std::endl;

        /* construct histogram of pulse heights */
	stats->cd();
	std::string histname = "PH2";
	std::string info = "Pulse Height for PMT-" + thisrun.pmtnumber[1] + thisrun.labels[1] + ", " + thisrun.pmtmanufacturer[1] + " " + thisrun.pmtpartnumber[1] + " at " + std::to_string(thisrun.voltages[1]) + "V";
	TH1D *h1 = new TH1D(histname.c_str(), info.c_str(), 3000, -1, 2000);
	for(int i = 0; i < num ; i++)
	{
		h1->Fill(data[i].fitmax[1]); //use fit  maximums
	}
	h1->GetXaxis()->SetTitle("Maximum Waveform Height (mV)");
	h1->GetYaxis()->SetTitle("Number of waveforms");
	h1->Fit("gaus");
	h1->Write();
	delete h1;
	
	histname = "PH1";
	info = "Pulse Height for PMT-" + thisrun.pmtnumber[0] + thisrun.labels[0] + ", " + thisrun.pmtmanufacturer[0] + " " + thisrun.pmtpartnumber[0] + " at " + std::to_string(thisrun.voltages[0]) + "V";
	TH1D *h2 = new TH1D(histname.c_str(), info.c_str(), 3000, -1, 2000);
	for(int i = 0; i < num ; i++)
	{
		h2->Fill(data[i].fitmax[0]);
	}

	h2->GetXaxis()->SetTitle("Maximum Waveform Height (mV)");
	h2->GetYaxis()->SetTitle("Number of waveforms");
	h2->Fit("gaus");
	h2->Write();
	delete h2;

	std::cout << "Scatter plots" << std::endl;

	/* Scatterplot of amplitudes of PMT's */
	std::vector<double> max1;
	std::vector<double> max2;

	for(unsigned int i = 0; i < data.size(); i++)
	{
		max1.push_back(data[i].fitmax[0]);
		max2.push_back(data[i].fitmax[1]);
	}
	
	TCanvas *c1 = new TCanvas("XY Plot", "PMT maximums XY plot");
	TGraph * g1 = new TGraph(data.size(),&max1[0],&max2[0]);
	std::string pmta = "PMT-" + thisrun.pmtnumber[0] + thisrun.labels[0] + " maximum [mV]";
	std::string pmtb = "PMT-" + thisrun.pmtnumber[1] + thisrun.labels[1] + " maximum [mV]";
	g1->GetXaxis()->SetTitle(pmta.c_str());
	g1->GetXaxis()->SetTitle("PMT Timing Difference");
	g1->GetYaxis()->SetTitle(pmtb.c_str());
        g1->SetMarkerStyle(7);
	g1->SetDrawOption("AP");
	g1->Draw("AP");
	c1->Write();
	
	delete c1;
	delete g1;

	/* Histogram of threshold timing difference */	

	std::cout << "Threshold Timing plots" << std::endl;
	
	timing->cd();
	for(unsigned int j =0; j < thisrun.thresholds.size(); j++)
	{
		std::string histname3 = "Timing Difference " + std::to_string(thisrun.thresholds[j]) + "mV";
		std::string info3 = "Timing difference of PMT's at" + std::to_string(thisrun.thresholds[j]) + "mV";
		TH1D *h4 = new TH1D(histname3.c_str(), info3.c_str(), 400, -100, 100);
		for(int i = 0; i < num; i++)
		{
			h4->Fill(data[i].differences[j]);
		}
		h4->GetXaxis()->SetTitle("Timing difference PMT1B-PMT1A (ns)");
		h4->GetYaxis()->SetTitle("Number of waveforms");
		h4->Fit("gaus");
		h4->Write();
		TF1 *tfff = h4->GetFunction("gaus");
		thisrun.fitmean[j] = tfff->GetParameter(1); 
		thisrun.fitsigma[j] = tfff->GetParameter(2);
		thisrun.fitmaxX[j] = tfff->GetMaximumX();
		delete tfff;
		delete h4;
	}
	

	/* graph of the mean from the threshold plots */
	std::vector<double> zero(0,6);
	TCanvas *c10 = new TCanvas("Mean Plot", "PMT Timing means");
	TGraph *g10 = new TGraph(thisrun.fitmean.size(),&thisrun.thresholds[0],&thisrun.fitmaxX[0]);
	g10->GetXaxis()->SetTitle("Threshold");
	g10->GetYaxis()->SetTitle("Mean of timing fit");
	g10->SetMarkerStyle(7);
	g10->SetMarkerColor(4);
	g10->SetTitle("Mean of difference thresholds");
	g10->SetDrawOption("AP");
	g10->Draw("AP");
	c10->Write();
	delete c10;
	delete g10;

	/* Scatterplot of mean and sigmas from threshold plots above */
	TCanvas *c5 = new TCanvas("XY Plot Timing", "PMT Timing sigmas vs means");
	TGraph * g5 = new TGraph(thisrun.fitmean.size(),&thisrun.fitmean[0],&thisrun.fitsigma[0]);
	g5->GetXaxis()->SetTitle("Means");
	g5->GetYaxis()->SetTitle("Sigmas");
	g5->SetMarkerStyle(7);
	g5->SetMarkerColor(4);
	g5->SetDrawOption("AP");
	g5->Draw("AP");
	c5->Write();
	
	delete c5;
	delete g5;

	file->Write();	
	file->Close();
	delete file;


	return 0;
}


