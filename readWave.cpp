#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <root/TGraph.h>
#include <root/TApplication.h>

int main(int argc, char* argv[])
{
	if (argc < 2)
	{
		std::cout << "Usage: readwave <filename>" << std::endl;
		return 0;
	}

	std::ifstream file;
		
	file.open(argv[1]);
	
	std::vector<double> time;
	std::vector<double> ch1;
	std::vector<double> ch2;
	std::vector<double> ch3;


	for(int i = 0; i < 7; i++)
	{	
		std::string line;
		std::getline(file, line);
	}

	for (std::string line; std::getline(file, line); )
	{
		std::stringstream ss(line);
		double d, a, b, c;
		ss >> d >> std::ws >> a >> std::ws >> b >> std::ws >> c;
		time.push_back(d);
		ch1.push_back(a);
		ch2.push_back(b);
		ch3.push_back(c);
	}
	file.close();

	TApplication* rootapp = new TApplication("PMT timing", &argc, argv);
	TGraph *gr1 = new TGraph(time.size(), &time[0], &ch1[0]);
	gr1->Draw();
	rootapp->Run();	

	return 0;
}
