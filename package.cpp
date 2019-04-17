#include "package.h"

void package::addLine(std::vector<double> &line)
{	
	for(int i =0; i < line.size(); i++)
	{
		switch(i)
		{
			case 0: time.push_back(line[0]);
				break;
			case 1: ch1.push_back(line[1]);
				break;
			case 2: ch2.push_back(line[2]);
				break;
			case 3: ch3.push_back(line[3]);
				break;
			case 4: ch4.push_back(line[4]);
				break;
			default: break;
		}
	}
}
