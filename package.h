#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <deque>
#include <algorithm>
#include <numeric>
#include <cmath>

class package
{
	private:

	public:
	std::vector<double> time;
	std::vector<double> ch1;
	std::vector<double> ch2;
	std::vector<double> ch3;
	std::vector<double> ch4;
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
	void addLine(std::vector<double> &line)
	void findmaxfiltered(int width);
	void findmax();
	void findmin();



};
