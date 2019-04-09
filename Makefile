ROOT=`root-config --libs --glibs --incdir'
FILES=readwave.cpp countfiles.cpp
INCROOT="/usr/include/root/ROOT"
INCROOT2="/usr/include/root"

default: 
	g++ -g `root-config --libs --glibs` -I/cvmfs/sft.cern.ch/lcg/releases/LCG_94/ROOT/6.14.04/x86_64-slc6-gcc62-opt/include -I/cvmfs/sft.cern.ch/lcg/releases/LCG_94/ROOT/6.14.04/x86_64-slc6-gcc62-opt/include/ROOT readWave.cpp countfiles.cpp fermidirac.cpp -o readWave

fast: 
	g++ -Ofast -march=native `root-config --libs --glibs` -I/cvmfs/sft.cern.ch/lcg/releases/LCG_94/ROOT/6.14.04/x86_64-slc6-gcc62-opt/include -I/cvmfs/sft.cern.ch/lcg/releases/LCG_94/ROOT/6.14.04/x86_64-slc6-gcc62-opt/include/ROOT readWave.cpp countfiles.cpp fermidirac.cpp -o readWave

run:
	./readWave ../PMT_led_trigger_waveforms 
	#root PMT_fit.root

clean:
	rm -f readWave *.png
