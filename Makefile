ROOT=`root-config --libs --glibs'
FILES=readwave.cpp countfiles.cpp
INCROOT="/usr/include/root/ROOT"
INCROOT2="/usr/include/root"

default: 
	g++ -g -Wall `root-config --libs --glibs ` -I/usr/include/root -I/usr/include/root/ROOT readWave.cpp countfiles.cpp fermidirac.cpp -o readWave

fast: 
	g++ -Ofast -march=native `root-config --libs --glibs ` -I/usr/include/root -I/usr/include/root/ROOT readWave.cpp countfiles.cpp fermidirac.cpp -o readWave

run:
	./readWave ../PMT_led_trigger_waveforms 
	#root PMT_fit.root

clean:
	rm -f readWave *.png
