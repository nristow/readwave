ROOT=`root-config --libs --glibs'
FILES=readwave.cpp countfiles.cpp
INCROOT="/usr/include/root/ROOT"
INCROOT2="/usr/include/root"

default: 
	g++ -g `root-config --libs --glibs ` -I/usr/include/root -I/usr/include/root/ROOT readWave.cpp countfiles.cpp -o readWave

run:
	./readWave ../PMT_led_trigger_waveforms 
	#root PMT_fit.root

clean:
	rm -f readWave *.png
