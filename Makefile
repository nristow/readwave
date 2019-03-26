ROOT=`root-config --libs --glibs'
INCROOT="/usr/include/root/ROOT"
INCROOT2="/usr/include/root"

default: 
	g++ -g `root-config --libs --glibs ` -I/usr/include/root -I/usr/include/root/ROOT readWave.cpp -o readWave

run:
	#`readWave ../waveform_1553105709/` 
	root PMT_fit.root

clean:
	rm -f readWave *.png
