outname = readWave
ROOT != root-config --libs --glibs
FILES=readWave.cpp countfiles.cpp fermidirac.cpp
INCROOT != root-config --incdir

default: 
	g++ -g -Wall $(ROOT) -I$(INCROOT) $(FILES) -o$(outname)

fast: 
	g++ -Ofast -march=native $(ROOT) -I$(INCROOT) $(FILES) -o$(outname)

run:
	./readWave ../PMT_led_trigger_waveforms 
	#root PMT_fit.root

clean:
	rm -f readWave *.png
