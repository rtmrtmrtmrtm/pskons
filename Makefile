
CXX = c++ -O 

FLAGS = -std=c++17 -I/opt/local/include -I/usr/local/include
LIBS = -L/opt/local/lib -L/usr/local/lib -lfftw3 -lsndfile -lportaudio

# CXX += -DAIRSPYHF
# LIBS += -lairspyhf -lliquid -lusb

SRC = ui.cc util.cc demod.cc fft.cc snd.cc taps.cc bench.cc mod.cc

pskons: $(SRC) demod.h fft.h snd.h util.h bench.h 
	$(CXX) $(FLAGS) $(SRC) -o pskons $(LIBS) -pthread

clean:
	rm -f pskons 

