
CXX = clang++ -O 
# CXX = g++-mp-10 -O
# CXX += -g -fsanitize=address
CXX += -Wall

# CXX = g++9 -O3

FLAGS = -std=c++17 -I/opt/local/include -I/usr/local/include
LIBS = -L/opt/local/lib -L/usr/local/lib -lfftw3 -lsndfile -lportaudio

SRC = ui.cc util.cc demod.cc fft.cc snd.cc taps.cc bench.cc mod.cc

pskons: $(SRC) demod.h fft.h snd.h util.h bench.h 
	$(CXX) $(FLAGS) $(SRC) -o pskons $(LIBS) -pthread

db: $(SRC) demod.h fft.h
	$(CXX) -g -fsanitize=address $(FLAGS) $(SRC) -o db $(LIBS) -pthread

clean:
	rm -f pskons db

