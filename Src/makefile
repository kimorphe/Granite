all: ascan bscan
ascan: ascan.o fft.o wave1d.o
	g++ -o ascan  ascan.o fft.o wave1d.o
bscan: bscan.o fft.o wave1d.o
	g++ -o bscan  bscan.o fft.o wave1d.o

bscan.o: bscan.cpp
	g++ -c  bscan.cpp
ascan.o: ascan.cpp
	g++ -c  ascan.cpp
fft.o: fft.cpp
	g++ -c fft.cpp
wave1d.o: wave1d.cpp
	g++ -c wave1d.cpp


