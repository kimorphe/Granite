all: ascan bscan mean_awv
ascan: ascan.o fft.o wave1d.o voigt.o
	g++ -o ascan  ascan.o fft.o wave1d.o voigt.o
bscan: bscan.o fft.o wave1d.o voigt.o
	g++ -o bscan  bscan.o fft.o wave1d.o voigt.o
mean_awv: mean_awv.o fft.o wave1d.o
	g++ -o mean_awv mean_awv.o fft.o wave1d.o

mean_awv.o: mean_awv.cpp
	g++ -c mean_awv.cpp
bscan.o: bscan.cpp
	g++ -c  bscan.cpp
ascan.o: ascan.cpp
	g++ -c  ascan.cpp
fft.o: fft.cpp
	g++ -c fft.cpp
wave1d.o: wave1d.cpp
	g++ -c wave1d.cpp
voigt.o: voigt.cpp
	g++ -c voigt.cpp


