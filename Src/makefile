all: ascan bscan mean_awv win_test tmax
ascan: ascan.o fft.o wave1d.o voigt.o
	g++ -o ascan  ascan.o fft.o wave1d.o voigt.o
	cp ascan ../
bscan: bscan.o fft.o wave1d.o voigt.o
	g++ -o bscan  bscan.o fft.o wave1d.o voigt.o
	cp bscan ../
mean_awv: mean_awv.o fft.o wave1d.o
	g++ -o mean_awv mean_awv.o fft.o wave1d.o
	cp mean_awv ../
win_test: win_test.o fft.o wave1d.o
	g++ -o win_test  win_test.o fft.o wave1d.o
	cp win_test ../
tmax: tmax.o fft.o wave1d.o
	g++ -o tmax  tmax.o fft.o wave1d.o
	cp tmax ../

tmax.o: tmax.cpp
	g++ -c tmax.cpp
win_test.o: win_test.cpp
	g++ -c win_test.cpp
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


