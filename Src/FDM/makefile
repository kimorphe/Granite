all: fdm

fdm:: domain.o main.o field.o inwv.o cpbank.o recs.o 
	g++ -o fdm domain.o main.o field.o inwv.o cpbank.o recs.o

inwv.o: inwv.cpp fdm2d.h
	g++ -c inwv.cpp
field.o:field.cpp fdm2d.h
	g++ -c field.cpp
domain.o: domain.cpp fdm2d.h
	g++ -c domain.cpp
main.o: main.cpp fdm2d.h
	g++ -c main.cpp

cpbank.o: cpbank.cpp fdm2d.h
	g++ -c cpbank.cpp
recs.o: recs.cpp fdm2d.h
	g++ -c recs.cpp

