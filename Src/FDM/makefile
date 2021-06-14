all: fdm

fdm:: domain.o main.o field.o inwv.o cpbank.o recs.o
	g++ -o fdm domain.o main.o field.o inwv.o cpbank.o recs.o

inwv.o: inwv.cpp
	g++ -c inwv.cpp
field.o:field.cpp
	g++ -c field.cpp
domain.o: domain.cpp
	g++ -c domain.cpp
main.o: main.cpp
	g++ -c main.cpp

cpbank.o: cpbank.cpp
	g++ -c cpbank.cpp
recs.o: recs.cpp
	g++ -c recs.cpp

