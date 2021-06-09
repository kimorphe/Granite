all: fdm

fdm:: domain.o main.o field.o inwv.o
	g++ -o fdm domain.o main.o field.o inwv.o

inwv.o: inwv.cpp
	g++ -c inwv.cpp
field.o:field.cpp
	g++ -c field.cpp
domain.o: domain.cpp
	g++ -c domain.cpp
main.o: main.cpp
	g++ -c main.cpp

