all: fdm

fdm:: domain.o main.o
	g++ -o fdm domain.o main.o

domain.o: domain.cpp
	g++ -c domain.cpp
main.o: main.cpp
	g++ -c main.cpp
