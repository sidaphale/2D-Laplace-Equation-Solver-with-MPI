all : Run1

Run1 : parallel.cpp grid.hpp cgkrylov.hpp
	mpiicpc -Wall -std=c++11 -o Run1 parallel.cpp

clean :
	rm *.out Run1
