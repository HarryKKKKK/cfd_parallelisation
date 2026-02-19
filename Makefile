CXX      := g++
CXXFLAGS := -O3 -std=c++17 -Iheader
OMPFLAGS := -fopenmp

SRC  := $(wildcard src/*.cpp)
MAIN := main.cpp

.PHONY: all serial omp clean

all: serial omp

serial:
	$(CXX) $(CXXFLAGS) $(MAIN) $(SRC) -o serial.exe

omp:
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) $(MAIN) $(SRC) -o omp.exe

clean:
	rm -f serial.exe omp.exe
