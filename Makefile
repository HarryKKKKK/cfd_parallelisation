CXX := g++
CXXFLAGS := -O3 -std=c++17 -Iheader
OMPFLAGS := -fopenmp

COMMON_SRC := \
  src/init.cpp \
  src/physics.cpp \
  src/types.cpp \
  src/utils.cpp

MAIN_SRC := main.cpp
SOLVER_SRC := src/solver.cpp

.PHONY: all serial omp clean

all: serial omp

serial: serial.exe
omp: omp.exe

serial.exe: $(MAIN_SRC) $(COMMON_SRC) $(SOLVER_SRC)
	$(CXX) $(CXXFLAGS) $^ -o $@

omp.exe: $(MAIN_SRC) $(COMMON_SRC) $(SOLVER_SRC)
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) $^ -o $@

clean:
	rm -f serial.exe omp.exe
