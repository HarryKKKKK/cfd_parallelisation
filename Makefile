CXX      := g++
MPICXX   := mpicxx

CXXFLAGS := -O3 -std=c++17 -Iheader
OMPFLAGS := -fopenmp

COMMON_SRC := \
  src/init.cpp \
  src/physics.cpp \
  src/types.cpp \
  src/utils.cpp \
  src/solver.cpp

# MPI-only sources (add these files)
MPI_SRC := \
  src/mpi_solver.cpp 

MAIN_SERIAL := main_output.cpp
MAIN_MPI    := main_mpi.cpp

.PHONY: all serial omp mpi clean

all: serial omp mpi

serial: serial.exe
omp: omp.exe
mpi: mpi.exe

serial.exe: $(MAIN_SERIAL) $(COMMON_SRC)
	$(CXX) $(CXXFLAGS) $(MAIN_SERIAL) $(COMMON_SRC) -o $@

# omp.exe: $(MAIN_SERIAL) $(COMMON_SRC)
# 	$(CXX) $(CXXFLAGS) $(OMPFLAGS) $(MAIN_SERIAL) $(COMMON_SRC) -o $@

# mpi.exe: $(MAIN_MPI) $(COMMON_SRC) $(MPI_SRC)
# 	$(MPICXX) $(CXXFLAGS) $(MAIN_MPI) $(COMMON_SRC) $(MPI_SRC) -o $@

clean:
	rm -f serial.exe omp.exe mpi.exe
