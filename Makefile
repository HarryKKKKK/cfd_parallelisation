CXX      := g++
MPICXX   := mpicxx

CXXFLAGS := -O3 -std=c++17 -Iheader
OMPFLAGS := -fopenmp

# sources used by ALL variants
BASE_SRC := \
  src/init.cpp \
  src/physics.cpp \
  src/types.cpp \
  src/utils.cpp

SERIAL_SOLVER_SRC := src/solver.cpp
MPI_SOLVER_SRC    := src/mpi_solver.cpp

MAIN_SERIAL := main_output.cpp
MAIN_MPI    := main_mpi.cpp

.PHONY: all serial omp mpi clean

all: serial mpi omp

serial: serial.exe
mpi: mpi.exe
# omp: omp.exe

serial.exe: $(MAIN_SERIAL) $(BASE_SRC) $(SERIAL_SOLVER_SRC)
	$(CXX) $(CXXFLAGS) $^ -o $@

# omp.exe: $(MAIN_SERIAL) $(BASE_SRC) $(SERIAL_SOLVER_SRC)
# 	$(CXX) $(CXXFLAGS) $(OMPFLAGS) $^ -o $@

mpi.exe: $(MAIN_MPI) $(BASE_SRC) $(MPI_SOLVER_SRC)
	$(MPICXX) $(CXXFLAGS) $^ -o $@

clean:
	rm -f serial.exe mpi.exe omp.exe