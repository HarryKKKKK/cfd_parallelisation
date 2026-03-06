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

MAIN_SERIAL      := main_output.cpp
MAIN_MPI_OUTPUT  := main_mpi.cpp
MAIN_MPI_STRONG  := main_strong.cpp
MAIN_MPI_WEAK    := main_weak.cpp

.PHONY: all serial omp mpi mpi_strong mpi_weak clean

all: serial omp mpi mpi_strong mpi_weak

# ---- serial / omp use the serial solver.cpp ----
serial: serial.exe
omp: omp.exe

serial.exe: $(MAIN_SERIAL) $(BASE_SRC) $(SERIAL_SOLVER_SRC)
	$(CXX) $(CXXFLAGS) $^ -o $@

omp.exe: $(MAIN_SERIAL) $(BASE_SRC) $(SERIAL_SOLVER_SRC)
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) $^ -o $@

# ---- mpi variants use mpi_solver.cpp (do NOT link solver.cpp) ----
mpi: mpi.exe
mpi_strong: mpi_strong.exe
mpi_weak: mpi_weak.exe

mpi.exe: $(MAIN_MPI_OUTPUT) $(BASE_SRC) $(MPI_SOLVER_SRC)
	$(MPICXX) $(CXXFLAGS) $^ -o $@

mpi_strong.exe: $(MAIN_MPI_STRONG) $(BASE_SRC) $(MPI_SOLVER_SRC)
	$(MPICXX) $(CXXFLAGS) $^ -o $@

mpi_weak.exe: $(MAIN_MPI_WEAK) $(BASE_SRC) $(MPI_SOLVER_SRC)
	$(MPICXX) $(CXXFLAGS) $^ -o $@

clean:
	rm -f serial.exe omp.exe mpi.exe mpi_strong.exe mpi_weak.exe