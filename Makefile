CXX      := g++
MPICXX   := mpicxx

CXXFLAGS := -O3 -std=c++17 -Iheader
OMPFLAGS := -fopenmp

# =========================
# Common sources
# =========================
BASE_SRC := \
  src/init.cpp \
  src/physics.cpp \
  src/types.cpp \
  src/utils.cpp

SERIAL_SOLVER_SRC := src/solver.cpp
MPI_SOLVER_SRC    := src/mpi_solver.cpp

# =========================
# Main files
# =========================
MAIN_OUTPUT        := main_output.cpp
MAIN_OUTPUT_MPI    := main_output_mpi.cpp
MAIN_BASELINE      := main_baseline.cpp
MAIN_BASELINE_MPI  := main_baseline_mpi.cpp

# =========================
# Phony targets
# =========================
.PHONY: all output scaling clean \
        serial_output omp_output mpi_output \
        serial_scaling omp_scaling mpi_scaling

# Build everything
all: output scaling

# =========================
# Group targets
# =========================
output: serial_output.exe omp_output.exe mpi_output.exe
scaling: serial_scaling.exe omp_scaling.exe mpi_scaling.exe

serial_output: serial_output.exe
omp_output: omp_output.exe
mpi_output: mpi_output.exe

serial_scaling: serial_scaling.exe
omp_scaling: omp_scaling.exe
mpi_scaling: mpi_scaling.exe

# =========================
# Output executables
# =========================
serial_output.exe: $(MAIN_OUTPUT) $(BASE_SRC) $(SERIAL_SOLVER_SRC)
	$(CXX) $(CXXFLAGS) $^ -o $@

omp_output.exe: $(MAIN_OUTPUT) $(BASE_SRC) $(SERIAL_SOLVER_SRC)
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) $^ -o $@

mpi_output.exe: $(MAIN_OUTPUT_MPI) $(BASE_SRC) $(MPI_SOLVER_SRC)
	$(MPICXX) $(CXXFLAGS) $^ -o $@

# =========================
# Scaling / baseline executables
# =========================
serial_scaling.exe: $(MAIN_BASELINE) $(BASE_SRC) $(SERIAL_SOLVER_SRC)
	$(CXX) $(CXXFLAGS) $^ -o $@

omp_scaling.exe: $(MAIN_BASELINE) $(BASE_SRC) $(SERIAL_SOLVER_SRC)
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) $^ -o $@

mpi_scaling.exe: $(MAIN_BASELINE_MPI) $(BASE_SRC) $(MPI_SOLVER_SRC)
	$(MPICXX) $(CXXFLAGS) $^ -o $@

# =========================
# Clean
# =========================
clean:
	rm -f \
	  serial_output.exe omp_output.exe mpi_output.exe \
	  serial_scaling.exe omp_scaling.exe mpi_scaling.exe