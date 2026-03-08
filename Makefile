CXX      := g++
MPICXX   := mpicxx

CXXFLAGS :=  -O3 -DNDEBUG -march=native -std=c++17 -Iheader
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
MAIN_OUTPUT            := main_output.cpp
MAIN_OUTPUT_MPI        := main_output_mpi.cpp
MAIN_BASE              := main_baseline.cpp
MAIN_BASE_MPI          := main_baseline_mpi.cpp
MAIN_BASE_SCALING      := main_base_scaling.cpp
MAIN_BASE_SCALING_MPI  := main_base_scaling_mpi.cpp

# =========================
# Phony targets
# =========================
.PHONY: all output base scaling clean \
        serial_output omp_output mpi_output \
        serial_base omp_base mpi_base \
        serial_scaling omp_scaling mpi_scaling

# Build everything
all: output base scaling

# =========================
# Group targets
# =========================
output: serial_output.exe omp_output.exe mpi_output.exe
base:   serial_base.exe   omp_base.exe   mpi_base.exe
scaling: serial_scaling.exe omp_scaling.exe mpi_scaling.exe

serial_output: serial_output.exe
omp_output: omp_output.exe
mpi_output: mpi_output.exe

serial_base: serial_base.exe
omp_base: omp_base.exe
mpi_base: mpi_base.exe

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
# Baseline executables
# =========================
serial_base.exe: $(MAIN_BASE) $(BASE_SRC) $(SERIAL_SOLVER_SRC)
	$(CXX) $(CXXFLAGS) $^ -o $@

omp_base.exe: $(MAIN_BASE) $(BASE_SRC) $(SERIAL_SOLVER_SRC)
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) $^ -o $@

mpi_base.exe: $(MAIN_BASE_MPI) $(BASE_SRC) $(MPI_SOLVER_SRC)
	$(MPICXX) $(CXXFLAGS) $^ -o $@

# =========================
# Scaling executables
# =========================
serial_scaling.exe: $(MAIN_BASE_SCALING) $(BASE_SRC) $(SERIAL_SOLVER_SRC)
	$(CXX) $(CXXFLAGS) $^ -o $@

omp_scaling.exe: $(MAIN_BASE_SCALING) $(BASE_SRC) $(SERIAL_SOLVER_SRC)
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) $^ -o $@

mpi_scaling.exe: $(MAIN_BASE_SCALING_MPI) $(BASE_SRC) $(MPI_SOLVER_SRC)
	$(MPICXX) $(CXXFLAGS) $^ -o $@

# =========================
# Clean
# =========================
clean:
	rm -f \
	  serial_output.exe omp_output.exe mpi_output.exe \
	  serial_base.exe omp_base.exe mpi_base.exe \
	  serial_scaling.exe omp_scaling.exe mpi_scaling.exe