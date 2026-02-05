CXX      = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra

SRC = serial.cpp \
      init.cpp \
      physics.cpp \
      solver.cpp \
      types.cpp \
      utils.cpp

TARGET = serial

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET)

clean:
	rm -f $(TARGET)
