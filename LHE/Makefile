CXX = g++
CXXFLAGS = -std=c++17 -Wall -Iinclude -fopenmp
LDFLAGS = -lyaml-cpp -lm

SRC = source/main.cpp \
      source/LHE.cpp \
      source/Monte-Carlo.cpp \
      source/Physics.cpp 

BIN = bin/Eco+

all: $(BIN)

$(BIN): $(SRC)
	mkdir -p bin
	$(CXX) $(CXXFLAGS) $(SRC) -o $(BIN) $(LDFLAGS)

clean:
	rm -rf bin
