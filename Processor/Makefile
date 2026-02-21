CXX = g++

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)

CXXFLAGS = -Wall -Wextra -I include $(ROOTCFLAGS)

BIN_DIR = bin

PLOT_SRC = source/main_plot.cpp \
      	   source/plot.C \

COMPARE_SRC = source/main_compare.cpp \
			  source/compare.C

PLOT_BIN = $(BIN_DIR)/Analyzer
COMPARE_BIN = $(BIN_DIR)/Comparator


all: $(PLOT_BIN) $(COMPARE_BIN)

$(PLOT_BIN) : $(PLOT_SRC)
	mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(PLOT_SRC) -o $(PLOT_BIN) $(ROOTLIBS)

$(COMPARE_BIN) : $(COMPARE_SRC)
	mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(COMPARE_SRC) -o $(COMPARE_BIN) $(ROOTLIBS)

clean:
	rm -rf $(BIN_DIR)

re: clean all
