CXX ?= g++
CXXFLAGS := -std=c++11 -O3 -fopenmp -Wall -Wextra
LDFLAGS := -fopenmp -lz

SRC_DIR := ../../source

OBJS := \
	$(SRC_DIR)/SlamQuant.cpp \
	$(SRC_DIR)/SlamSolver.cpp \
	$(SRC_DIR)/SlamCompat.cpp \
	$(SRC_DIR)/SlamVarianceAnalysis.cpp \
	$(SRC_DIR)/SlamReadBuffer.cpp \
	$(SRC_DIR)/SlamDump.cpp \
	$(SRC_DIR)/SlamQcOutput.cpp

HTSLIB := $(SRC_DIR)/htslib/libhts.a

all: slam_requant

slam_requant: slam_requant.cpp $(OBJS) $(HTSLIB)
	$(CXX) $(CXXFLAGS) -I$(SRC_DIR) $^ $(LDFLAGS) -o $@

$(HTSLIB):
	$(MAKE) -C $(SRC_DIR)/htslib

clean:
	rm -f slam_requant

.PHONY: all clean
