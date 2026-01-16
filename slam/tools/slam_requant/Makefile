CXX ?= g++
CXXFLAGS := -std=c++11 -O3 -fopenmp -Wall -Wextra
LDFLAGS := -fopenmp -lz

ROOT_DIR := $(abspath ../../../..)
SLAM_SRC_DIR := $(ROOT_DIR)/slam/source
LEGACY_SRC_DIR := $(ROOT_DIR)/core/legacy/source
SRC_DIR := $(SLAM_SRC_DIR)

OBJS := \
	$(SRC_DIR)/SlamQuant.cpp \
	$(SRC_DIR)/SlamSolver.cpp \
	$(SRC_DIR)/SlamCompat.cpp \
	$(SRC_DIR)/SlamVarianceAnalysis.cpp \
	$(SRC_DIR)/SlamReadBuffer.cpp \
	$(SRC_DIR)/SlamDump.cpp \
	$(SRC_DIR)/SlamQcOutput.cpp

HTSLIB := $(LEGACY_SRC_DIR)/htslib/libhts.a

all: slam_requant

slam_requant: slam_requant.cpp $(OBJS) $(HTSLIB)
	$(CXX) $(CXXFLAGS) -I$(SRC_DIR) -I$(LEGACY_SRC_DIR) $^ $(LDFLAGS) -o $@

$(HTSLIB):
	$(MAKE) -C $(LEGACY_SRC_DIR)/htslib lib-static

clean:
	rm -f slam_requant

.PHONY: all clean
