# Makefile for run_flexfilter_mex standalone tool
# Moved from reference/tools/Makefile.flexfilter

ROOT_DIR := $(abspath ../../..)
LEGACY_SRC_DIR := $(ROOT_DIR)/core/legacy/source
LIBFLEX_DIR := $(ROOT_DIR)/flex/source/libflex
LIBFLEX := $(LIBFLEX_DIR)/libflex.a
LIBSCRNA_DIR := $(ROOT_DIR)/core/features/libscrna
LIBSCRNA := $(LIBSCRNA_DIR)/libscrna.a
LIBSCRNA_INC := $(LIBSCRNA_DIR)/include
HTSLIB_DIR := $(LEGACY_SRC_DIR)/htslib
HTSLIB := $(HTSLIB_DIR)/libhts.a

CXX = g++
# Keep -fopenmp for OccupancyGuard.cpp which still uses OpenMP (now in libscrna)
CXXFLAGS = -O3 -std=c++11 -Wall -Wextra -fopenmp \
           -I$(ROOT_DIR)/flex/source -I$(LIBFLEX_DIR) \
           -I$(LIBSCRNA_INC) \
           -I$(LEGACY_SRC_DIR) -I$(HTSLIB_DIR)

# Link deps: libflex.a, libscrna.a (provides EmptyDrops/OrdMag/Occupancy), MexWriter.o, htslib
# NOTE: libscrna.a must come after libflex.a since FlexFilter depends on libscrna symbols
LDFLAGS = $(LIBFLEX) $(LIBSCRNA) $(LEGACY_SRC_DIR)/MexWriter.o \
          -L$(HTSLIB_DIR) -lhts \
          -pthread -fopenmp

TARGET = run_flexfilter_mex
SRC = run_flexfilter_mex.cpp
OBJ = $(SRC:.cpp=.o)

# Always rebuild libflex to avoid stale linkage when headers change
.PHONY: libflex_rebuild libscrna_rebuild
libflex_rebuild:
	$(MAKE) -C $(LIBFLEX_DIR) clean
	$(MAKE) -C $(LIBFLEX_DIR) libflex.a

libscrna_rebuild:
	$(MAKE) -C $(LIBSCRNA_DIR) libscrna.a

all: $(TARGET)

$(TARGET): libflex_rebuild libscrna_rebuild $(OBJ) $(LIBFLEX) $(LIBSCRNA) $(HTSLIB)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJ) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)

$(HTSLIB):
	$(MAKE) -C $(HTSLIB_DIR) lib-static

.PHONY: all clean
