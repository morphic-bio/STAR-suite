# Makefile for run_flexfilter_mex standalone tool
# Moved from reference/tools/Makefile.flexfilter

LIBFLEX_DIR = ../../source/libflex
LIBFLEX = $(LIBFLEX_DIR)/libflex.a
HTSLIB_DIR = ../../source/htslib

CXX = g++
# Keep -fopenmp for OccupancyGuard.cpp which still uses OpenMP
# Use pkg-config for glib-2.0 instead of hard-coded paths
GLIB_CFLAGS = $(shell pkg-config --cflags glib-2.0)
GLIB_LIBS = $(shell pkg-config --libs glib-2.0)

CXXFLAGS = -O3 -std=c++11 -Wall -Wextra -fopenmp \
           -I../../source -I../../source/libflex \
           -I../../reference/tools/flex_debug/flex/third_party \
           -I$(HTSLIB_DIR) \
           $(GLIB_CFLAGS)

# Link deps: libflex.a, MexWriter.o, htslib, glib, pthread, fopenmp
LDFLAGS = $(LIBFLEX) ../../source/MexWriter.o \
          -L$(HTSLIB_DIR) -lhts \
          $(GLIB_LIBS) -pthread -fopenmp

TARGET = run_flexfilter_mex
SRC = run_flexfilter_mex.cpp
OBJ = $(SRC:.cpp=.o)

# Always rebuild libflex to avoid stale linkage when headers change
.PHONY: libflex_rebuild
libflex_rebuild:
	$(MAKE) -C $(LIBFLEX_DIR) clean
	$(MAKE) -C $(LIBFLEX_DIR) libflex.a

all: $(TARGET)

$(TARGET): libflex_rebuild $(OBJ) $(LIBFLEX)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJ) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)

.PHONY: all clean
