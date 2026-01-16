# Makefile for libflex - Flex pipeline shared library

# Compiler and flags (from parent Makefile, but don't include it to avoid Depend.list issues)
CXX ?= g++
# Keep -fopenmp for OccupancyGuard.cpp which still uses OpenMP
CXXFLAGS_common := -std=c++11 -fopenmp
CXXFLAGS_main := -O3 $(CXXFLAGS_common)
CXXFLAGS ?= -pipe -Wall -Wextra
CPPFLAGS ?=
CFLAGS ?= -pipe -Wall -Wextra -O3
AR ?= ar

# Object files for libflex (active flexfilter pipeline only)
LIBFLEX_OBJECTS = \
	OrdMagStage.o \
	EmptyDropsMultinomial.o \
	EmptyDropsCRSampler.o \
	OccupancyGuard.o \
	FlexFilter.o \
	FlexFilterAdapters.o \
	FlexFilterIO.o

# Build static library (default target)
all: libflex.a

libflex.a: $(LIBFLEX_OBJECTS)
	$(AR) rcs $@ $^
	ranlib $@

# Compile rules for libflex objects
%.o: %.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $(CXXFLAGS_common) -I.. -I../htslib -I. $< -o $@

%.o: %.c
	$(CXX) -c $(CPPFLAGS) $(CFLAGS) -I.. -I../htslib -I. $< -o $@

clean:
	rm -f $(LIBFLEX_OBJECTS) libflex.a

.PHONY: clean all
