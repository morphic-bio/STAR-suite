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

# libscrna include path (for adapter headers in Stage 2+)
# From flex/source/libflex/ go up to STAR-suite root (../../..)
ROOT_DIR := $(abspath ../../..)
LIBSCRNA_INC := $(ROOT_DIR)/core/features/libscrna/include

# Object files for libflex (Flex-specific orchestration only)
# NOTE: OrdMagStage, EmptyDropsMultinomial, EmptyDropsCRSampler, OccupancyGuard
#       are now provided by libscrna (linked via core/legacy/source/Makefile).
#       Only Flex-specific orchestration code remains here.
LIBFLEX_OBJECTS = \
	FlexFilter.o \
	FlexFilterAdapters.o \
	FlexFilterIO.o

# Build static library (default target)
all: libflex.a

# NOTE: Delete archive first to avoid stale objects from previous builds.
# ar rcs does not remove old members, so without rm -f, migrated-away objects
# (OrdMagStage.o, EmptyDropsMultinomial.o, etc.) could persist and win the link.
libflex.a: $(LIBFLEX_OBJECTS)
	rm -f $@
	$(AR) rcs $@ $^
	ranlib $@

# Compile rules for libflex objects
%.o: %.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $(CXXFLAGS_common) -I.. -I../htslib -I. -I$(LIBSCRNA_INC) $< -o $@

%.o: %.c
	$(CXX) -c $(CPPFLAGS) $(CFLAGS) -I.. -I../htslib -I. $< -o $@

clean:
	rm -f $(LIBFLEX_OBJECTS) libflex.a

.PHONY: clean all
