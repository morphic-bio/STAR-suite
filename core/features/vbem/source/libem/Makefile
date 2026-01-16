# Makefile for libem - EM/VB quantification library

CXX ?= g++
CXXFLAGS_common := -std=c++11 -O3 -march=native -fopenmp -Wall -Wextra
CXXFLAGS ?= -pipe
AR ?= ar

# OpenBLAS by default (not used in v1, but ready for future)
LDFLAGS_BLAS ?= -lopenblas

LIBEM_OBJECTS = ec_loader.o em_engine.o vb_engine.o alignment_filter.o ec_builder.o extended_pruning.o alignment_model.o gc_bias.o effective_length.o effective_length_wrapper.o fld_accumulator.o slam_snp_em.o

all: libem.a

libem.a: $(LIBEM_OBJECTS)
	$(AR) rcs $@ $^
	ranlib $@

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $(CXXFLAGS_common) -I.. $< -o $@

clean:
	rm -f $(LIBEM_OBJECTS) libem.a

.PHONY: all clean
