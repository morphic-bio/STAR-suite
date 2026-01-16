# Makefile for em_quant CLI tool

CXX ?= g++
CXXFLAGS := -std=c++11 -O3 -march=native -fopenmp -Wall -Wextra
LDFLAGS := -fopenmp

LIBEM := ../../source/libem/libem.a

all: em_quant

em_quant: em_quant.cpp $(LIBEM)
	$(CXX) $(CXXFLAGS) -I../../source $< $(LIBEM) $(LDFLAGS) -o $@

$(LIBEM):
	$(MAKE) -C ../../source/libem

test: em_quant
	@./run_parity_test.sh

clean:
	rm -f em_quant

.PHONY: all clean test
