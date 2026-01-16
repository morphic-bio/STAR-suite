# Makefile for em_quant CLI tool

CXX ?= g++
CXXFLAGS := -std=c++11 -O3 -march=native -fopenmp -Wall -Wextra
LDFLAGS := -fopenmp

ROOT_DIR := $(abspath ../../../../..)
VBEM_SRC_DIR := $(ROOT_DIR)/core/features/vbem/source
LIBEM := $(VBEM_SRC_DIR)/libem/libem.a

all: em_quant

em_quant: em_quant.cpp $(LIBEM)
	$(CXX) $(CXXFLAGS) -I$(VBEM_SRC_DIR) $< $(LIBEM) $(LDFLAGS) -o $@

$(LIBEM):
	$(MAKE) -C $(VBEM_SRC_DIR)/libem

test: em_quant
	@./run_parity_test.sh

clean:
	rm -f em_quant

.PHONY: all clean test
