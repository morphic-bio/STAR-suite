CXX ?= g++
CXXFLAGS := -std=c++11 -O3 -Wall -Wextra

SCRIPT_DIR := $(shell dirname $(abspath $(lastword $(MAKEFILE_LIST))))
LIBEM := $(SCRIPT_DIR)/../../source/libem/libem.a

all: compute_gc_bias

compute_gc_bias: compute_gc_bias.cpp $(LIBEM)
	$(CXX) $(CXXFLAGS) -I$(SCRIPT_DIR)/../../source compute_gc_bias.cpp $(LIBEM) -o compute_gc_bias

$(LIBEM):
	$(MAKE) -C $(SCRIPT_DIR)/../../source/libem

clean:
	rm -f compute_gc_bias

.PHONY: all clean
