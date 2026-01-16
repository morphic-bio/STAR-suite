CXX ?= g++
CXXFLAGS := -std=c++11 -O3 -Wall -Wextra

ROOT_DIR := $(abspath ../../../../..)
VBEM_SRC_DIR := $(ROOT_DIR)/core/features/vbem/source
LIBEM := $(VBEM_SRC_DIR)/libem/libem.a

all: compute_gc_bias

compute_gc_bias: compute_gc_bias.cpp $(LIBEM)
	$(CXX) $(CXXFLAGS) -I$(VBEM_SRC_DIR) compute_gc_bias.cpp $(LIBEM) -o compute_gc_bias

$(LIBEM):
	$(MAKE) -C $(VBEM_SRC_DIR)/libem

clean:
	rm -f compute_gc_bias

.PHONY: all clean
