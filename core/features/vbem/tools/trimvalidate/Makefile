# Makefile for trimvalidate - CLI validator tool

CXX ?= g++
CXXFLAGS := -std=c++11 -O3 -Wall -Wextra

ROOT_DIR := $(abspath ../../../../..)
VBEM_SRC_DIR := $(ROOT_DIR)/core/features/vbem/source
LIBTRIM_DIR := $(VBEM_SRC_DIR)/libtrim
LIBTRIM := $(LIBTRIM_DIR)/libtrim.a

all: trimvalidate

trimvalidate: trimvalidate.cpp $(LIBTRIM)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LIBTRIM)

$(LIBTRIM):
	$(MAKE) -C $(LIBTRIM_DIR)

test: trimvalidate
	@./run_parity_test.sh

clean:
	rm -f trimvalidate

.PHONY: all clean test
