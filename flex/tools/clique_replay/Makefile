# Makefile for clique replay tool and unit tests
# Standalone — only depends on UMICorrector.{h,cpp} and klib/khash.h

ROOT_DIR := $(abspath ../../..)
FLEX_SRC := $(ROOT_DIR)/flex/source

CXX = g++
CXXFLAGS = -O3 -std=c++11 -Wall -Wextra -I$(FLEX_SRC)

UMI_CORRECTOR_SRC = $(FLEX_SRC)/UMICorrector.cpp

REPLAY_TARGET  = clique_replay
TEST_TARGET    = test_clique
FANOUT_TARGET  = test_fanout

.DEFAULT_GOAL := all

all: $(REPLAY_TARGET) $(TEST_TARGET) $(FANOUT_TARGET)

$(REPLAY_TARGET): clique_replay.cpp $(UMI_CORRECTOR_SRC) $(FLEX_SRC)/UMICorrector.h
	$(CXX) $(CXXFLAGS) -o $@ clique_replay.cpp $(UMI_CORRECTOR_SRC)

$(TEST_TARGET): test_clique.cpp $(UMI_CORRECTOR_SRC) $(FLEX_SRC)/UMICorrector.h
	$(CXX) $(CXXFLAGS) -o $@ test_clique.cpp $(UMI_CORRECTOR_SRC)

$(FANOUT_TARGET): test_fanout.cpp $(UMI_CORRECTOR_SRC) $(FLEX_SRC)/UMICorrector.h $(FLEX_SRC)/hash_shims_cpp_compat.h
	$(CXX) $(CXXFLAGS) -o $@ test_fanout.cpp $(UMI_CORRECTOR_SRC)

test: $(TEST_TARGET) $(FANOUT_TARGET)
	./$(TEST_TARGET)
	./$(FANOUT_TARGET)

clean:
	rm -f $(REPLAY_TARGET) $(TEST_TARGET) $(FANOUT_TARGET)

.PHONY: all test clean
