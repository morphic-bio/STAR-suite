# Makefile for hash screen replay tool and unit tests.
# Standalone — no STAR dependency.

CXX = g++
CXXFLAGS = -O3 -std=c++11 -Wall -Wextra

REPLAY_TARGET = hash_screen_replay
TEST_TARGET   = test_hash_screen

COMMON_HDRS = FlexHashScreenCommon.h FlexHashScreenFlat.h FlexHashScreenTiered.h
COMMON_SRCS = FlexHashScreenFlat.cpp FlexHashScreenTiered.cpp

.DEFAULT_GOAL := all

all: $(REPLAY_TARGET) $(TEST_TARGET) flex_hash_cache_pack

CORE_SRC = ../../../core/legacy/source
STORAGE_SRC = $(CORE_SRC)/FlexHashCacheStorage.cpp $(CORE_SRC)/FlexProbePairKhash.cpp
STORAGE_HDR = $(CORE_SRC)/FlexHashCacheStorage.h $(CORE_SRC)/FlexProbePairKhash.h $(CORE_SRC)/FlexHashScreenDecision.h ../../source/klib/khash.h

flex_hash_cache_pack: flex_hash_cache_pack.cpp $(STORAGE_SRC) $(STORAGE_HDR)
	$(CXX) $(CXXFLAGS) -I$(CORE_SRC) -I../../source -o $@ flex_hash_cache_pack.cpp $(STORAGE_SRC)

$(REPLAY_TARGET): hash_screen_replay.cpp $(COMMON_SRCS) $(COMMON_HDRS)
	$(CXX) $(CXXFLAGS) -o $@ hash_screen_replay.cpp $(COMMON_SRCS)

$(TEST_TARGET): test_hash_screen.cpp $(COMMON_SRCS) $(COMMON_HDRS)
	$(CXX) $(CXXFLAGS) -o $@ test_hash_screen.cpp $(COMMON_SRCS)

test: $(TEST_TARGET)
	./$(TEST_TARGET)

clean:
	rm -f $(REPLAY_TARGET) $(TEST_TARGET) flex_hash_cache_pack

.PHONY: all test clean
