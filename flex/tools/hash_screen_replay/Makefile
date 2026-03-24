# Makefile for hash screen replay tool and unit tests.
# Standalone — no STAR dependency.

CXX = g++
CXXFLAGS = -O3 -std=c++11 -Wall -Wextra

REPLAY_TARGET = hash_screen_replay
TEST_TARGET   = test_hash_screen

COMMON_HDRS = FlexHashScreenCommon.h FlexHashScreenFlat.h FlexHashScreenTiered.h
COMMON_SRCS = FlexHashScreenFlat.cpp FlexHashScreenTiered.cpp

.DEFAULT_GOAL := all

all: $(REPLAY_TARGET) $(TEST_TARGET)

$(REPLAY_TARGET): hash_screen_replay.cpp $(COMMON_SRCS) $(COMMON_HDRS)
	$(CXX) $(CXXFLAGS) -o $@ hash_screen_replay.cpp $(COMMON_SRCS)

$(TEST_TARGET): test_hash_screen.cpp $(COMMON_SRCS) $(COMMON_HDRS)
	$(CXX) $(CXXFLAGS) -o $@ test_hash_screen.cpp $(COMMON_SRCS)

test: $(TEST_TARGET)
	./$(TEST_TARGET)

clean:
	rm -f $(REPLAY_TARGET) $(TEST_TARGET)

.PHONY: all test clean
