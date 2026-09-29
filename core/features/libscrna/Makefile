# libscrna - Shared single-cell RNA analysis library
# Contains: EmptyDrops, OrdMag (Simple EmptyDrops), Occupancy filtering

CXX := g++
AR := ar
CXXFLAGS := -std=c++11 -O3 -Wall -Wextra -fPIC
INCLUDES := -I./include

# Source files
SRCS := src/EmptyDropsCRSampler.cpp \
        src/EmptyDropsMultinomial.cpp \
        src/OrdMagStage.cpp \
        src/SimpleEDCaller.cpp \
        src/OccupancyGuard.cpp \
        src/scrna_api.cpp

# Object files
OBJS := $(SRCS:.cpp=.o)

# Output library
TARGET := libscrna.a

# Public headers (for dependent projects)
PUBLIC_HEADERS := include/scrna_types.h \
                  include/SampleMatrixData.h \
                  include/AdaptiveAmbientWindow.h \
                  include/FlexTagGroup.h \
                  include/EmptyDropsCRSampler.h \
                  include/EmptyDropsMultinomial.h \
                  include/OrdMagStage.h \
                  include/ParallelTasks.h \
                  include/ThreadPermits.h \
                  include/OrdMagRank.h \
                  include/SimpleEDCaller.h \
                  include/MitochondrialRankMask.h \
                  include/OccupancyGuard.h \
                  include/CRLogProb.h \
                  include/scrna_api.h

.PHONY: all clean lib test tools

all: lib

lib: $(TARGET)

# Simple EmptyDrops CLI (wrapper around libscrna C API)
TOOLS_DIR := tools
BIN_DIR := bin
SIMPLEED_TOOL := $(BIN_DIR)/scrna_simpleed

tools: lib $(SIMPLEED_TOOL)

$(SIMPLEED_TOOL): $(TOOLS_DIR)/scrna_simpleed.cpp $(TARGET)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(TARGET)
	@echo "Built $(SIMPLEED_TOOL)"

$(TARGET): $(OBJS)
	rm -f $@
	$(AR) rcs $@ $^
	@echo "Built $(TARGET)"

src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
	rm -f tests/*.o tests/test_*
	rm -f $(SIMPLEED_TOOL)

# Dependency on headers
src/EmptyDropsCRSampler.o src/OrdMagStage.o src/OccupancyGuard.o: include/pcg_extras.hpp include/pcg_uint128.hpp
src/EmptyDropsCRSampler.o: include/EmptyDropsCRSampler.h include/scrna_types.h include/pcg_random.hpp include/ParallelTasks.h include/ThreadPermits.h
src/EmptyDropsMultinomial.o: include/SparseCountView.h include/EmptyDropsMultinomial.h include/EmptyDropsCRSampler.h include/CRLogProb.h include/scrna_types.h include/SimpleGoodTuring/sgt.h
src/OrdMagStage.o: include/OrdMagStage.h include/OrdMagRank.h include/ParallelTasks.h include/ThreadPermits.h include/scrna_types.h include/pcg_random.hpp
src/scrna_api.o: include/ScrnaTrace.h include/SparseCountView.h include/scrna_api.h include/OrdMagStage.h include/OrdMagRank.h include/EmptyDropsMultinomial.h include/scrna_types.h
src/OccupancyGuard.o: include/OccupancyGuard.h include/scrna_types.h include/pcg_random.hpp

# Test targets (tests are in ../../../tests/emptydrops/)
# Run: make -C ../../../tests/emptydrops test
test: lib
	@echo "Tests moved to tests/emptydrops/"
	@echo "Run: make -C ../../../tests/emptydrops test"
	$(MAKE) -C ../../../tests/emptydrops test

src/SimpleEDCaller.o: include/SparseCountView.h include/SimpleEDCaller.h include/OrdMagStage.h include/OrdMagRank.h include/AdaptiveAmbientWindow.h include/EmptyDropsMultinomial.h include/scrna_api.h
