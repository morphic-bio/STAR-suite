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
        src/OccupancyGuard.cpp \
        src/scrna_api.cpp

# Object files
OBJS := $(SRCS:.cpp=.o)

# Output library
TARGET := libscrna.a

# Public headers (for dependent projects)
PUBLIC_HEADERS := include/scrna_types.h \
                  include/SampleMatrixData.h \
                  include/EmptyDropsCRSampler.h \
                  include/EmptyDropsMultinomial.h \
                  include/OrdMagStage.h \
                  include/OccupancyGuard.h \
                  include/CRLogProb.h \
                  include/scrna_api.h

.PHONY: all clean lib test

all: lib

lib: $(TARGET)

$(TARGET): $(OBJS)
	$(AR) rcs $@ $^
	@echo "Built $(TARGET)"

src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
	rm -f tests/*.o tests/test_*

# Dependency on headers
src/EmptyDropsCRSampler.o: include/EmptyDropsCRSampler.h include/scrna_types.h include/pcg_random.hpp
src/EmptyDropsMultinomial.o: include/EmptyDropsMultinomial.h include/EmptyDropsCRSampler.h include/CRLogProb.h include/scrna_types.h include/SimpleGoodTuring/sgt.h
src/OrdMagStage.o: include/OrdMagStage.h include/scrna_types.h include/pcg_random.hpp
src/OccupancyGuard.o: include/OccupancyGuard.h include/scrna_types.h include/pcg_random.hpp

# Test targets (tests are in ../../tests/emptydrops/)
# Run: make -C ../../tests/emptydrops test
test: lib
	@echo "Tests moved to tests/emptydrops/"
	@echo "Run: make -C ../../tests/emptydrops test"
	$(MAKE) -C ../../tests/emptydrops test
